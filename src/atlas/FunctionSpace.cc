/*
 * (C) Copyright 1996-2015 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include <cassert>
#include <iostream>
#include <sstream>
#include <limits>

#include "eckit/exception/Exceptions.h"
#include "eckit/types/Types.h"

#include "atlas/atlas_config.h"
#include "atlas/runtime/ErrorHandling.h"
#include "atlas/FunctionSpace.h"
#include "atlas/Field.h"
#include "atlas/actions/BuildParallelFields.h"
#include "atlas/util/Debug.h"
#include "atlas/util/Bitflags.h"
#include "atlas/util/DataType.h"

namespace atlas {
namespace next {

namespace {

void assert_shared(const eckit::Owned* owned)
{
  if( owned->owners() == 0 )
  {
    throw eckit::SeriousBug("Cannot create shared_ptr from stack allocated or scoped_ptr",Here());
  }
}

}

eckit::SharedPtr<FunctionSpace const> FunctionSpace::shared_from_this() const
{
  assert_shared(this);
  return eckit::SharedPtr<FunctionSpace const>(this);
}

eckit::SharedPtr<FunctionSpace> FunctionSpace::shared_from_this()
{
  assert_shared(this);
  return eckit::SharedPtr<FunctionSpace>(this);
}

eckit::SharedPtr<FunctionSpace const> FunctionSpace::ptr() const
{
  assert_shared(this);
  return eckit::SharedPtr<FunctionSpace const>(this);
}

eckit::SharedPtr<FunctionSpace const> FunctionSpace::cptr() const
{
  ASSERT(owners()!=0);
  return eckit::SharedPtr<FunctionSpace const>(this);
}

eckit::SharedPtr<FunctionSpace> FunctionSpace::ptr()
{
  ASSERT(owners()!=0);
  return eckit::SharedPtr<FunctionSpace>(this);
}

}
}


using atlas::util::Topology;

#ifdef ATLAS_HAVE_FORTRAN
#define REMOTE_IDX_BASE 1
#else
#define REMOTE_IDX_BASE 0
#endif

namespace atlas {

template<class T>
inline std::ostream& operator<<(std::ostream& s,const std::vector<T>& v)
{
	return eckit::__print_list(s,v);
}

//----------------------------------------------------------------------------------------------------------------------

#if !DEPRECATE_OLD_FUNCTIONSPACE
FunctionSpace::FunctionSpace(const std::string& name,
                             const std::string& shape_func,
                             const std::vector<size_t>& shape,
                             Mesh& mesh) :
  name_(name),
  shape_(shape),
  gather_scatter_(new mpl::GatherScatter()),
  fullgather_(new mpl::GatherScatter()),
  halo_exchange_(new mpl::HaloExchange()),
  checksum_(new mpl::Checksum()),
  mesh_(&mesh)
{
  //std::cout << "C++ : shape Constructor" << std::endl;
  dof_ = 1;
  size_t extsize = shape_.size();
  shapef_.resize(extsize);
  for (size_t i=0; i<extsize; ++i)
  {
    shapef_[extsize-1-i] = shape_[i];
    if( shape_[i] != FunctionSpace::UNDEF_VARS )
      dof_ *= shape_[i];
  }
  glb_dof_ = dof_;
}

/// TEMPORARY CONSTRUCTOR, JUST FOR EVOLUTIONARY STEP TO NEXT DESIGN
FunctionSpace::FunctionSpace(const std::string& name,
                             const std::vector<size_t>& shape) :
  name_(name),
  shape_(shape),
  gather_scatter_(new mpl::GatherScatter()),
  fullgather_(new mpl::GatherScatter()),
  halo_exchange_(new mpl::HaloExchange()),
  checksum_(new mpl::Checksum()),
  mesh_(0)
{
  //std::cout << "C++ : shape Constructor" << std::endl;
  dof_ = 1;
  size_t extsize = shape_.size();
  shapef_.resize(extsize);
  for (size_t i=0; i<extsize; ++i)
  {
    shapef_[extsize-1-i] = shape_[i];
    if( shape_[i] != FunctionSpace::UNDEF_VARS )
      dof_ *= shape_[i];
  }
  glb_dof_ = dof_;
}


FunctionSpace::~FunctionSpace()
{
}

void FunctionSpace::resize(const std::vector<size_t>& shape)
{
  if (shape.size() != shape_.size() )
    throw eckit::BadParameter("Cannot resize shape: shape sizes don't match.",Here());

  size_t extsize = shape_.size();

  for (size_t i=1; i<extsize; ++i)
  {
    if (shape[i] != shape_[i])
      throw eckit::BadParameter("Only the first extent can be resized for now!",Here());
  }

  shape_ = shape;
  shapef_.resize(extsize);
  for (size_t i=0; i<extsize; ++i)
  {
    shapef_[extsize-1-i] = shape_[i];
  }

  dof_ = 1;
  for (size_t i=0; i<extsize; ++i)
  {
    if( shape_[i] != FunctionSpace::UNDEF_VARS )
      dof_ *= shape_[i];
  }

  for( size_t f=0; f<nb_fields(); ++f)
  {
    std::vector< size_t > field_shape(extsize);
    for (size_t i=0; i<extsize; ++i)
    {
      if( shape_[i] == FunctionSpace::UNDEF_VARS )
        field_shape[i] = field(f).shape(i);
      else
        field_shape[i] = shape_[i];
    }
    field(f).resize(field_shape);
  }
}

namespace {


  Field* check_if_exists( const FunctionSpace& fs,
                const std::string& name,
                const std::vector<size_t>&  shape,
                size_t nb_vars,
                CreateBehavior b )
  {
    using namespace eckit;

    if( fs.has_field(name) )
    {
      if( b == IF_EXISTS_FAIL )
      {
        std::ostringstream msg; msg << "field with name " << name << " already exists" << std::endl;
        throw eckit::Exception( msg.str(), Here() );
      }

      const Field& f= fs.field(name);

      if( f.shape() != shape )
      {
        std::ostringstream msg; msg << "field exists with name " << name << " has unexpected shape ";
        __print_list(msg, f.shape());
        msg << " instead of ";
        __print_list(msg, shape);
        msg << std::endl;
        throw eckit::Exception(msg.str(),Here());
      }
      return const_cast<Field*>(&f);
    }

    return 0;
  }


  template <typename T>
  Field& on_create_field(FunctionSpace& fs, const std::string& name, size_t nb_vars, CreateBehavior b )
  {
    Field* field = 0;

    size_t rank = fs.shape().size();
    std::vector< size_t > field_shape(rank);
    for (size_t i=0; i<rank; ++i)
    {
      if( fs.shape()[i] == FunctionSpace::UNDEF_VARS )
        field_shape[i] = nb_vars;
      else
        field_shape[i] = fs.shape()[i];
    }

    if( (field = check_if_exists(fs, name, field_shape, nb_vars, b )) )
      return *field;

    field = Field::create<T>(name,field_shape);
    fs.add(field);

    return *field;
  }


}

Field& FunctionSpace::add( Field* field )
{
  fields_.insert( field->name(), Field::Ptr(field) );
  fields_.sort();
  return *field;
}

template<> Field& FunctionSpace::create_field<double>(const std::string& name, size_t nb_vars, CreateBehavior b )
{
  return on_create_field<double>(*this,name,nb_vars,b);
}

template<> Field& FunctionSpace::create_field<float>(const std::string& name, size_t nb_vars, CreateBehavior b )
{
  return on_create_field<float>(*this,name,nb_vars,b);
}

template<> Field& FunctionSpace::create_field<int>(const std::string& name, size_t nb_vars, CreateBehavior b )
{
  return on_create_field<int>(*this,name,nb_vars,b);
}

template<> Field& FunctionSpace::create_field<long>(const std::string& name, size_t nb_vars, CreateBehavior b )
{
  return on_create_field<long>(*this,name,nb_vars,b);
}



void FunctionSpace::remove_field(const std::string& name)
{
  NOTIMP; ///< @todo DenseMap needs to have erase() function

//  if( has_field(name) )
//  {
//    fields_.erase(name);
//  }
//  else
//  {
//    std::stringstream msg;
//    msg << "Could not find field \"" << name << "\" in FunctionSpace \"" << name_ << "\"";
//    throw eckit::OutOfRange(msg.str(),Here());
//  }
}

const Field& FunctionSpace::field( size_t idx ) const
{
  return *fields_.at( idx );
}

Field& FunctionSpace::field( size_t idx )
{
  return *fields_.at( idx );
}


const Field& FunctionSpace::field(const std::string& name) const
{
  if( has_field(name) )
  {
    return *fields_.get(name);
  }
  else
  {
    std::stringstream msg;
    msg << "Could not find field \"" << name << "\" in FunctionSpace \"" << name_ << "\"";
    throw eckit::OutOfRange(msg.str(),Here());
  }
}

Field& FunctionSpace::field(const std::string& name)
{
  if( has_field(name) )
  {
    return *fields_.get(name);
  }
  else
  {
    std::stringstream msg;
    msg << "Could not find field \"" << name << "\" in FunctionSpace \"" << name_ << "\"";
    throw eckit::OutOfRange(msg.str(),Here());
  }
}

void FunctionSpace::parallelise(const int part[], const int remote_idx[], const gidx_t glb_idx[], size_t parsize)
{
  halo_exchange_->setup(part,remote_idx,REMOTE_IDX_BASE,parsize);
  gather_scatter_->setup(part,remote_idx,REMOTE_IDX_BASE,glb_idx,-1,parsize);
  fullgather_->setup(part,remote_idx,REMOTE_IDX_BASE,glb_idx,-1,parsize,true);
  checksum_->setup(part,remote_idx,REMOTE_IDX_BASE,glb_idx,-1,parsize);
  glb_dof_ = gather_scatter_->glb_dof();
  for( int b=shapef_.size()-2; b>=0; --b)
  {
    if( shapef_[b] != FunctionSpace::UNDEF_VARS )
      glb_dof_ *= shapef_[b];
  }
}

void FunctionSpace::parallelise(FunctionSpace& other_shape)
{
  halo_exchange_ = mpl::HaloExchange::Ptr( &other_shape.halo_exchange() );
    gather_scatter_ = mpl::GatherScatter::Ptr( &other_shape.gather_scatter() );
}

void FunctionSpace::parallelise()
{
  Field& ridx = field("remote_idx");
  Field& part = field("partition");
  Field& gidx = field("glb_idx");

  if( name() == "nodes")
  {
    ArrayView<int,1> flags ( field("flags") );
    std::vector<int> mask(shape(0));
    for( size_t j=0; j<mask.size(); ++j )
    {
      mask[j] = Topology::check(flags(j),Topology::GHOST) ? 1 : 0;
    }
    halo_exchange_->setup(part.data<int>(),ridx.data<int>(),REMOTE_IDX_BASE,shape(0));
    gather_scatter_->setup(part.data<int>(),ridx.data<int>(),REMOTE_IDX_BASE,gidx.data<gidx_t>(),mask.data(),shape(0));
    fullgather_->setup(part.data<int>(),ridx.data<int>(),REMOTE_IDX_BASE,gidx.data<gidx_t>(),-1,shape(0),true);
    checksum_->setup(part.data<int>(),ridx.data<int>(),REMOTE_IDX_BASE,gidx.data<gidx_t>(),mask.data(),shape(0));
    glb_dof_ = gather_scatter_->glb_dof();
    for( int b=shapef_.size()-2; b>=0; --b)
    {
      if( shapef_[b] != FunctionSpace::UNDEF_VARS )
        glb_dof_ *= shapef_[b];
    }
  }
  else
  {
    halo_exchange_->setup(part.data<int>(),ridx.data<int>(),REMOTE_IDX_BASE,shape(0));
    gather_scatter_->setup(part.data<int>(),ridx.data<int>(),REMOTE_IDX_BASE,gidx.data<gidx_t>(),-1,shape(0));
    fullgather_->setup(part.data<int>(),ridx.data<int>(),REMOTE_IDX_BASE,gidx.data<gidx_t>(),-1,shape(0),true);
    checksum_->setup(part.data<int>(),ridx.data<int>(),REMOTE_IDX_BASE,gidx.data<gidx_t>(),-1,shape(0));
  }
  glb_dof_ = gather_scatter_->glb_dof();
  for( int b=shapef_.size()-2; b>=0; --b)
  {
    if( shapef_[b] != FunctionSpace::UNDEF_VARS )
      glb_dof_ *= shapef_[b];
  }
}

void FunctionSpace::print(std::ostream& os, bool dump) const
{
    os << "FunctionSpace[name=" << name() << ",";
    for(size_t i = 0; i < nb_fields(); ++i)
    {
        if(dump)
            field(i).dump(os);
        else
            os << field(i) << ",";
    }
    os << "]";
}
#endif

//----------------------------------------------------------------------------------------------------------------------

// C wrapper interfaces to C++ routines

void atlas__NextFunctionSpace__delete (next::FunctionSpace* This)
{
  ATLAS_ERROR_HANDLING(
    ASSERT( This );
    delete This;
    This = 0;
  );
}

const char* atlas__NextFunctionSpace__name (next::FunctionSpace* This) {
  ATLAS_ERROR_HANDLING(
    ASSERT( This );
    return This->name().c_str();
  );
  return 0;
}

#if !DEPRECATE_OLD_FUNCTIONSPACE

Metadata* atlas__FunctionSpace__metadata (FunctionSpace* This)
{
  ASSERT( This );
  return &This->metadata();
}

int atlas__FunctionSpace__dof (FunctionSpace* This) {
  ASSERT( This );
  return This->dof();
}

int atlas__FunctionSpace__glb_dof (FunctionSpace* This) {
  ASSERT( This );
  return This->glb_dof();
}

void atlas__FunctionSpace__create_field_double (FunctionSpace* This, char* name, int nb_vars) {
  ATLAS_ERROR_HANDLING(
      ASSERT( This );
      This->create_field<double>( std::string(name), nb_vars )
  );
}

void atlas__FunctionSpace__create_field_float (FunctionSpace* This, char* name, int nb_vars) {
  ASSERT( This );
  ATLAS_ERROR_HANDLING( This->create_field<float>( std::string(name), nb_vars ) );
}

void atlas__FunctionSpace__create_field_int (FunctionSpace* This, char* name, int nb_vars) {
  ASSERT( This );
  ATLAS_ERROR_HANDLING( This->create_field<int>( std::string(name), nb_vars ) );
}

void atlas__FunctionSpace__create_field_long (FunctionSpace* This, char* name, int nb_vars) {
  ASSERT( This );
  ATLAS_ERROR_HANDLING( This->create_field<long>( std::string(name), nb_vars ) );
}

void atlas__FunctionSpace__remove_field (FunctionSpace* This, char* name ) {
  ASSERT( This );
  ATLAS_ERROR_HANDLING( This->remove_field( std::string(name) ) );
}

int atlas__FunctionSpace__has_field (FunctionSpace* This, char* name) {
  ASSERT( This );
  return This->has_field( std::string(name) );
}

const char* atlas__FunctionSpace__name (FunctionSpace* This) {
  ASSERT( This );
  return This->name().c_str();
}

void atlas__FunctionSpace__shapef (FunctionSpace* This, int* &shape, int &rank) {
  ASSERT( This );
  shape = const_cast<int*>(&(This->shapef()[0]));
  rank = This->shapef().size();
}

Field* atlas__FunctionSpace__field (FunctionSpace* This, char* name) {
  ASSERT( This );
  ATLAS_ERROR_HANDLING( return &This->field( std::string(name) ) );
  return 0;
}

void atlas__FunctionSpace__parallelise (FunctionSpace* This) {
  ASSERT( This );
  ATLAS_ERROR_HANDLING( This->parallelise() );
}

void atlas__FunctionSpace__halo_exchange_int (FunctionSpace* This, int field_data[], int field_size) {
  ASSERT( This );
  ATLAS_ERROR_HANDLING(
    This->halo_exchange(field_data,field_size) );
}

void atlas__FunctionSpace__halo_exchange_float (FunctionSpace* This, float field_data[], int field_size) {
  ASSERT( This );
  ATLAS_ERROR_HANDLING( This->halo_exchange(field_data,field_size) );
}

void atlas__FunctionSpace__halo_exchange_double (FunctionSpace* This, double field_data[], int field_size) {
  ASSERT( This );
  ATLAS_ERROR_HANDLING( This->halo_exchange(field_data,field_size) );
}

void atlas__FunctionSpace__gather_int (FunctionSpace* This, int field_data[], int field_size, int glbfield_data[], int glbfield_size) {
  ASSERT( This );
  ATLAS_ERROR_HANDLING(
    This->gather(field_data,field_size, glbfield_data,glbfield_size) );
}

void atlas__FunctionSpace__gather_float (FunctionSpace* This, float field_data[], int field_size, float glbfield_data[], int glbfield_size) {
  ASSERT( This );
  ATLAS_ERROR_HANDLING(
    This->gather(field_data,field_size, glbfield_data,glbfield_size) );
}

void atlas__FunctionSpace__gather_double (FunctionSpace* This, double field_data[], int field_size, double glbfield_data[], int glbfield_size) {
  ASSERT( This );
  ATLAS_ERROR_HANDLING(
    This->gather(field_data,field_size, glbfield_data,glbfield_size) );
}

mpl::HaloExchange* atlas__FunctionSpace__halo_exchange (FunctionSpace* This) {
  ASSERT( This );
  ATLAS_ERROR_HANDLING( return &This->halo_exchange() );
  return 0;
}

mpl::GatherScatter* atlas__FunctionSpace__gather (FunctionSpace* This) {
  ASSERT( This );
  ATLAS_ERROR_HANDLING( return &This->gather_scatter() );
  return 0;
}

mpl::Checksum* atlas__FunctionSpace__checksum (FunctionSpace* This) {
  ASSERT( This );
  ATLAS_ERROR_HANDLING( return &This->checksum() );
  return 0;
}

void atlas__FunctionSpace__delete (FunctionSpace* This) {
  ASSERT( This );
  delete This;
}
#endif
// ------------------------------------------------------------------

} // namespace atlas

