/*
 * (C) Copyright 1996-2017 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "atlas/functionspace/StructuredColumns.h"

#include "eckit/utils/MD5.h"

#include "atlas/mesh/Mesh.h"
#include "atlas/field/FieldSet.h"
#include "atlas/util/Checksum.h"
#include "atlas/runtime/ErrorHandling.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/parallel/GatherScatter.h"
#include "atlas/array/MakeView.h"

namespace atlas {
namespace functionspace {
namespace detail {

namespace {
void set_field_metadata(const eckit::Parametrisation& config, Field& field)
{
  bool global(false);
  if( config.get("global",global) )
  {
    if( global )
    {
      size_t owner(0);
      config.get("owner",owner);
      field.metadata().set("owner",owner);
    }
  }
  field.metadata().set("global",global);
}
}

size_t StructuredColumns::config_size(const eckit::Parametrisation& config) const
{
  size_t size = this->size();
  bool global(false);
  if( config.get("global",global) )
  {
    if( global )
    {
      size_t owner(0);
      config.get("owner",owner);
      size = (parallel::mpi::comm().rank() == owner ? grid_.size() : 0);
    }
  }
  return size;
}


// ----------------------------------------------------------------------------
// Constructor
// ----------------------------------------------------------------------------
StructuredColumns::StructuredColumns(const Grid& grid) :
  StructuredColumns::StructuredColumns(grid, grid::Partitioner() ) {
}

StructuredColumns::StructuredColumns( const Grid& grid, const grid::Partitioner& p ) :
  grid_(grid)
{
    if ( not grid_ )
    {
      throw eckit::BadCast("Grid is not a grid::Structured type", Here());
    }

    grid::Partitioner partitioner( p );
    if( not partitioner ) {
      if( grid_.domain().global() ) {
        if( grid::Partitioner::exists("trans") )
          partitioner = grid::Partitioner("trans");
        else
          partitioner = grid::Partitioner("equal_regions");
      } else {
        partitioner = grid::Partitioner("checkerboard");
      }
    }

    grid::Distribution distribution(grid,partitioner);

    int mpi_rank = parallel::mpi::comm().rank();
    size_t ny( 0 );
    size_t j_begin( std::numeric_limits<size_t>::max() ) ;
    std::vector<size_t> i_begin( grid_.ny(), std::numeric_limits<size_t>::max() );
    std::vector<size_t> nx( grid_.ny(), 0 );
    size_t c( 0 );
    for( size_t j=0; j<grid_.ny(); ++j ) {
      bool j_used(false);
      for( size_t i=0; i<grid_.nx(j); ++i ) {
         if( distribution.partition(c++) == mpi_rank ) {
           j_used = true;
           j_begin = std::min( j_begin, j );
           i_begin[j-j_begin] = std::min( i_begin[j-j_begin], i );
           ++nx[j-j_begin];
         }
      }
      if( j_used )
        ++ny;
    }

    // export
    ny_ = ny;
    j_begin_ = j_begin;
    i_begin_.assign(i_begin.begin(),i_begin.begin()+ny);
    nx_.assign(nx.begin(),nx.begin()+ny);
    npts_ = 0;
    for( size_t& n : nx_ ) npts_ += n;

    gather_scatter_ = new parallel::GatherScatter();

    std::vector<int>    part(npts_,mpi_rank);
    std::vector<int>    remote_idx; remote_idx.reserve(npts_);
    std::vector<gidx_t> global_idx; global_idx.reserve(npts_);

    size_t gidx_begin=1;
    for( size_t j=0; j<j_begin_; ++j ) {
      gidx_begin += grid_.nx(j);
    }
    c = 0;
    for( size_t j=0; j<ny_; ++j ) {
      for( size_t i=0; i<nx_[j]; ++i ) {
         global_idx.push_back( gidx_begin + i_begin[j] + i );
         remote_idx.push_back( c++ );
      }
      gidx_begin += grid_.nx(j_begin_+j);
    }
    gather_scatter_->setup(part.data(), remote_idx.data(), 0, global_idx.data(), npts_ );
}
// ----------------------------------------------------------------------------



// ----------------------------------------------------------------------------
// Destructor
// ----------------------------------------------------------------------------
StructuredColumns::~StructuredColumns()
{
}
// ----------------------------------------------------------------------------


size_t StructuredColumns::footprint() const {
  size_t size = sizeof(*this);
  // TODO
  return size;
}

// ----------------------------------------------------------------------------
// Create Field
// ----------------------------------------------------------------------------
Field StructuredColumns::createField(const std::string& name, array::DataType datatype, const eckit::Parametrisation& options ) const
{
    size_t npts = config_size(options);
    Field field = Field(name, datatype, array::make_shape(npts));
    field.set_functionspace(this);
    set_field_metadata(options,field);
    return field;
}
// ----------------------------------------------------------------------------



// ----------------------------------------------------------------------------
// Create Field with vertical levels
// ----------------------------------------------------------------------------
Field StructuredColumns::createField(
    const std::string& name, array::DataType datatype,
    size_t levels, const eckit::Parametrisation& options) const
{
    size_t npts = config_size(options);
    Field field = Field(
                    name, datatype, array::make_shape(npts, levels));

    field.set_functionspace(this);
    field.set_levels(levels);
    set_field_metadata(options,field);
    return field;
}
// ----------------------------------------------------------------------------


// ----------------------------------------------------------------------------
// Gather FieldSet
// ----------------------------------------------------------------------------
void StructuredColumns::gather(
    const FieldSet& local_fieldset,
    FieldSet& global_fieldset ) const
{
    ASSERT(local_fieldset.size() == global_fieldset.size());

    for( size_t f=0; f<local_fieldset.size(); ++f ) {

      const Field& loc = local_fieldset[f];
      Field& glb = global_fieldset[f];
      const size_t nb_fields = 1;
      size_t root(0);
      glb.metadata().get("owner",root);

      if     ( loc.datatype() == array::DataType::kind<int>() ) {
        parallel::Field<int const> loc_field( array::make_storageview<int>(loc).data(),loc.stride(0));
        parallel::Field<int      > glb_field( array::make_storageview<int>(glb).data(),glb.stride(0));
        gather_scatter_->gather( &loc_field, &glb_field, nb_fields, root );
      }
      else if( loc.datatype() == array::DataType::kind<long>() ) {
        parallel::Field<long const> loc_field( array::make_storageview<long>(loc).data(),loc.stride(0));
        parallel::Field<long      > glb_field( array::make_storageview<long>(glb).data(),glb.stride(0));
        gather_scatter_->gather( &loc_field, &glb_field, nb_fields, root );
      }
      else if( loc.datatype() == array::DataType::kind<float>() ) {
        parallel::Field<float const> loc_field( array::make_storageview<float>(loc).data(),loc.stride(0));
        parallel::Field<float      > glb_field( array::make_storageview<float>(glb).data(),glb.stride(0));
        gather_scatter_->gather( &loc_field, &glb_field, nb_fields, root );
      }
      else if( loc.datatype() == array::DataType::kind<double>() ) {
        parallel::Field<double const> loc_field( array::make_storageview<double>(loc).data(),loc.stride(0));
        parallel::Field<double      > glb_field( array::make_storageview<double>(glb).data(),glb.stride(0));
        gather_scatter_->gather( &loc_field, &glb_field, nb_fields, root );
      }
      else throw eckit::Exception("datatype not supported",Here());
    }
}
// ----------------------------------------------------------------------------



// ----------------------------------------------------------------------------
// Gather Field
// ----------------------------------------------------------------------------
void StructuredColumns::gather(
    const Field& local,
    Field& global) const
{
    FieldSet local_fields;
    FieldSet global_fields;
    local_fields.add(local);
    global_fields.add(global);
    gather(local_fields,global_fields);
}
// ----------------------------------------------------------------------------



// ----------------------------------------------------------------------------
// Scatter FieldSet
// ----------------------------------------------------------------------------
void StructuredColumns::scatter(
    const FieldSet& global_fieldset,
    FieldSet& local_fieldset) const
{
  ASSERT(local_fieldset.size() == global_fieldset.size());

  for( size_t f=0; f<local_fieldset.size(); ++f ) {

    const Field& glb = global_fieldset[f];
    Field& loc = local_fieldset[f];
    const size_t nb_fields = 1;
    size_t root(0);
    glb.metadata().get("owner",root);

    if     ( loc.datatype() == array::DataType::kind<int>() ) {
      parallel::Field<int const> glb_field( array::make_storageview<int>(glb).data(),glb.stride(0));
      parallel::Field<int      > loc_field( array::make_storageview<int>(loc).data(),loc.stride(0));
      gather_scatter_->scatter( &glb_field, &loc_field, nb_fields, root );
    }
    else if( loc.datatype() == array::DataType::kind<long>() ) {
      parallel::Field<long const> glb_field( array::make_storageview<long>(glb).data(),glb.stride(0));
      parallel::Field<long      > loc_field( array::make_storageview<long>(loc).data(),loc.stride(0));
      gather_scatter_->scatter( &glb_field, &loc_field, nb_fields, root );
    }
    else if( loc.datatype() == array::DataType::kind<float>() ) {
      parallel::Field<float const> glb_field( array::make_storageview<float>(glb).data(),glb.stride(0));
      parallel::Field<float      > loc_field( array::make_storageview<float>(loc).data(),loc.stride(0));
      gather_scatter_->scatter( &glb_field, &loc_field, nb_fields, root );
    }
    else if( loc.datatype() == array::DataType::kind<double>() ) {
      parallel::Field<double const> glb_field( array::make_storageview<double>(glb).data(),glb.stride(0));
      parallel::Field<double      > loc_field( array::make_storageview<double>(loc).data(),loc.stride(0));
      gather_scatter_->scatter( &glb_field, &loc_field, nb_fields, root );
    }
    else throw eckit::Exception("datatype not supported",Here());

    glb.metadata().broadcast(loc.metadata(),root);
    loc.metadata().set("global",false);
  }
}
// ----------------------------------------------------------------------------



// ----------------------------------------------------------------------------
// Scatter Field
// ----------------------------------------------------------------------------
void StructuredColumns::scatter(
    const Field& global,
    Field& local) const
{
    FieldSet global_fields;
    FieldSet local_fields;
    global_fields.add(global);
    local_fields.add(local);
    scatter(global_fields, local_fields);
}
// ----------------------------------------------------------------------------






// ----------------------------------------------------------------------------
// Retrieve Global index from Local one
// ----------------------------------------------------------------------------
double StructuredColumns::y(
    size_t j) const
{
  return grid_.y(j+j_begin_);
}
// ----------------------------------------------------------------------------



// ----------------------------------------------------------------------------
// Retrieve Global index from Local one
// ----------------------------------------------------------------------------
double StructuredColumns::x(
    size_t i,
    size_t j) const
{
  return grid_.x(i+i_begin_[j],j+j_begin_);
}
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// Checksum FieldSet
// ----------------------------------------------------------------------------
std::string StructuredColumns::checksum(
    const FieldSet& fieldset) const
{
    eckit::MD5 md5;
    NOTIMP;
    return md5;
}
// ----------------------------------------------------------------------------



// ----------------------------------------------------------------------------
// Checksum Field
// ----------------------------------------------------------------------------
std::string StructuredColumns::checksum(
    const Field& field) const
{
    FieldSet fieldset;
    fieldset.add(field);
    return checksum(fieldset);
}
// ----------------------------------------------------------------------------

} // namespace detail

// ----------------------------------------------------------------------------

StructuredColumns::StructuredColumns() :
  FunctionSpace(),
  functionspace_(nullptr) {
}

StructuredColumns::StructuredColumns( const FunctionSpace& functionspace ) :
  FunctionSpace( functionspace ),
  functionspace_( dynamic_cast<const detail::StructuredColumns*>( get() ) ) {
}

StructuredColumns::StructuredColumns( const Grid& grid ) :
  FunctionSpace( new detail::StructuredColumns(grid) ),
  functionspace_( dynamic_cast<const detail::StructuredColumns*>( get() ) ) {
}

StructuredColumns::StructuredColumns( const Grid& grid, const grid::Partitioner& partitioner ) :
  FunctionSpace( new detail::StructuredColumns(grid,partitioner) ),
  functionspace_( dynamic_cast<const detail::StructuredColumns*>( get() ) ) {
}


Field StructuredColumns::createField(const std::string& name, array::DataType datatype, const eckit::Parametrisation& options ) const
{
  return functionspace_->createField(name,datatype,options);
}

Field StructuredColumns::createField(
    const std::string& name, array::DataType datatype,
    size_t levels, const eckit::Parametrisation& options) const
{
  return functionspace_->createField(name,datatype,levels,options);
}

void StructuredColumns::gather( const FieldSet& local, FieldSet& global ) const {
  functionspace_->gather(local,global);
}

void StructuredColumns::gather( const Field& local, Field& global ) const {
  functionspace_->gather(local,global);
}

void StructuredColumns::scatter( const FieldSet& global, FieldSet& local ) const {
  functionspace_->scatter(global,local);
}

void StructuredColumns::scatter( const Field& global, Field& local ) const {
  functionspace_->scatter(global,local);
}

std::string StructuredColumns::checksum( const FieldSet& fieldset ) const {
  return functionspace_->checksum(fieldset);
}

std::string StructuredColumns::checksum( const Field& field) const {
  return functionspace_->checksum(field);
}


// ----------------------------------------------------------------------------
// Fortran interfaces
// ----------------------------------------------------------------------------
extern "C"
{

const detail::StructuredColumns* atlas__functionspace__StructuredColumns__new__grid (const Grid::Implementation* grid)
{
  ATLAS_ERROR_HANDLING(
    return new detail::StructuredColumns( Grid(grid) );
  );
  return 0;
}

void atlas__functionspace__StructuredColumns__delete (detail::StructuredColumns* This)
{
  ATLAS_ERROR_HANDLING(
    ASSERT(This);
    delete This;
  );
}

field::FieldImpl* atlas__fs__StructuredColumns__create_field_name_kind (const detail::StructuredColumns* This, const char* name, int kind, const eckit::Parametrisation* options)
{
  ATLAS_ERROR_HANDLING(
    ASSERT(This);
    field::FieldImpl* field;
    {
      Field f = This->createField(std::string(name),array::DataType(kind),*options);
      field = f.get();
      field->attach();
    }
    field->detach();
    return field;
  );
  return 0;
}

field::FieldImpl* atlas__fs__StructuredColumns__create_field_name_kind_lev (const detail::StructuredColumns* This, const char* name, int kind, int levels, const eckit::Parametrisation* options)
{
  ATLAS_ERROR_HANDLING(
    ASSERT(This);
    field::FieldImpl* field;
    {
      Field f = This->createField(std::string(name),array::DataType(kind),levels,*options);
      field = f.get();
      field->attach();
    }
    field->detach();
    return field;
  );
  return 0;
}

void atlas__functionspace__StructuredColumns__gather (const detail::StructuredColumns* This, const field::FieldImpl* local, field::FieldImpl* global)
{
  ATLAS_ERROR_HANDLING(
    ASSERT(This);
    ASSERT(global);
    ASSERT(local);
    const Field l(local);
    Field g(global);
    This->gather(l,g);
  );
}

void atlas__functionspace__StructuredColumns__scatter (const detail::StructuredColumns* This, const field::FieldImpl* global, field::FieldImpl* local)
{
  ATLAS_ERROR_HANDLING(
    ASSERT(This);
    ASSERT(global);
    ASSERT(local);
    const Field g(global);
    Field l(local);
    This->scatter(g,l);
  );
}

void atlas__fs__StructuredColumns__checksum_fieldset(const detail::StructuredColumns* This, const field::FieldSetImpl* fieldset, char* &checksum, int &size, int &allocated)
{
  ASSERT(This);
  ASSERT(fieldset);
  ATLAS_ERROR_HANDLING(
    std::string checksum_str (This->checksum(fieldset));
    size = checksum_str.size();
    checksum = new char[size+1]; allocated = true;
    strcpy(checksum,checksum_str.c_str());
  );
}

void atlas__fs__StructuredColumns__checksum_field(const detail::StructuredColumns* This, const field::FieldImpl* field, char* &checksum, int &size, int &allocated)
{
  ASSERT(This);
  ASSERT(field);
  ATLAS_ERROR_HANDLING(
    std::string checksum_str (This->checksum(field));
    size = checksum_str.size();
    checksum = new char[size+1]; allocated = true;
    strcpy(checksum,checksum_str.c_str());
  );
}

}
// ----------------------------------------------------------------------------

} // namespace functionspace
} // namespace atlas

