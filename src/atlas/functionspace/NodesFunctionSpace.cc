/*
 * (C) Copyright 1996-2014 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include <algorithm>
#include "atlas/atlas_config.h"
#include "atlas/mpi/Collectives.h"
#include "atlas/Mesh.h"
#include "atlas/Nodes.h"
#include "atlas/mpl/HaloExchange.h"
#include "atlas/mpl/GatherScatter.h"
#include "atlas/util/IsGhost.h"
#include "atlas/functionspace/NodesFunctionSpace.h"
#include "atlas/field/FieldT.h"
#include "atlas/actions/BuildParallelFields.h"
#include "atlas/actions/BuildHalo.h"
#include "atlas/actions/BuildPeriodicBoundaries.h"
#include "atlas/atlas_omp.h"
#include "atlas/runtime/ErrorHandling.h"

#ifdef ATLAS_HAVE_FORTRAN
#define REMOTE_IDX_BASE 1
#else
#define REMOTE_IDX_BASE 0
#endif


namespace atlas {
namespace functionspace {

const std::string& nodes_str() { static std::string str("nodes"); return str; }

namespace {

template <typename T>
ArrayView<T,3> leveled_view(const Field &field)
{
  if( field.has_levels() )
    return ArrayView<T,3> ( field.data<T>(), make_shape(field.shape(0),field.shape(1),field.stride(1)) );
  else
    return ArrayView<T,3> ( field.data<T>(), make_shape(field.shape(0),1,field.stride(0)) );
}

template <typename T>
ArrayView<T,2> surface_view(const Field &field)
{
  return ArrayView<T,2> ( field.data<T>(), make_shape(field.shape(0),field.stride(0)) );
}

template <typename T>
ArrayView<T,2> leveled_scalar_view(const Field &field)
{
  if( field.has_levels() )
    return ArrayView<T,2> ( field.data<T>(), make_shape(field.shape(0),field.shape(1)) );
  else
    return ArrayView<T,2> ( field.data<T>(), make_shape(field.shape(0),1) );
}

template <typename T>
ArrayView<T,1> surface_scalar_view(const Field &field)
{
  return ArrayView<T,1> ( field.data<T>(), make_shape(field.size()) );
}


}


NodesFunctionSpace::NodesFunctionSpace(const std::string& name, Mesh& mesh, const Halo& halo)
  : next::FunctionSpace(name),
    mesh_(mesh),
    nodes_(mesh_.nodes()),
    halo_(halo),
    nb_nodes_(0),
    nb_nodes_global_(0),
    nb_nodes_global_foreach_rank_()
{
  actions::build_nodes_parallel_fields( mesh_.nodes() );
  actions::build_periodic_boundaries(mesh_);

  if( ! mesh_.halo_exchange().has(halo_name()) && halo_.size() > 0)
  {
    // Create new halo-exchange
    mpl::HaloExchange* halo_exchange = new mpl::HaloExchange( halo_name() );

    actions::build_halo(mesh_,halo_.size());

    actions::renumber_nodes_glb_idx(mesh_.nodes());

    Field& ridx = mesh_.nodes().remote_index();
    Field& part = mesh_.nodes().partition();

    std::stringstream ss;
    ss << "nb_nodes_including_halo["<<halo_.size()<<"]";
    mesh.metadata().get(ss.str(),nb_nodes_);

    halo_exchange->setup(part.data<int>(),ridx.data<int>(),REMOTE_IDX_BASE,nb_nodes_);

    // Store it in the mesh
    mesh_.halo_exchange().add(halo_exchange);
  }
  if( !nb_nodes_ ) {
    std::stringstream ss;
    ss << "nb_nodes_including_halo["<<halo_.size()<<"]";
    if( ! mesh.metadata().get(ss.str(),nb_nodes_) ) {
      nb_nodes_ = mesh_.nodes().metadata().get<size_t>("nb_owned");
    }
  }

  if( !mesh_.gather_scatter().has(gather_scatter_name()) )
  {
    // Create new gather_scatter
    mpl::GatherScatter* gather_scatter = new mpl::GatherScatter( gather_scatter_name() );

    Field& ridx = mesh_.nodes().remote_index();
    Field& part = mesh_.nodes().partition();
    Field& gidx = mesh_.nodes().global_index();

    util::IsGhost is_ghost(mesh_.nodes());
    std::vector<int> mask(mesh_.nodes().size());
    atlas_omp_parallel_for( size_t n=0; n<mask.size(); ++n ) {
      mask[n] = is_ghost(n) ? 1 : 0;

      // --> This would add periodic west-bc to the gather, but means that global-sums, means, etc are computed wrong
      //if( mask[j] == 1 && util::Topology::check(flags(j),util::Topology::BC) ) {
      //  mask[j] = 0;
      //}
    }
    gather_scatter->setup(part.data<int>(),ridx.data<int>(),REMOTE_IDX_BASE,gidx.data<gidx_t>(),mask.data(),nb_nodes_);

    // Store it in the mesh
    mesh_.gather_scatter().add(gather_scatter);
  }

  if( !mesh_.checksum().has(checksum_name()) )
  {
    // Create new checksum
    mpl::Checksum* checksum = new mpl::Checksum( checksum_name() );

    Field& ridx = mesh_.nodes().remote_index();
    Field& part = mesh_.nodes().partition();
    Field& gidx = mesh_.nodes().global_index();
    util::IsGhost is_ghost(mesh_.nodes());
    std::vector<int> mask(mesh_.nodes().size());
    atlas_omp_parallel_for( size_t n=0; n<mask.size(); ++n ) {
      mask[n] = is_ghost(n) ? 1 : 0;
    }
    checksum->setup(part.data<int>(),ridx.data<int>(),REMOTE_IDX_BASE,gidx.data<gidx_t>(),mask.data(),nb_nodes_);

    // Store it in the mesh
    mesh_.checksum().add(checksum);
  }

  nb_nodes_global_ = mesh_.gather_scatter().get(gather_scatter_name()).glb_dof();
  const std::vector<int>& glb_dofs = mesh_.gather_scatter().get(gather_scatter_name()).glb_dofs();
  nb_nodes_global_foreach_rank_.assign( glb_dofs.begin(), glb_dofs.end() );
}
NodesFunctionSpace::~NodesFunctionSpace() {}

size_t NodesFunctionSpace::nb_nodes() const
{
  return nb_nodes_;
}

size_t NodesFunctionSpace::nb_nodes_global() const
{
  return nb_nodes_global_;
}


std::vector<size_t> NodesFunctionSpace::nb_nodes_global_foreach_rank() const
{
  return nb_nodes_global_foreach_rank_;
}

std::string NodesFunctionSpace::halo_name() const
{
  std::stringstream ss; ss << "nodes_" << halo_.size();
  return ss.str();
}

std::string NodesFunctionSpace::gather_scatter_name() const
{
  return "nodes_gather_scatter";
}

std::string NodesFunctionSpace::checksum_name() const
{
  return "nodes_checksum";
}

Field* NodesFunctionSpace::createField(DataType datatype) const {
  Field* field = new Field(datatype,make_shape(nb_nodes()));
  field->set_functionspace(nodes_str());
  return field;
}

Field* NodesFunctionSpace::createField(DataType datatype, size_t levels) const {
  Field* field = new Field(datatype,make_shape(nb_nodes(),levels));
  field->set_levels(levels);
  field->set_functionspace(nodes_str());
  return field;
}

Field* NodesFunctionSpace::createField(const std::string& name,DataType datatype) const {
  Field* field = new Field(name,datatype,make_shape(nb_nodes()));
  field->set_functionspace(nodes_str());
  return field;
}

Field* NodesFunctionSpace::createField(const std::string& name,DataType datatype, size_t levels) const {
  Field* field = new Field(name,datatype,make_shape(nb_nodes(),levels));
  field->set_levels(levels);
  field->set_functionspace(nodes_str());
  return field;
}

Field* NodesFunctionSpace::createField(DataType datatype, const std::vector<size_t>& variables) const {
  std::vector<size_t> shape(1,nb_nodes());
  for( size_t i=0; i<variables.size(); ++i ) shape.push_back(variables[i]);
  Field* field = new Field(datatype,shape);
  field->set_functionspace(nodes_str());
  return field;
}

Field* NodesFunctionSpace::createField(DataType datatype, size_t levels, const std::vector<size_t>& variables) const {
  std::vector<size_t> shape(1,nb_nodes()); shape.push_back(levels);
  for( size_t i=0; i<variables.size(); ++i ) shape.push_back(variables[i]);
  Field* field = new Field(datatype,shape);
  field->set_levels(levels);
  field->set_functionspace(nodes_str());
  return field;
}

Field* NodesFunctionSpace::createField(const std::string& name,DataType datatype, const std::vector<size_t>& variables) const {
  std::vector<size_t> shape(1,nb_nodes());
  for( size_t i=0; i<variables.size(); ++i ) shape.push_back(variables[i]);
  Field* field = new Field(name,datatype,shape);
  field->set_functionspace(nodes_str());
  return field;
}

Field* NodesFunctionSpace::createField(const std::string& name, DataType datatype, size_t levels, const std::vector<size_t>& variables) const {
  std::vector<size_t> shape(1,nb_nodes()); shape.push_back(levels);
  for( size_t i=0; i<variables.size(); ++i ) shape.push_back(variables[i]);
  Field* field = new Field(name,datatype,shape);
  field->set_levels(levels);
  field->set_functionspace(nodes_str());
  return field;
}

Field* NodesFunctionSpace::createField(const Field& other) const {
  ArrayShape shape = other.shape();
  shape[0] = nb_nodes();
  Field* field = new Field(other.datatype(),shape);
  if( other.has_levels() )
    field->set_levels(field->shape(1));
  field->set_functionspace(nodes_str());
  return field;
}

Field* NodesFunctionSpace::createField(const std::string& name, const Field& other) const {
  ArrayShape shape = other.shape();
  shape[0] = nb_nodes();
  Field* field = new Field(name,other.datatype(),shape);
  if( other.has_levels() )
    field->set_levels(field->shape(1));
  field->set_functionspace(nodes_str());
  return field;
}

Field* NodesFunctionSpace::createGlobalField(DataType datatype) const {
  Field* field = new Field(datatype,make_shape(nb_nodes_global()));
  field->set_functionspace(nodes_str());
  return field;
}

Field* NodesFunctionSpace::createGlobalField(DataType datatype, size_t levels) const {
  Field* field = new Field(datatype,make_shape(nb_nodes_global(),levels));
  field->set_levels(levels);
  field->set_functionspace(nodes_str());
  return field;
}

Field* NodesFunctionSpace::createGlobalField(const std::string& name,DataType datatype) const {
  Field* field = new Field(name,datatype,make_shape(nb_nodes_global()));
  field->set_functionspace(nodes_str());
  return field;
}

Field* NodesFunctionSpace::createGlobalField(const std::string& name, DataType datatype, size_t levels) const {
  Field* field = new Field(name,datatype,make_shape(nb_nodes_global(),levels));
  field->set_levels(levels);
  field->set_functionspace(nodes_str());
  return field;
}

Field* NodesFunctionSpace::createGlobalField(DataType datatype, const std::vector<size_t>& variables) const {
  std::vector<size_t> shape(1,nb_nodes_global());
  for( size_t i=0; i<variables.size(); ++i ) shape.push_back(variables[i]);
  Field* field = new Field(datatype,shape);
  field->set_functionspace(nodes_str());
  return field;
}

Field* NodesFunctionSpace::createGlobalField(DataType datatype, size_t levels, const std::vector<size_t>& variables) const {
  std::vector<size_t> shape(1,nb_nodes_global()); shape.push_back(levels);
  for( size_t i=0; i<variables.size(); ++i ) shape.push_back(variables[i]);
  Field* field = new Field(datatype,shape);
  field->set_levels(levels);
  field->set_functionspace(nodes_str());
  return field;
}

Field* NodesFunctionSpace::createGlobalField(const std::string& name, DataType datatype, const std::vector<size_t>& variables) const {
  std::vector<size_t> shape(1,nb_nodes_global());
  for( size_t i=0; i<variables.size(); ++i ) shape.push_back(variables[i]);
  Field* field = new Field(name,datatype,shape);
  field->set_functionspace(nodes_str());
  return field;
}

Field* NodesFunctionSpace::createGlobalField(const std::string& name, DataType datatype, size_t levels, const std::vector<size_t>& variables) const {
  std::vector<size_t> shape(1,nb_nodes_global()); shape.push_back(levels);
  for( size_t i=0; i<variables.size(); ++i ) shape.push_back(variables[i]);
  Field* field = new Field(name,datatype,shape);
  field->set_levels(levels);
  field->set_functionspace(nodes_str());
  return field;
}

Field* NodesFunctionSpace::createGlobalField(const Field& other) const {
  ArrayShape shape = other.shape();
  shape[0] = nb_nodes_global();
  Field* field = new Field(other.datatype(),shape);
  if( other.has_levels() )
    field->set_levels(field->shape(1));
  field->set_functionspace(nodes_str());
  return field;
}

Field* NodesFunctionSpace::createGlobalField(const std::string& name,const Field& other) const {
  ArrayShape shape = other.shape();
  shape[0] = nb_nodes_global();
  Field* field = new Field(name,other.datatype(),shape);
  if( other.has_levels() )
    field->set_levels(field->shape(1));
  field->set_functionspace(nodes_str());
  return field;
}

void NodesFunctionSpace::haloExchange( FieldSet& fieldset ) const
{
  if( halo_.size() ) {
    const mpl::HaloExchange& halo_exchange = mesh_.halo_exchange().get(halo_name());
    for( size_t f=0; f<fieldset.size(); ++f ) {
      const Field& field = fieldset[f];
      ArrayStrides strides = make_strides(field.stride(0),1);
      ArrayShape   shape   = make_shape(field.shape(0),field.stride(0));
      if     ( field.datatype() == DataType::kind<int>() ) {
        ArrayView<int,2> view(field.data<int>(),strides.data(),shape.data());
        halo_exchange.execute( view );
      }
      else if( field.datatype() == DataType::kind<long>() ) {
        ArrayView<long,2> view(field.data<long>(),strides.data(),shape.data());
        halo_exchange.execute( view );
      }
      else if( field.datatype() == DataType::kind<float>() ) {
        ArrayView<float,2> view(field.data<float>(),strides.data(),shape.data());
        halo_exchange.execute( view );
      }
      else if( field.datatype() == DataType::kind<double>() ) {
        ArrayView<double,2> view(field.data<double>(),strides.data(),shape.data());
        halo_exchange.execute( view );
      }
      else throw eckit::Exception("datatype not supported",Here());
    }
  }
}
void NodesFunctionSpace::haloExchange( Field& field ) const
{
  if( halo_.size() ) {
    FieldSet fieldset;
    fieldset.add(field);
    haloExchange(fieldset);
  }
}
const mpl::HaloExchange& NodesFunctionSpace::halo_exchange() const
{
  return mesh_.halo_exchange().get(halo_name());
}


void NodesFunctionSpace::gather( const FieldSet& local_fieldset, FieldSet& global_fieldset ) const
{
  const mpl::GatherScatter& gather_scatter = mesh_.gather_scatter().get(gather_scatter_name());

  ASSERT(local_fieldset.size() == global_fieldset.size());

  for( size_t f=0; f<local_fieldset.size(); ++f ) {

    const Field& loc = local_fieldset[f];
    Field& glb = global_fieldset[f];

    if     ( loc.datatype() == DataType::kind<int>() ) {
      mpl::Field<int const> loc_field(loc.data<int>(),loc.stride(0));
      mpl::Field<int      > glb_field(glb.data<int>(),glb.stride(0));
      gather_scatter.gather( &loc_field, &glb_field, 1 );
    }
    else if( loc.datatype() == DataType::kind<long>() ) {
      mpl::Field<long const> loc_field(loc.data<long>(),loc.stride(0));
      mpl::Field<long      > glb_field(glb.data<long>(),glb.stride(0));
      gather_scatter.gather( &loc_field, &glb_field, 1 );
    }
    else if( loc.datatype() == DataType::kind<float>() ) {
      mpl::Field<float const> loc_field(loc.data<float>(),loc.stride(0));
      mpl::Field<float      > glb_field(glb.data<float>(),glb.stride(0));
      gather_scatter.gather( &loc_field, &glb_field, 1 );
    }
    else if( loc.datatype() == DataType::kind<double>() ) {
      mpl::Field<double const> loc_field(loc.data<double>(),loc.stride(0));
      mpl::Field<double      > glb_field(glb.data<double>(),glb.stride(0));
      gather_scatter.gather( &loc_field, &glb_field, 1 );
    }
    else throw eckit::Exception("datatype not supported",Here());
  }
}
void NodesFunctionSpace::gather( const Field& local, Field& global ) const
{
  FieldSet local_fields;
  FieldSet global_fields;
  local_fields.add(local);
  global_fields.add(global);
  gather(local_fields,global_fields);
}
const mpl::GatherScatter& NodesFunctionSpace::gather() const
{
  return mesh_.gather_scatter().get(gather_scatter_name());
}
const mpl::GatherScatter& NodesFunctionSpace::scatter() const
{
  return mesh_.gather_scatter().get(gather_scatter_name());
}


void NodesFunctionSpace::scatter( const FieldSet& global_fieldset, FieldSet& local_fieldset ) const
{
  const mpl::GatherScatter& gather_scatter = mesh_.gather_scatter().get(gather_scatter_name());

  ASSERT(local_fieldset.size() == global_fieldset.size());

  for( size_t f=0; f<local_fieldset.size(); ++f ) {

    const Field& glb = global_fieldset[f];
    Field& loc = local_fieldset[f];

    if     ( loc.datatype() == DataType::kind<int>() ) {
      mpl::Field<int const> glb_field(glb.data<int>(),glb.stride(0));
      mpl::Field<int      > loc_field(loc.data<int>(),loc.stride(0));
      gather_scatter.scatter( &glb_field, &loc_field, 1 );
    }
    else if( loc.datatype() == DataType::kind<long>() ) {
      mpl::Field<long const> glb_field(glb.data<long>(),glb.stride(0));
      mpl::Field<long      > loc_field(loc.data<long>(),loc.stride(0));
      gather_scatter.scatter( &glb_field, &loc_field, 1 );
    }
    else if( loc.datatype() == DataType::kind<float>() ) {
      mpl::Field<float const> glb_field(glb.data<float>(),glb.stride(0));
      mpl::Field<float      > loc_field(loc.data<float>(),loc.stride(0));
      gather_scatter.scatter( &glb_field, &loc_field, 1 );
    }
    else if( loc.datatype() == DataType::kind<double>() ) {
      mpl::Field<double const> glb_field(glb.data<double>(),glb.stride(0));
      mpl::Field<double      > loc_field(loc.data<double>(),loc.stride(0));
      gather_scatter.scatter( &glb_field, &loc_field, 1 );
    }
    else throw eckit::Exception("datatype not supported",Here());
  }
}
void NodesFunctionSpace::scatter( const Field& global, Field& local ) const
{
  FieldSet global_fields;
  FieldSet local_fields;
  global_fields.add(global);
  local_fields.add(local);
  scatter(global_fields,local_fields);
}

namespace {
template <typename T>
std::string checksum_3d_field(const mpl::Checksum& checksum, const Field& field )
{
  ArrayView<T,3> values = leveled_view<T>(field);
  Array<T> surface_field ( make_shape(values.shape(0),values.shape(2) ) );
  ArrayView<T,2> surface(surface_field);
  for( size_t n=0; n<values.shape(0); ++n ) {
    for( size_t j=0; j<surface.shape(1); ++j )
    {
      surface(n,j) = 0.;
      for( size_t l=0; l<values.shape(1);++l )
        surface(n,j) += values(n,l,j);
    }
  }
  return checksum.execute( surface.data(), surface.stride(0) );
}
}

std::string NodesFunctionSpace::checksum( const FieldSet& fieldset ) const {
  const mpl::Checksum& checksum = mesh().checksum().get(checksum_name());
  eckit::MD5 md5;
  for( size_t f=0; f<fieldset.size(); ++f ) {
    const Field& field=fieldset[f];
    if     ( field.datatype() == DataType::kind<int>() )
      md5 << checksum_3d_field<int>(checksum,field);
    else if( field.datatype() == DataType::kind<long>() )
      md5 << checksum_3d_field<long>(checksum,field);
    else if( field.datatype() == DataType::kind<float>() )
      md5 << checksum_3d_field<float>(checksum,field);
    else if( field.datatype() == DataType::kind<double>() )
      md5 << checksum_3d_field<double>(checksum,field);
    else throw eckit::Exception("datatype not supported",Here());
  }
  return md5;
}
std::string NodesFunctionSpace::checksum( const Field& field ) const {
  FieldSet fieldset;
  fieldset.add(field);
  return checksum(fieldset);
}

const mpl::Checksum& NodesFunctionSpace::checksum() const
{
  return mesh_.checksum().get(checksum_name());
}



//std::string NodesFunctionSpace::checksum( const FieldSet& fieldset ) const {
//  const mpl::Checksum& checksum = mesh_.checksum().get(checksum_name());

//  eckit::MD5 md5;
//  for( size_t f=0; f<fieldset.size(); ++f ) {
//    const Field& field=fieldset[f];
//    if     ( field.datatype() == DataType::kind<int>() )
//      md5 << checksum.execute( field.data<int>(), field.stride(0) );
//    else if( field.datatype() == DataType::kind<long>() )
//      md5 << checksum.execute( field.data<long>(), field.stride(0) );
//    else if( field.datatype() == DataType::kind<float>() )
//      md5 << checksum.execute( field.data<float>(), field.stride(0) );
//    else if( field.datatype() == DataType::kind<double>() )
//      md5 << checksum.execute( field.data<double>(), field.stride(0) );
//    else throw eckit::Exception("datatype not supported",Here());
//  }
//  return md5;
//}
//std::string NodesFunctionSpace::checksum( const Field& field ) const {
//  FieldSet fieldset;
//  fieldset.add(field);
//  return checksum(fieldset);
//}

namespace { inline double sqr(const double& val) { return val*val; } }

namespace detail { // Collectives implementation



template< typename T >
void dispatch_sum( const NodesFunctionSpace& fs, const Field& field, T& result, size_t& N )
{
  const util::IsGhost is_ghost(fs.nodes());
  const ArrayView<T,2> arr = leveled_scalar_view<T>( field );
  T local_sum = 0;
  atlas_omp_pragma( omp parallel for default(shared) reduction(+:local_sum) )
  for( size_t n=0; n<arr.shape(0); ++n ) {
    if( ! is_ghost(n) ) {
      for( size_t l=0; l<arr.shape(1); ++l )
        local_sum += arr(n,l);
    }
  }
  eckit::mpi::all_reduce(local_sum,result,eckit::mpi::sum());
  N = fs.nb_nodes_global_foreach_rank()[0] * arr.shape(1);
}

template< typename T >
void sum( const NodesFunctionSpace& fs , const Field& field, T& result, size_t& N )
{
  if( field.datatype() == DataType::kind<T>() ) {
    return dispatch_sum(fs,field,result,N);
  }
  else
  {
    switch( field.datatype().kind() )
    {
      case DataType::KIND_INT32 : {
        int tmp;
        dispatch_sum(fs,field,tmp,N);
        result = tmp;
        return;
      }
      case DataType::KIND_INT64 : {
        long tmp;
        dispatch_sum(fs,field,tmp,N);
        result = tmp;
        return;
      }
      case DataType::KIND_REAL32 : {
        float tmp;
        dispatch_sum(fs,field,tmp,N);
        result = tmp;
        return;
      }
      case DataType::KIND_REAL64 : {
        double tmp;
        dispatch_sum(fs,field,tmp,N);
        result = tmp;
        return;
      }
      default: throw eckit::Exception("datatype not supported",Here());
    }
  }
}



template< typename T >
void dispatch_sum( const NodesFunctionSpace& fs, const Field& field, std::vector<T>& result, size_t& N )
{
  const ArrayView<T,3> arr = leveled_view<T>(field);
  const util::IsGhost is_ghost(fs.nodes());
  const size_t nvar = arr.shape(2);
  std::vector<T> local_sum(nvar,0);

  atlas_omp_parallel
  {
    std::vector<T> local_sum_private(nvar,0);
    atlas_omp_for( size_t n=0; n<arr.shape(0); ++n )
    {
      if( ! is_ghost(n) ) {
        for( size_t l=0; l<arr.shape(1); ++l ) {
          for( size_t j=0; j<arr.shape(2); ++j ) {
            local_sum_private[j] += arr(n,l,j);
          }
        }
      }
    }
    atlas_omp_critical
    {
      for( size_t j=0; j<nvar; ++j ) {
        local_sum[j] += local_sum_private[j];
      }
    }
  }
  eckit::mpi::all_reduce(local_sum,result,eckit::mpi::sum());
  N = fs.nb_nodes_global_foreach_rank()[0] * arr.shape(1);
}

template< typename T >
void sum( const NodesFunctionSpace& fs , const Field& field, std::vector<T>& result, size_t& N )
{
  if( field.datatype() == DataType::kind<T>() ) {
    return dispatch_sum(fs,field,result,N);
  }
  else
  {
    switch( field.datatype().kind() )
    {
      case DataType::KIND_INT32 : {
        std::vector<int> tmp;
        dispatch_sum(fs,field,tmp,N);
        result.assign(tmp.begin(),tmp.end());
        return;
      }
      case DataType::KIND_INT64 : {
        std::vector<long> tmp;
        dispatch_sum(fs,field,tmp,N);
        result.assign(tmp.begin(),tmp.end());
        return;
      }
      case DataType::KIND_REAL32 : {
        std::vector<float> tmp;
        dispatch_sum(fs,field,tmp,N);
        result.assign(tmp.begin(),tmp.end());
        return;
      }
      case DataType::KIND_REAL64 : {
        std::vector<double> tmp;
        dispatch_sum(fs,field,tmp,N);
        result.assign(tmp.begin(),tmp.end());
        return;
      }
      default: throw eckit::Exception("datatype not supported",Here());
    }
  }
}


template< typename T >
void dispatch_sum_per_level( const NodesFunctionSpace& fs, const Field& field, Field& sum, size_t& N )
{
  ArrayShape shape;
  shape.reserve(field.rank()-1);
  for( size_t j=1; j<field.rank(); ++j )
    shape.push_back(field.shape(j));
  sum.resize(shape);

  const ArrayView<T,3> arr = leveled_view<T>(field);
  ArrayView<T,2> sum_per_level( sum.data<T>(), make_shape(sum.shape(0),sum.stride(0)) );
  sum_per_level = 0;
  const util::IsGhost is_ghost(fs.nodes());

  atlas_omp_parallel
  {
    Array<T> sum_per_level_private(sum_per_level.shape(0),sum_per_level.shape(1));
    ArrayView<T> sum_per_level_private_view(sum_per_level_private); sum_per_level_private_view = 0.;
    atlas_omp_for( size_t n=0; n<arr.shape(0); ++n )
    {
      if( ! is_ghost(n) ) {
        for( size_t l=0; l<arr.shape(1); ++l ) {
          for( size_t j=0; j<arr.shape(2); ++j ) {
            sum_per_level_private(l,j) += arr(n,l,j);
          }
        }
      }
    }
    atlas_omp_critical
    {
      for( size_t l=0; l<sum_per_level_private.shape(0); ++l ) {
        for( size_t j=0; j<sum_per_level_private.shape(1); ++j ) {
          sum_per_level(l,j) += sum_per_level_private(l,j);
        }
      }
    }
  }
  eckit::mpi::all_reduce(sum_per_level.data(),sum.size(),eckit::mpi::sum());
  N = fs.nb_nodes_global_foreach_rank()[0];
}

void sum_per_level( const NodesFunctionSpace& fs, const Field& field, Field& sum, size_t& N )
{
  if( field.datatype() != sum.datatype() ) {
    throw eckit::Exception("Field and sum are not of same datatype.",Here());
  }
  switch( field.datatype().kind() )
  {
    case DataType::KIND_INT32 :
      return dispatch_sum_per_level<int>(fs,field,sum,N);
    case DataType::KIND_INT64 :
      return dispatch_sum_per_level<long>(fs,field,sum,N);
    case DataType::KIND_REAL32 :
      return dispatch_sum_per_level<float>(fs,field,sum,N);
    case DataType::KIND_REAL64 :
      return dispatch_sum_per_level<double>(fs,field,sum,N);
    default: throw eckit::Exception("datatype not supported",Here());
  }
}


template< typename DATATYPE >
void dispatch_order_independent_sum_2d( const NodesFunctionSpace& fs , const Field& field, DATATYPE& result, size_t& N )
{
  size_t root = 0;
  Field::Ptr global( fs.createGlobalField(field) );
  fs.gather(field,*global);
  result = std::accumulate(global->data<DATATYPE>(),global->data<DATATYPE>()+global->size(),0.);
  eckit::mpi::broadcast(result,root);
  N = fs.nb_nodes_global_foreach_rank()[0];
}

template< typename T >
void dispatch_order_independent_sum( const NodesFunctionSpace& fs , const Field& field, T& result, size_t& N )
{
  if( field.has_levels() )
  {
    const ArrayView<T,2> arr = leveled_scalar_view<T>(field);

    Field::Ptr surface_field( fs.createField<T>() );
    ArrayView<T,1> surface = surface_scalar_view<T>( *surface_field );

    for( size_t n=0; n<arr.shape(0); ++n ) {
      for( size_t l=0; l<arr.shape(1); ++l ) {
          surface(n) += arr(n,l);
      }
    }
    dispatch_order_independent_sum_2d( fs, *surface_field, result, N );
    N *= arr.shape(1);
  }
  else
  {
    dispatch_order_independent_sum_2d( fs, field, result, N );
  }
}

template< typename T >
void order_independent_sum( const NodesFunctionSpace& fs , const Field& field, T& result, size_t& N )
{
  if( field.datatype() == DataType::kind<T>() ) {
    return dispatch_order_independent_sum(fs,field,result,N);
  }
  else
  {
    switch( field.datatype().kind() )
    {
      case DataType::KIND_INT32 : {
        int tmp;
        dispatch_order_independent_sum(fs,field,tmp,N);
        result = tmp;
        return;
      }
      case DataType::KIND_INT64 : {
        long tmp;
        dispatch_order_independent_sum(fs,field,tmp,N);
        result = tmp;
        return;
      }
      case DataType::KIND_REAL32 : {
        float tmp;
        dispatch_order_independent_sum(fs,field,tmp,N);
        result = tmp;
        return;
      }
      case DataType::KIND_REAL64 : {
        double tmp;
        dispatch_order_independent_sum(fs,field,tmp,N);
        result = tmp;
        return;
      }
      default: throw eckit::Exception("datatype not supported",Here());
    }
  }
}

template< typename DATATYPE >
void dispatch_order_independent_sum_2d( const NodesFunctionSpace& fs, const Field& field, std::vector<DATATYPE>& result, size_t& N )
{
  size_t nvar = field.stride(0);
  result.resize(nvar);
  for( size_t j=0; j<nvar; ++j ) result[j] = 0.;
  size_t root = 0;
  Field::Ptr global( fs.createGlobalField(field) );
  fs.gather(field,*global);
  if( eckit::mpi::rank() == 0 ) {
    const ArrayView<DATATYPE,2> glb( global->data<DATATYPE>(), make_shape(global->shape(0),global->stride(0)) );
    for( size_t n=0; n<fs.nb_nodes_global(); ++n ) {
      for( size_t j=0; j<nvar; ++j ) {
        result[j] += glb(n,j);
      }
    }
  }
  eckit::mpi::broadcast(result,root);
  N = fs.nb_nodes_global_foreach_rank()[0];
}

template< typename T >
void dispatch_order_independent_sum( const NodesFunctionSpace& fs, const Field& field, std::vector<T>& result, size_t& N )
{
  if( field.has_levels() )
  {
    const size_t nvar = field.stride(1);
    const ArrayView<T,3> arr = leveled_view<T>(field);

    Field::Ptr surface_field( fs.createField<T>(make_shape(nvar)) );
    ArrayView<T,2> surface = surface_view<T>( *surface_field );

    for( size_t n=0; n<arr.shape(0); ++n ) {
      for( size_t l=0; l<arr.shape(1); ++l ) {
        for( size_t j=0; j<arr.shape(2); ++j ) {
          surface(n,j) += arr(n,l,j);
        }
      }
    }

    dispatch_order_independent_sum_2d( fs, *surface_field, result, N );
    N *= arr.shape(1);
  }
  else
  {
    dispatch_order_independent_sum_2d( fs, field, result, N );
  }
}

template< typename T >
void order_independent_sum( const NodesFunctionSpace& fs, const Field& field, std::vector<T>& result, size_t& N )
{
  if( field.datatype() == DataType::kind<T>() ) {
    return dispatch_order_independent_sum(fs,field,result,N);
  }
  else
  {
    switch( field.datatype().kind() )
    {
      case DataType::KIND_INT32 : {
          std::vector<int> tmp;
          dispatch_order_independent_sum(fs,field,tmp,N);
          result.assign(tmp.begin(),tmp.end());
          return;
      }
      case DataType::KIND_INT64 : {
        std::vector<long> tmp;
        dispatch_order_independent_sum(fs,field,tmp,N);
        result.assign(tmp.begin(),tmp.end());
        return;
      }
      case DataType::KIND_REAL32 : {
        std::vector<float> tmp;
        dispatch_order_independent_sum(fs,field,tmp,N);
        result.assign(tmp.begin(),tmp.end());
        return;
      }
      case DataType::KIND_REAL64 : {
        std::vector<double> tmp;
        dispatch_order_independent_sum(fs,field,tmp,N);
        result.assign(tmp.begin(),tmp.end());
        return;
      }
      default: throw eckit::Exception("datatype not supported",Here());
    }
  }
}


template< typename T >
void dispatch_order_independent_sum_per_level( const NodesFunctionSpace& fs, const Field& field, Field& sumfield, size_t& N )
{
  ArrayShape shape;
  shape.reserve(field.rank()-1);
  for( size_t j=1; j<field.rank(); ++j )
    shape.push_back(field.shape(j));
  sumfield.resize(shape);

  ArrayView<T,2> sum ( sumfield.data<T>(), make_shape(sumfield.shape(0),sumfield.stride(0)) );
  sum = 0.;

  eckit::Log::info() << field << std::endl;
  eckit::Log::info() << sumfield << std::endl;

  size_t root = 0;
  Field::Ptr global( fs.createGlobalField(field) );

  eckit::Log::info() << *global << std::endl;

  fs.gather(field,*global);
  if( eckit::mpi::rank() == 0 ) {
    const ArrayView<T,3> glb = leveled_view<T>(*global);

    DEBUG_VAR(glb.shape(0));
    DEBUG_VAR(glb.shape(1));
    DEBUG_VAR(glb.shape(2));
    DEBUG_VAR(sum.shape(0));
    DEBUG_VAR(sum.shape(1));

    for( size_t n=0; n<glb.shape(0); ++n ) {
      for( size_t l=0; l<glb.shape(1); ++l ) {
        for( size_t j=0; j<glb.shape(2); ++j ) {
          sum(l,j) += glb(n,l,j);
        }
      }
    }
  }
  eckit::mpi::broadcast(sumfield.data<T>(),sumfield.size(),root);
  N = fs.nb_nodes_global_foreach_rank()[0];
}

void order_independent_sum_per_level( const NodesFunctionSpace& fs, const Field& field, Field& sum, size_t& N )
{
  if( field.datatype() != sum.datatype() ) {
    throw eckit::Exception("Field and sum are not of same datatype.",Here());
  }
  switch( field.datatype().kind() )
  {
    case DataType::KIND_INT32 :
      return dispatch_order_independent_sum_per_level<int>(fs,field,sum,N);
    case DataType::KIND_INT64 :
      return dispatch_order_independent_sum_per_level<long>(fs,field,sum,N);
    case DataType::KIND_REAL32 :
      return dispatch_order_independent_sum_per_level<float>(fs,field,sum,N);
    case DataType::KIND_REAL64 :
      return dispatch_order_independent_sum_per_level<double>(fs,field,sum,N);
    default: throw eckit::Exception("datatype not supported",Here());
  }
}

template< typename T >
void dispatch_minimum( const NodesFunctionSpace& fs, const Field& field, std::vector<T>& min )
{
  const ArrayView<T,3> arr = leveled_view<T>(field);
  const size_t nvar = arr.shape(2);
  min.resize(nvar);
  std::vector<T> local_minimum(nvar,std::numeric_limits<T>::max());
  atlas_omp_parallel
  {
    std::vector<T> local_minimum_private(nvar,std::numeric_limits<T>::max());
    atlas_omp_for( size_t n=0; n<arr.shape(0); ++n ) {
      for( size_t l=0; l<arr.shape(1); ++l ) {
        for( size_t j=0; j<arr.shape(2); ++j ) {
          local_minimum_private[j] = std::min(arr(n,l,j),local_minimum_private[j]);
        }
      }
    }
    atlas_omp_critical
    {
      for( size_t j=0; j<arr.shape(2); ++j ) {
        local_minimum[j] = std::min(local_minimum_private[j],local_minimum[j]);
      }
    }
  }
  eckit::mpi::all_reduce(local_minimum,min,eckit::mpi::min());
}

template< typename T >
void minimum( const NodesFunctionSpace& fs, const Field& field, std::vector<T>& min )
{
  if( field.datatype() == DataType::kind<T>() ) {
    return dispatch_minimum(fs,field,min);
  }
  else
  {
    switch( field.datatype().kind() )
    {
      case DataType::KIND_INT32 : {
        std::vector<int> tmp;
        dispatch_minimum(fs,field,tmp);
        min.assign(tmp.begin(),tmp.end());
        return;
      }
      case DataType::KIND_INT64 : {
        std::vector<long> tmp;
        dispatch_minimum(fs,field,tmp);
        min.assign(tmp.begin(),tmp.end());
        return;
      }
      case DataType::KIND_REAL32 : {
        std::vector<float> tmp;
        dispatch_minimum(fs,field,tmp);
        min.assign(tmp.begin(),tmp.end());
        return;
      }
      case DataType::KIND_REAL64 : {
        std::vector<double> tmp;
        dispatch_minimum(fs,field,tmp);
        min.assign(tmp.begin(),tmp.end());
        return;
      }
      default: throw eckit::Exception("datatype not supported",Here());
    }
  }
}

template< typename T >
void dispatch_maximum( const NodesFunctionSpace& fs, const Field& field, std::vector<T>& max )
{
  const ArrayView<T,3> arr = leveled_view<T>(field);
  const size_t nvar = arr.shape(2);
  max.resize(nvar);
  std::vector<T> local_maximum(nvar,-std::numeric_limits<T>::max());
  atlas_omp_parallel
  {
    std::vector<T> local_maximum_private(nvar,-std::numeric_limits<T>::max());
    atlas_omp_for( size_t n=0; n<arr.shape(0); ++n ) {
      for( size_t l=0; l<arr.shape(1); ++l ) {
        for( size_t j=0; j<nvar; ++j ) {
          local_maximum_private[j] = std::max(arr(n,l,j),local_maximum_private[j]);
        }
      }
    }
    atlas_omp_critical
    {
      for( size_t j=0; j<nvar; ++j ) {
        local_maximum[j] = std::max(local_maximum_private[j],local_maximum[j]);
      }
    }
  }
  eckit::mpi::all_reduce(local_maximum,max,eckit::mpi::max());
}

template< typename T >
void maximum( const NodesFunctionSpace& fs, const Field& field, std::vector<T>& max )
{
  if( field.datatype() == DataType::kind<T>() ) {
    return dispatch_maximum(fs,field,max);
  }
  else
  {
    switch( field.datatype().kind() )
    {
      case DataType::KIND_INT32 : {
        std::vector<int> tmp;
        dispatch_maximum(fs,field,tmp);
        max.assign(tmp.begin(),tmp.end());
        return;
      }
      case DataType::KIND_INT64 : {
        std::vector<long> tmp;
        dispatch_maximum(fs,field,tmp);
        max.assign(tmp.begin(),tmp.end());
        return;
      }
      case DataType::KIND_REAL32 : {
        std::vector<float> tmp;
        dispatch_maximum(fs,field,tmp);
        max.assign(tmp.begin(),tmp.end());
        return;
      }
      case DataType::KIND_REAL64 : {
        std::vector<double> tmp;
        dispatch_maximum(fs,field,tmp);
        max.assign(tmp.begin(),tmp.end());
        return;
      }
      default: throw eckit::Exception("datatype not supported",Here());
    }
  }
}

template< typename T >
void minimum( const NodesFunctionSpace& fs, const Field& field, T& min )
{
  std::vector<T> v;
  minimum(fs,field,v);
  min = v[0];
}

template< typename T >
void maximum( const NodesFunctionSpace& fs, const Field& field, T& max )
{
  std::vector<T> v;
  maximum(fs,field,v);
  max = v[0];
}

template< typename T >
void dispatch_minimum_per_level( const NodesFunctionSpace& fs, const Field& field, Field& min_field )
{
  ArrayShape shape;
  shape.reserve(field.rank()-1);
  for( size_t j=1; j<field.rank(); ++j )
    shape.push_back(field.shape(j));
  min_field.resize(shape);
  const size_t nvar = field.stride(1);
  ArrayView<T,2> min( min_field.data<T>(), make_shape(min_field.shape(0),min_field.stride(0)) );
  min = std::numeric_limits<T>::max();
  const ArrayView<T,3> arr = leveled_view<T>(field);
  atlas_omp_parallel
  {
    Array<T> min_private(min.shape(0),min.shape(1));
    ArrayView<T> min_private_view(min_private); min_private_view = std::numeric_limits<T>::max();
    atlas_omp_for( size_t n=0; n<arr.shape(0); ++n ) {
      for( size_t l=0; l<arr.shape(1); ++l ) {
        for( size_t j=0; j<nvar; ++j ) {
          min_private(l,j) = std::min(arr(n,l,j),min_private(l,j));
        }
      }
    }
    atlas_omp_critical
    {
      for( size_t l=0; l<arr.shape(1); ++l ) {
        for( size_t j=0; j<nvar; ++j ) {
          min(l,j) = std::min(min_private(l,j),min(l,j));
        }
      }
    }
  }
  eckit::mpi::all_reduce(min_field.data<T>(),min_field.size(),eckit::mpi::min());
}

void minimum_per_level( const NodesFunctionSpace& fs, const Field& field, Field& min )
{
  if( field.datatype() != min.datatype() ) {
    throw eckit::Exception("Field and min are not of same datatype.",Here());
  }
  switch( field.datatype().kind() )
  {
    case DataType::KIND_INT32 :
      return dispatch_minimum_per_level<int>(fs,field,min);
    case DataType::KIND_INT64 :
      return dispatch_minimum_per_level<long>(fs,field,min);
    case DataType::KIND_REAL32 :
      return dispatch_minimum_per_level<float>(fs,field,min);
    case DataType::KIND_REAL64 :
      return dispatch_minimum_per_level<double>(fs,field,min);
    default: throw eckit::Exception("datatype not supported",Here());
  }
}

template< typename T >
void dispatch_maximum_per_level( const NodesFunctionSpace& fs, const Field& field, Field& max_field )
{
  ArrayShape shape;
  shape.reserve(field.rank()-1);
  for( size_t j=1; j<field.rank(); ++j )
    shape.push_back(field.shape(j));
  max_field.resize(shape);
  const size_t nvar = field.stride(1);
  ArrayView<T,2> max = ArrayView<T,2>( max_field.data<T>(), make_shape(max_field.shape(0),max_field.stride(0)) );
  max = -std::numeric_limits<T>::max();
  const ArrayView<T,3> arr = leveled_view<T>(field);
  atlas_omp_parallel
  {
    Array<T> max_private(max.shape(0),max.shape(1));
    ArrayView<T> max_private_view(max_private); max_private_view = -std::numeric_limits<T>::max();
    atlas_omp_for( size_t n=0; n<arr.shape(0); ++n ) {
      for( size_t l=0; l<arr.shape(1); ++l ) {
        for( size_t j=0; j<nvar; ++j ) {
          max_private(l,j) = std::max(arr(n,l,j),max_private(l,j));
        }
      }
    }
    atlas_omp_critical
    {
      for( size_t l=0; l<arr.shape(1); ++l ) {
        for( size_t j=0; j<nvar; ++j ) {
          max(l,j) = std::max(max_private(l,j),max(l,j));
        }
      }
    }
  }
  eckit::mpi::all_reduce(max_field.data<T>(),max_field.size(),eckit::mpi::max());
}

void maximum_per_level( const NodesFunctionSpace& fs, const Field& field, Field& max )
{
  if( field.datatype() != max.datatype() ) {
    throw eckit::Exception("Field and max are not of same datatype.",Here());
  }
  switch( field.datatype().kind() )
  {
    case DataType::KIND_INT32 :
      return dispatch_maximum_per_level<int>(fs,field,max);
    case DataType::KIND_INT64 :
      return dispatch_maximum_per_level<long>(fs,field,max);
    case DataType::KIND_REAL32 :
      return dispatch_maximum_per_level<float>(fs,field,max);
    case DataType::KIND_REAL64 :
      return dispatch_maximum_per_level<double>(fs,field,max);
    default: throw eckit::Exception("datatype not supported",Here());
  }
}


template< typename T >
void dispatch_minimum_and_location( const NodesFunctionSpace& fs, const Field& field, std::vector<T>& min, std::vector<gidx_t>& glb_idx, std::vector<size_t>& level )
{
  ArrayView<T,3> arr = leveled_view<T>(field);
  size_t nvar = arr.shape(2);
  min.resize(nvar);
  glb_idx.resize(nvar);
  level.resize(nvar);
  std::vector<T> local_minimum(nvar,std::numeric_limits<T>::max());
  std::vector<size_t> loc_node(nvar);
  std::vector<size_t> loc_level(nvar);
  atlas_omp_parallel
  {
    std::vector<T> local_minimum_private(nvar,std::numeric_limits<T>::max());
    std::vector<size_t> loc_node_private(nvar);
    std::vector<size_t> loc_level_private(nvar);
    atlas_omp_for( size_t n=0; n<arr.shape(0); ++n ) {
      for( size_t l=0; l<arr.shape(1); ++l ) {
        for( size_t j=0; j<nvar; ++j ) {
          if( arr(n,l,j) < local_minimum_private[j] ) {
            local_minimum_private[j] = arr(n,l,j);
            loc_node_private[j] = n;
            loc_level_private[j] = l;
          }
        }
      }
    }
    atlas_omp_critical_ordered
    {
      for( size_t l=0; l<arr.shape(1); ++l ) {
        for( size_t j=0; j<nvar; ++j ) {
          if( local_minimum_private[j] < local_minimum[j] ) {
            local_minimum[j] = local_minimum_private[j];
            loc_node[j] = loc_node_private[j];
            loc_level[j] = loc_level_private[j];
          }
        }
      }
    }
  }
  std::vector< std::pair<T,int> > min_and_gidx_loc(nvar);
  std::vector< std::pair<T,int> > min_and_level_loc(nvar);
  std::vector< std::pair<T,int> > min_and_gidx_glb(nvar);
  std::vector< std::pair<T,int> > min_and_level_glb(nvar);
  const gidx_t* global_index = fs.nodes().global_index().data<gidx_t>();
  for( size_t j=0; j<nvar; ++j ) {
    gidx_t glb_idx = global_index[loc_node[j]];
    ASSERT( glb_idx < std::numeric_limits<int>::max() ); // pairs with 64bit integer for second not implemented
    min_and_gidx_loc[j] = std::make_pair(local_minimum[j],glb_idx);
    min_and_level_loc[j] = std::make_pair(local_minimum[j],loc_level[j]);
  }
  eckit::mpi::all_reduce(min_and_gidx_loc, min_and_gidx_glb, eckit::mpi::minloc());
  eckit::mpi::all_reduce(min_and_level_loc,min_and_level_glb,eckit::mpi::minloc());
  for( size_t j=0; j<nvar; ++j ) {
    min[j]     = min_and_gidx_glb[j].first;
    glb_idx[j] = min_and_gidx_glb[j].second;
    level[j]   = min_and_level_glb[j].second;
  }
}

template< typename T >
void minimum_and_location( const NodesFunctionSpace& fs, const Field& field, std::vector<T>& min, std::vector<gidx_t>& glb_idx, std::vector<size_t>& level )
{
  if( field.datatype() == DataType::kind<T>() ) {
    return dispatch_minimum_and_location(fs,field,min,glb_idx,level);
  }
  else
  {
    switch( field.datatype().kind() )
    {
      case DataType::KIND_INT32 : {
        std::vector<int> tmp;
        dispatch_minimum_and_location(fs,field,tmp,glb_idx,level);
        min.assign(tmp.begin(),tmp.end());
        return;
      }
      case DataType::KIND_INT64 : {
        std::vector<long> tmp;
        dispatch_minimum_and_location(fs,field,tmp,glb_idx,level);
        min.assign(tmp.begin(),tmp.end());
        return;
      }
      case DataType::KIND_REAL32 : {
        std::vector<float> tmp;
        dispatch_minimum_and_location(fs,field,tmp,glb_idx,level);
        min.assign(tmp.begin(),tmp.end());
        return;
      }
      case DataType::KIND_REAL64 : {
        std::vector<double> tmp;
        dispatch_minimum_and_location(fs,field,tmp,glb_idx,level);
        min.assign(tmp.begin(),tmp.end());
        return;
      }
      default: throw eckit::Exception("datatype not supported",Here());
    }
  }
}


template< typename T >
void dispatch_maximum_and_location( const NodesFunctionSpace& fs, const Field& field, std::vector<T>& max, std::vector<gidx_t>& glb_idx, std::vector<size_t>& level )
{
  ArrayView<T,3> arr = leveled_view<T>(field);
  size_t nvar = arr.shape(2);
  max.resize(nvar);
  glb_idx.resize(nvar);
  level.resize(nvar);
  std::vector<T> local_maximum(nvar,-std::numeric_limits<T>::max());
  std::vector<size_t> loc_node(nvar);
  std::vector<size_t> loc_level(nvar);
  atlas_omp_parallel
  {
    std::vector<T> local_maximum_private(nvar,-std::numeric_limits<T>::max());
    std::vector<size_t> loc_node_private(nvar);
    std::vector<size_t> loc_level_private(nvar);
    atlas_omp_for( size_t n=0; n<arr.shape(0); ++n ) {
      for( size_t l=0; l<arr.shape(1); ++l ) {
        for( size_t j=0; j<nvar; ++j ) {
          if( arr(n,l,j) > local_maximum_private[j] ) {
            local_maximum_private[j] = arr(n,l,j);
            loc_node_private[j] = n;
            loc_level_private[j] = l;
          }
        }
      }
    }
    atlas_omp_critical_ordered
    {
      for( size_t l=0; l<arr.shape(1); ++l ) {
        for( size_t j=0; j<nvar; ++j ) {
          if( local_maximum_private[j] > local_maximum[j] ) {
            local_maximum[j] = local_maximum_private[j];
            loc_node[j] = loc_node_private[j];
            loc_level[j] = loc_level_private[j];
          }
        }
      }
    }
  }
  std::vector< std::pair<T,int> > max_and_gidx_loc(nvar);
  std::vector< std::pair<T,int> > max_and_level_loc(nvar);
  std::vector< std::pair<T,int> > max_and_gidx_glb(nvar);
  std::vector< std::pair<T,int> > max_and_level_glb(nvar);
  const gidx_t* global_index = fs.nodes().global_index().data<gidx_t>();
  for( size_t j=0; j<nvar; ++j ) {
    gidx_t glb_idx = global_index[loc_node[j]];
    ASSERT( glb_idx < std::numeric_limits<int>::max() ); // pairs with 64bit integer for second not implemented
    max_and_gidx_loc[j] = std::make_pair(local_maximum[j],glb_idx);
    max_and_level_loc[j] = std::make_pair(local_maximum[j],loc_level[j]);
  }
  eckit::mpi::all_reduce(max_and_gidx_loc, max_and_gidx_glb, eckit::mpi::maxloc());
  eckit::mpi::all_reduce(max_and_level_loc,max_and_level_glb,eckit::mpi::maxloc());
  for( size_t j=0; j<nvar; ++j ) {
    max[j]     = max_and_gidx_glb[j].first;
    glb_idx[j] = max_and_gidx_glb[j].second;
    level[j]   = max_and_level_glb[j].second;
  }
}

template< typename T >
void maximum_and_location( const NodesFunctionSpace& fs, const Field& field, std::vector<T>& max, std::vector<gidx_t>& glb_idx, std::vector<size_t>& level )
{
  if( field.datatype() == DataType::kind<T>() ) {
    return dispatch_maximum_and_location(fs,field,max,glb_idx,level);
  }
  else
  {
    switch( field.datatype().kind() )
    {
      case DataType::KIND_INT32 : {
        std::vector<int> tmp;
        dispatch_maximum_and_location(fs,field,tmp,glb_idx,level);
        max.assign(tmp.begin(),tmp.end());
        return;
      }
      case DataType::KIND_INT64 : {
        std::vector<long> tmp;
        dispatch_maximum_and_location(fs,field,tmp,glb_idx,level);
        max.assign(tmp.begin(),tmp.end());
        return;
      }
      case DataType::KIND_REAL32 : {
        std::vector<float> tmp;
        dispatch_maximum_and_location(fs,field,tmp,glb_idx,level);
        max.assign(tmp.begin(),tmp.end());
        return;
      }
      case DataType::KIND_REAL64 : {
        std::vector<double> tmp;
        dispatch_maximum_and_location(fs,field,tmp,glb_idx,level);
        max.assign(tmp.begin(),tmp.end());
        return;
      }
      default: throw eckit::Exception("datatype not supported",Here());
    }
  }
}

template< typename T >
void minimum_and_location( const NodesFunctionSpace& fs, const Field& field, std::vector<T>& min, std::vector<gidx_t>& glb_idx)
{
  std::vector<size_t> level;
  minimum_and_location(fs,field,min,glb_idx,level);
}

template< typename T >
void maximum_and_location( const NodesFunctionSpace& fs, const Field& field, std::vector<T>& max, std::vector<gidx_t>& glb_idx)
{
  std::vector<size_t> level;
  maximum_and_location(fs,field,max,glb_idx,level);
}

template< typename T >
void minimum_and_location( const NodesFunctionSpace& fs, const Field& field, T& min, gidx_t& glb_idx, size_t& level)
{
  std::vector<T> minv;
  std::vector<gidx_t> gidxv;
  std::vector<size_t> levelv;
  minimum_and_location(fs,field,minv,gidxv,levelv);
  min = minv[0];
  glb_idx = gidxv[0];
  level = levelv[0];
}

template< typename T >
void maximum_and_location( const NodesFunctionSpace& fs, const Field& field, T& max, gidx_t& glb_idx, size_t& level)
{
  std::vector<T> maxv;
  std::vector<gidx_t> gidxv;
  std::vector<size_t> levelv;
  maximum_and_location(fs,field,maxv,gidxv,levelv);
  max = maxv[0];
  glb_idx = gidxv[0];
  level = levelv[0];
}

template< typename T >
void minimum_and_location( const NodesFunctionSpace& fs, const Field& field, T& min, gidx_t& glb_idx)
{
  size_t level;
  minimum_and_location(fs,field,min,glb_idx,level);
}

template< typename T >
void maximum_and_location( const NodesFunctionSpace& fs, const Field& field, T& max, gidx_t& glb_idx)
{
  size_t level;
  maximum_and_location(fs,field,max,glb_idx,level);
}

template< typename T >
void dispatch_minimum_and_location_per_level( const NodesFunctionSpace& fs, const Field& field, Field& min_field, Field& glb_idx_field )
{
  const ArrayView<T,3> arr = leveled_view<T>(field);
  ArrayShape shape;
  shape.reserve(field.rank()-1);
  for( size_t j=1; j<field.rank(); ++j )
    shape.push_back(field.shape(j));
  min_field.resize(shape);
  glb_idx_field.resize(shape);
  const size_t nvar = arr.shape(2);
  ArrayView<T,2> min( min_field.data<T>(), make_shape(min_field.shape(0),min_field.stride(0)) );
  min = std::numeric_limits<T>::max();
  ArrayView<gidx_t,2> glb_idx( glb_idx_field.data<gidx_t>(), make_shape(glb_idx_field.shape(0),glb_idx_field.stride(0)) );

  atlas_omp_parallel
  {
    Array<T> min_private(min.shape(0),min.shape(1));
    ArrayView<T> min_private_view(min_private); min_private_view = std::numeric_limits<T>::max();
    Array<gidx_t> glb_idx_private(glb_idx.shape(0),glb_idx.shape(1));
    atlas_omp_for( size_t n=0; n<arr.shape(0); ++n ) {
      for( size_t l=0; l<arr.shape(1); ++l ) {
        for( size_t j=0; j<nvar; ++j ) {
          if( arr(n,l,j) < min(l,j) ) {
            min_private(l,j) = arr(n,l,j);
            glb_idx_private(l,j) = n;
          }
        }
      }
    }
    atlas_omp_critical_ordered
    {
      for( size_t l=0; l<arr.shape(1); ++l ) {
        for( size_t j=0; j<nvar; ++j ) {
          if( min_private(l,j) < min(l,j) ) {
            min(l,j) = min_private(l,j);
            glb_idx(l,j) = glb_idx_private(l,j);
          }
        }
      }
    }
  }
  const size_t nlev = arr.shape(1);
  std::vector< std::pair<T,int> > min_and_gidx_loc(nlev*nvar);
  std::vector< std::pair<T,int> > min_and_gidx_glb(nlev*nvar);
  const gidx_t* global_index = fs.nodes().global_index().data<gidx_t>();
  atlas_omp_parallel_for( size_t l=0; l<nlev; ++l ) {
    for( size_t j=0; j<nvar; ++j ) {
      gidx_t gidx = global_index[glb_idx(l,j)];
      ASSERT( gidx < std::numeric_limits<int>::max() ); // pairs with 64bit integer for second not implemented
      min_and_gidx_loc[j+nvar*l] = std::make_pair(min(l,j),gidx);
    }
  }

  eckit::mpi::all_reduce(min_and_gidx_loc,min_and_gidx_glb,eckit::mpi::minloc());

  atlas_omp_parallel_for( size_t l=0; l<nlev; ++l ) {
    for( size_t j=0; j<nvar; ++j ) {
      min(l,j)     = min_and_gidx_glb[j+l*nvar].first;
      glb_idx(l,j) = min_and_gidx_glb[j+l*nvar].second;
    }
  }
}

void minimum_and_location_per_level( const NodesFunctionSpace& fs, const Field& field, Field& min, Field& glb_idx )
{
  if( field.datatype() != min.datatype() ) {
    throw eckit::Exception("Field and min are not of same datatype.",Here());
  }
  if( glb_idx.datatype() != DataType::kind<gidx_t>() ) {
    throw eckit::Exception("glb_idx Field is not of correct datatype",Here());
  }
  switch( field.datatype().kind() )
  {
    case DataType::KIND_INT32 :
      return dispatch_minimum_and_location_per_level<int>(fs,field,min,glb_idx);
    case DataType::KIND_INT64 :
      return dispatch_minimum_and_location_per_level<long>(fs,field,min,glb_idx);
    case DataType::KIND_REAL32 :
      return dispatch_minimum_and_location_per_level<float>(fs,field,min,glb_idx);
    case DataType::KIND_REAL64 :
      return dispatch_minimum_and_location_per_level<double>(fs,field,min,glb_idx);
    default: throw eckit::Exception("datatype not supported",Here());
  }
}


template< typename T >
void dispatch_maximum_and_location_per_level( const NodesFunctionSpace& fs, const Field& field, Field& max_field, Field& glb_idx_field )
{
  const ArrayView<T,3> arr = leveled_view<T>(field);
  ArrayShape shape;
  shape.reserve(field.rank()-1);
  for( size_t j=1; j<field.rank(); ++j )
    shape.push_back(field.shape(j));
  max_field.resize(shape);
  glb_idx_field.resize(shape);
  const size_t nvar = arr.shape(2);
  ArrayView<T,2> max( max_field.data<T>(), make_shape(max_field.shape(0),max_field.stride(0)) );
  max = -std::numeric_limits<T>::max();
  ArrayView<gidx_t,2> glb_idx( glb_idx_field.data<gidx_t>(), make_shape(glb_idx_field.shape(0),glb_idx_field.stride(0)) );

  atlas_omp_parallel
  {
    Array<T> max_private(max.shape(0),max.shape(1));
    ArrayView<T> max_private_view(max_private); max_private_view = -std::numeric_limits<T>::max();
    Array<gidx_t> glb_idx_private(glb_idx.shape(0),glb_idx.shape(1));
    atlas_omp_for( size_t n=0; n<arr.shape(0); ++n ) {
      for( size_t l=0; l<arr.shape(1); ++l ) {
        for( size_t j=0; j<nvar; ++j ) {
          if( arr(n,l,j) > max(l,j) ) {
            max_private(l,j) = arr(n,l,j);
            glb_idx(l,j) = n;
          }
        }
      }
    }
    atlas_omp_critical_ordered
    {
      for( size_t l=0; l<arr.shape(1); ++l ) {
        for( size_t j=0; j<nvar; ++j ) {
          if( max_private(l,j) > max(l,j) ) {
            max(l,j) = max_private(l,j);
            glb_idx(l,j) = glb_idx_private(l,j);
          }
        }
      }
    }
  }

  const size_t nlev = arr.shape(1);
  std::vector< std::pair<T,int> > max_and_gidx_loc(nlev*nvar);
  std::vector< std::pair<T,int> > max_and_gidx_glb(nlev*nvar);
  const gidx_t* global_index = fs.nodes().global_index().data<gidx_t>();
  atlas_omp_parallel_for( size_t l=0; l<nlev; ++l ) {
    for( size_t j=0; j<nvar; ++j ) {
      gidx_t gidx = global_index[glb_idx(l,j)];
      ASSERT( gidx < std::numeric_limits<int>::max() ); // pairs with 64bit integer for second not implemented
      max_and_gidx_loc[j+nvar*l] = std::make_pair(max(l,j),gidx);
    }
  }

  eckit::mpi::all_reduce(max_and_gidx_loc,max_and_gidx_glb,eckit::mpi::maxloc());

  atlas_omp_parallel_for( size_t l=0; l<nlev; ++l ) {
    for( size_t j=0; j<nvar; ++j ) {
      max(l,j)     = max_and_gidx_glb[j+l*nvar].first;
      glb_idx(l,j) = max_and_gidx_glb[j+l*nvar].second;
    }
  }
}

void maximum_and_location_per_level( const NodesFunctionSpace& fs, const Field& field, Field& max, Field& glb_idx )
{
  if( field.datatype() != max.datatype() ) {
    throw eckit::Exception("Field and max are not of same datatype.",Here());
  }
  if( glb_idx.datatype() != DataType::kind<gidx_t>() ) {
    throw eckit::Exception("glb_idx Field is not of correct datatype",Here());
  }
  switch( field.datatype().kind() )
  {
    case DataType::KIND_INT32 :
      return dispatch_maximum_and_location_per_level<int>(fs,field,max,glb_idx);
    case DataType::KIND_INT64 :
      return dispatch_maximum_and_location_per_level<long>(fs,field,max,glb_idx);
    case DataType::KIND_REAL32 :
      return dispatch_maximum_and_location_per_level<float>(fs,field,max,glb_idx);
    case DataType::KIND_REAL64 :
      return dispatch_maximum_and_location_per_level<double>(fs,field,max,glb_idx);
    default: throw eckit::Exception("datatype not supported",Here());
  }
}


template< typename T >
void mean( const NodesFunctionSpace& fs, const Field& field, T& result, size_t& N )
{
  sum(fs,field,result,N);
  result /= static_cast<double>(N);
}

template< typename T >
void mean( const NodesFunctionSpace& fs, const Field& field, std::vector<T>& result, size_t& N )
{
  sum(fs,field,result,N);
  for( size_t j=0; j<result.size(); ++j ) {
    result[j] /= static_cast<double>(N);
  }
}

template< typename T >
void dispatch_mean_per_level( const NodesFunctionSpace& fs, const Field& field, Field& mean, size_t& N )
{
  dispatch_sum_per_level<T>(fs,field,mean,N);
  T* rawdata = mean.data<T>();
  for( size_t j=0; j<mean.size(); ++j ) {
    rawdata[j] /= static_cast<double>(N);
  }
}


void mean_per_level( const NodesFunctionSpace& fs, const Field& field, Field& mean, size_t& N )
{
  if( field.datatype() != mean.datatype() ) {
    throw eckit::Exception("Field and sum are not of same datatype.",Here());
  }
  switch( field.datatype().kind() )
  {
    case DataType::KIND_INT32 :
      return dispatch_mean_per_level<int>(fs,field,mean,N);
    case DataType::KIND_INT64 :
      return dispatch_mean_per_level<long>(fs,field,mean,N);
    case DataType::KIND_REAL32 :
      return dispatch_mean_per_level<float>(fs,field,mean,N);
    case DataType::KIND_REAL64 :
      return dispatch_mean_per_level<double>(fs,field,mean,N);
    default: throw eckit::Exception("datatype not supported",Here());
  }
}

template< typename T >
void mean_and_standard_deviation( const NodesFunctionSpace& fs, const Field& field, T& mu, T& sigma, size_t& N )
{
  mean(fs,field,mu,N);
  Field::Ptr squared_diff_field( fs.createField(field) );
  ArrayView<T,2> squared_diff = leveled_scalar_view<T>( *squared_diff_field );
  ArrayView<T,2> values = leveled_scalar_view<T>( field );

  atlas_omp_parallel_for( size_t n=0; n<values.shape(0); ++n ) {
    for( size_t l=0; l<values.shape(1); ++l ) {
      squared_diff(n,l) = sqr( values(n,l) - mu );
    }
  }
  mean(fs,*squared_diff_field,sigma,N);
  sigma = std::sqrt(sigma);
}

template< typename T >
void mean_and_standard_deviation( const NodesFunctionSpace& fs, const Field& field, std::vector<T>& mu, std::vector<T>& sigma, size_t& N )
{
  mean(fs,field,mu,N);
  Field::Ptr squared_diff_field( fs.createField(field) );
  ArrayView<T,3> squared_diff = leveled_view<T>( *squared_diff_field );
  ArrayView<T,3> values = leveled_view<T>( field );

  atlas_omp_parallel_for( size_t n=0; n<values.shape(0); ++n ) {
    for( size_t l=0; l<values.shape(1); ++l ) {
      for( size_t j=0; j<values.shape(2); ++j ) {
        squared_diff(n,l,j) = sqr( values(n,l,j) - mu[j] );
      }
    }
  }
  mean(fs,*squared_diff_field,sigma,N);
  for( size_t j=0; j<sigma.size(); ++j ) {
    sigma[j] = std::sqrt(sigma[j]);
  }
}

template< typename T >
void dispatch_mean_and_standard_deviation_per_level( const NodesFunctionSpace& fs, const Field& field, Field& mean, Field& stddev, size_t& N )
{
  dispatch_mean_per_level<T>(fs,field,mean,N);
  Field::Ptr squared_diff_field( fs.createField(field) );
  ArrayView<T,3> squared_diff = leveled_view<T>( *squared_diff_field );
  ArrayView<T,3> values = leveled_view<T>( field );
  ArrayView<T,2> mu( mean.data<T>(), make_shape(values.shape(1),values.shape(2)) );

  atlas_omp_parallel_for( size_t n=0; n<values.shape(0); ++n ) {
    for( size_t l=0; l<values.shape(1); ++l ) {
      for( size_t j=0; j<values.shape(2); ++j ) {
        squared_diff(n,l,j) = sqr( values(n,l,j) - mu(l,j) );
      }
    }
  }
  dispatch_mean_per_level<T>(fs,*squared_diff_field,stddev,N);
  T* sigma = stddev.data<T>();
  atlas_omp_parallel_for( size_t j=0; j<stddev.size(); ++j ) {
    sigma[j] = std::sqrt(sigma[j]);
  }
}


void mean_and_standard_deviation_per_level( const NodesFunctionSpace& fs, const Field& field, Field& mean, Field& stddev, size_t& N )
{
  if( field.datatype() != mean.datatype() ) {
    throw eckit::Exception("Field and mean are not of same datatype.",Here());
  }
  if( field.datatype() != stddev.datatype() ) {
    throw eckit::Exception("Field and stddev are not of same datatype.",Here());
  }
  switch( field.datatype().kind() )
  {
    case DataType::KIND_INT32 :
      return dispatch_mean_and_standard_deviation_per_level<int>(fs,field,mean,stddev,N);
    case DataType::KIND_INT64 :
      return dispatch_mean_and_standard_deviation_per_level<long>(fs,field,mean,stddev,N);
    case DataType::KIND_REAL32 :
      return dispatch_mean_and_standard_deviation_per_level<float>(fs,field,mean,stddev,N);
    case DataType::KIND_REAL64 :
      return dispatch_mean_and_standard_deviation_per_level<double>(fs,field,mean,stddev,N);
    default: throw eckit::Exception("datatype not supported",Here());
  }
}


} // end collectives implementation

template<> void NodesFunctionSpace::sum( const Field& field, int&    result, size_t& N ) const { return detail::sum(*this,field,result,N); }
template<> void NodesFunctionSpace::sum( const Field& field, long&   result, size_t& N ) const { return detail::sum(*this,field,result,N); }
template<> void NodesFunctionSpace::sum( const Field& field, float&  result, size_t& N ) const { return detail::sum(*this,field,result,N); }
template<> void NodesFunctionSpace::sum( const Field& field, double& result, size_t& N ) const { return detail::sum(*this,field,result,N); }

template<> void NodesFunctionSpace::sum( const Field& field, std::vector<int>&    result, size_t& N ) const { return detail::sum(*this,field,result,N); }
template<> void NodesFunctionSpace::sum( const Field& field, std::vector<long>&   result, size_t& N ) const { return detail::sum(*this,field,result,N); }
template<> void NodesFunctionSpace::sum( const Field& field, std::vector<float>&  result, size_t& N ) const { return detail::sum(*this,field,result,N); }
template<> void NodesFunctionSpace::sum( const Field& field, std::vector<double>& result, size_t& N ) const { return detail::sum(*this,field,result,N); }


template<> void NodesFunctionSpace::orderIndependentSum( const Field& field, int&    result, size_t& N ) const { return detail::order_independent_sum(*this,field,result,N); }
template<> void NodesFunctionSpace::orderIndependentSum( const Field& field, long&   result, size_t& N ) const { return detail::order_independent_sum(*this,field,result,N); }
template<> void NodesFunctionSpace::orderIndependentSum( const Field& field, float&  result, size_t& N ) const { return detail::order_independent_sum(*this,field,result,N); }
template<> void NodesFunctionSpace::orderIndependentSum( const Field& field, double& result, size_t& N ) const { return detail::order_independent_sum(*this,field,result,N); }

template<> void NodesFunctionSpace::orderIndependentSum( const Field& field, std::vector<int>&    result, size_t& N ) const { return detail::order_independent_sum(*this,field,result,N); }
template<> void NodesFunctionSpace::orderIndependentSum( const Field& field, std::vector<long>&   result, size_t& N ) const { return detail::order_independent_sum(*this,field,result,N); }
template<> void NodesFunctionSpace::orderIndependentSum( const Field& field, std::vector<float>&  result, size_t& N ) const { return detail::order_independent_sum(*this,field,result,N); }
template<> void NodesFunctionSpace::orderIndependentSum( const Field& field, std::vector<double>& result, size_t& N ) const { return detail::order_independent_sum(*this,field,result,N); }

template<> void NodesFunctionSpace::minimum( const Field& field, int&    min ) const { return detail::minimum(*this,field,min); }
template<> void NodesFunctionSpace::minimum( const Field& field, long&   min ) const { return detail::minimum(*this,field,min); }
template<> void NodesFunctionSpace::minimum( const Field& field, float&  min ) const { return detail::minimum(*this,field,min); }
template<> void NodesFunctionSpace::minimum( const Field& field, double& min ) const { return detail::minimum(*this,field,min); }

template<> void NodesFunctionSpace::maximum( const Field& field, int&    max ) const { return detail::maximum(*this,field,max); }
template<> void NodesFunctionSpace::maximum( const Field& field, long&   max ) const { return detail::maximum(*this,field,max); }
template<> void NodesFunctionSpace::maximum( const Field& field, float&  max ) const { return detail::maximum(*this,field,max); }
template<> void NodesFunctionSpace::maximum( const Field& field, double& max ) const { return detail::maximum(*this,field,max); }

template<> void NodesFunctionSpace::minimum( const Field& field, std::vector<int>&    min ) const { return detail::minimum(*this,field,min); }
template<> void NodesFunctionSpace::minimum( const Field& field, std::vector<long>&   min ) const { return detail::minimum(*this,field,min); }
template<> void NodesFunctionSpace::minimum( const Field& field, std::vector<float>&  min ) const { return detail::minimum(*this,field,min); }
template<> void NodesFunctionSpace::minimum( const Field& field, std::vector<double>& min ) const { return detail::minimum(*this,field,min); }

template<> void NodesFunctionSpace::maximum( const Field& field, std::vector<int>&    max ) const { return detail::maximum(*this,field,max); }
template<> void NodesFunctionSpace::maximum( const Field& field, std::vector<long>&   max ) const { return detail::maximum(*this,field,max); }
template<> void NodesFunctionSpace::maximum( const Field& field, std::vector<float>&  max ) const { return detail::maximum(*this,field,max); }
template<> void NodesFunctionSpace::maximum( const Field& field, std::vector<double>& max ) const { return detail::maximum(*this,field,max); }

template<> void NodesFunctionSpace::maximumAndLocation( const Field& field, int&    max, gidx_t& glb_idx ) const { return detail::maximum_and_location(*this,field,max,glb_idx); }
template<> void NodesFunctionSpace::maximumAndLocation( const Field& field, long&   max, gidx_t& glb_idx ) const { return detail::maximum_and_location(*this,field,max,glb_idx); }
template<> void NodesFunctionSpace::maximumAndLocation( const Field& field, float&  max, gidx_t& glb_idx ) const { return detail::maximum_and_location(*this,field,max,glb_idx); }
template<> void NodesFunctionSpace::maximumAndLocation( const Field& field, double& max, gidx_t& glb_idx ) const { return detail::maximum_and_location(*this,field,max,glb_idx); }

template<> void NodesFunctionSpace::minimumAndLocation( const Field& field, int&    min, gidx_t& glb_idx ) const { return detail::minimum_and_location(*this,field,min,glb_idx); }
template<> void NodesFunctionSpace::minimumAndLocation( const Field& field, long&   min, gidx_t& glb_idx ) const { return detail::minimum_and_location(*this,field,min,glb_idx); }
template<> void NodesFunctionSpace::minimumAndLocation( const Field& field, float&  min, gidx_t& glb_idx ) const { return detail::minimum_and_location(*this,field,min,glb_idx); }
template<> void NodesFunctionSpace::minimumAndLocation( const Field& field, double& min, gidx_t& glb_idx ) const { return detail::minimum_and_location(*this,field,min,glb_idx); }

template<> void NodesFunctionSpace::maximumAndLocation( const Field& field, std::vector<int>&    max, std::vector<gidx_t>& glb_idx ) const { return detail::maximum_and_location(*this,field,max,glb_idx); }
template<> void NodesFunctionSpace::maximumAndLocation( const Field& field, std::vector<long>&   max, std::vector<gidx_t>& glb_idx ) const { return detail::maximum_and_location(*this,field,max,glb_idx); }
template<> void NodesFunctionSpace::maximumAndLocation( const Field& field, std::vector<float>&  max, std::vector<gidx_t>& glb_idx ) const { return detail::maximum_and_location(*this,field,max,glb_idx); }
template<> void NodesFunctionSpace::maximumAndLocation( const Field& field, std::vector<double>& max, std::vector<gidx_t>& glb_idx ) const { return detail::maximum_and_location(*this,field,max,glb_idx); }

template<> void NodesFunctionSpace::minimumAndLocation( const Field& field, std::vector<int>&    min, std::vector<gidx_t>& glb_idx ) const { return detail::minimum_and_location(*this,field,min,glb_idx); }
template<> void NodesFunctionSpace::minimumAndLocation( const Field& field, std::vector<long>&   min, std::vector<gidx_t>& glb_idx ) const { return detail::minimum_and_location(*this,field,min,glb_idx); }
template<> void NodesFunctionSpace::minimumAndLocation( const Field& field, std::vector<float>&  min, std::vector<gidx_t>& glb_idx ) const { return detail::minimum_and_location(*this,field,min,glb_idx); }
template<> void NodesFunctionSpace::minimumAndLocation( const Field& field, std::vector<double>& min, std::vector<gidx_t>& glb_idx ) const { return detail::minimum_and_location(*this,field,min,glb_idx); }


template<> void NodesFunctionSpace::minimumAndLocation( const Field& field, int&    min, gidx_t& glb_idx, size_t& level ) const { return detail::minimum_and_location(*this,field,min,glb_idx,level); }
template<> void NodesFunctionSpace::minimumAndLocation( const Field& field, long&   min, gidx_t& glb_idx, size_t& level ) const { return detail::minimum_and_location(*this,field,min,glb_idx,level); }
template<> void NodesFunctionSpace::minimumAndLocation( const Field& field, float&  min, gidx_t& glb_idx, size_t& level ) const { return detail::minimum_and_location(*this,field,min,glb_idx,level); }
template<> void NodesFunctionSpace::minimumAndLocation( const Field& field, double& min, gidx_t& glb_idx, size_t& level ) const { return detail::minimum_and_location(*this,field,min,glb_idx,level); }

template<> void NodesFunctionSpace::maximumAndLocation( const Field& field, int&    max, gidx_t& glb_idx, size_t& level ) const { return detail::maximum_and_location(*this,field,max,glb_idx,level); }
template<> void NodesFunctionSpace::maximumAndLocation( const Field& field, long&   max, gidx_t& glb_idx, size_t& level ) const { return detail::maximum_and_location(*this,field,max,glb_idx,level); }
template<> void NodesFunctionSpace::maximumAndLocation( const Field& field, float&  max, gidx_t& glb_idx, size_t& level ) const { return detail::maximum_and_location(*this,field,max,glb_idx,level); }
template<> void NodesFunctionSpace::maximumAndLocation( const Field& field, double& max, gidx_t& glb_idx, size_t& level ) const { return detail::maximum_and_location(*this,field,max,glb_idx,level); }

template<> void NodesFunctionSpace::minimumAndLocation( const Field& field, std::vector<int>&    min, std::vector<gidx_t>& glb_idx, std::vector<size_t>& level ) const { return detail::minimum_and_location(*this,field,min,glb_idx,level); }
template<> void NodesFunctionSpace::minimumAndLocation( const Field& field, std::vector<long>&   min, std::vector<gidx_t>& glb_idx, std::vector<size_t>& level ) const { return detail::minimum_and_location(*this,field,min,glb_idx,level); }
template<> void NodesFunctionSpace::minimumAndLocation( const Field& field, std::vector<float>&  min, std::vector<gidx_t>& glb_idx, std::vector<size_t>& level ) const { return detail::minimum_and_location(*this,field,min,glb_idx,level); }
template<> void NodesFunctionSpace::minimumAndLocation( const Field& field, std::vector<double>& min, std::vector<gidx_t>& glb_idx, std::vector<size_t>& level ) const { return detail::minimum_and_location(*this,field,min,glb_idx,level); }

template<> void NodesFunctionSpace::maximumAndLocation( const Field& field, std::vector<int>&    max, std::vector<gidx_t>& glb_idx, std::vector<size_t>& level ) const { return detail::maximum_and_location(*this,field,max,glb_idx,level); }
template<> void NodesFunctionSpace::maximumAndLocation( const Field& field, std::vector<long>&   max, std::vector<gidx_t>& glb_idx, std::vector<size_t>& level ) const { return detail::maximum_and_location(*this,field,max,glb_idx,level); }
template<> void NodesFunctionSpace::maximumAndLocation( const Field& field, std::vector<float>&  max, std::vector<gidx_t>& glb_idx, std::vector<size_t>& level ) const { return detail::maximum_and_location(*this,field,max,glb_idx,level); }
template<> void NodesFunctionSpace::maximumAndLocation( const Field& field, std::vector<double>& max, std::vector<gidx_t>& glb_idx, std::vector<size_t>& level ) const { return detail::maximum_and_location(*this,field,max,glb_idx,level); }



template<> void NodesFunctionSpace::mean( const Field& field, float&  result, size_t& N ) const { return detail::mean(*this,field,result,N); }
template<> void NodesFunctionSpace::mean( const Field& field, double& result, size_t& N ) const { return detail::mean(*this,field,result,N); }

template<> void NodesFunctionSpace::mean( const Field& field, std::vector<float>&  result, size_t& N ) const { return detail::mean(*this,field,result,N); }
template<> void NodesFunctionSpace::mean( const Field& field, std::vector<double>& result, size_t& N ) const { return detail::mean(*this,field,result,N); }

template<> void NodesFunctionSpace::meanAndStandardDeviation( const Field& field, float&  mu, float&  sigma, size_t& N ) const { return detail::mean_and_standard_deviation(*this,field,mu,sigma,N); }
template<> void NodesFunctionSpace::meanAndStandardDeviation( const Field& field, double& mu, double& sigma, size_t& N ) const { return detail::mean_and_standard_deviation(*this,field,mu,sigma,N); }

template<> void NodesFunctionSpace::meanAndStandardDeviation( const Field& field, std::vector<float>&  mu, std::vector<float>&  sigma, size_t& N ) const { return detail::mean_and_standard_deviation(*this,field,mu,sigma,N); }
template<> void NodesFunctionSpace::meanAndStandardDeviation( const Field& field, std::vector<double>& mu, std::vector<double>& sigma, size_t& N ) const { return detail::mean_and_standard_deviation(*this,field,mu,sigma,N); }


void NodesFunctionSpace::sumPerLevel(const Field &field, Field &sum, size_t &N) const { return detail::sum_per_level(*this,field,sum,N); }

void NodesFunctionSpace::orderIndependentSumPerLevel(const Field &field, Field &sum, size_t &N) const { return detail::order_independent_sum_per_level(*this,field,sum,N); }

void NodesFunctionSpace::minimumPerLevel(const Field &field, Field &min) const { return detail::minimum_per_level(*this,field,min); }

void NodesFunctionSpace::maximumPerLevel(const Field &field, Field &max) const { return detail::maximum_per_level(*this,field,max);}

void NodesFunctionSpace::minimumAndLocationPerLevel(const Field &field, Field &min, Field &glb_idx) const { detail::minimum_and_location_per_level(*this,field,min,glb_idx); }

void NodesFunctionSpace::maximumAndLocationPerLevel(const Field &field, Field &max, Field &glb_idx) const { detail::maximum_and_location_per_level(*this,field,max,glb_idx); }

void NodesFunctionSpace::meanPerLevel(const Field &field, Field &mean, size_t &N) const { return detail::mean_per_level(*this,field,mean,N); }

void NodesFunctionSpace::meanAndStandardDeviationPerLevel(const Field &field, Field &mean, Field &stddev, size_t &N) const { return detail::mean_and_standard_deviation_per_level(*this,field,mean,stddev,N); }


} // namespace functionspace
} // namespace atlas

