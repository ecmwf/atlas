/*
 * (C) Copyright 1996-2017 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include <algorithm>
#include <cmath>
#include <cstdarg>
#include <limits>

#include "eckit/os/BackTrace.h"
#include "eckit/utils/MD5.h"

#include "atlas/library/config.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/mesh/Nodes.h"
#include "atlas/mesh/actions/BuildParallelFields.h"
#include "atlas/mesh/actions/BuildHalo.h"
#include "atlas/mesh/actions/BuildPeriodicBoundaries.h"
#include "atlas/functionspace/NodeColumns.h"
#include "atlas/mesh/IsGhostNode.h"
#include "atlas/parallel/HaloExchange.h"
#include "atlas/parallel/GatherScatter.h"
#include "atlas/parallel/Checksum.h"
#include "atlas/parallel/omp/omp.h"
#include "atlas/runtime/ErrorHandling.h"
#include "atlas/runtime/Log.h"
#undef atlas_omp_critical_ordered
#define atlas_omp_critical_ordered atlas_omp_critical
#include "atlas/array/ArrayView.h"

#ifdef ATLAS_HAVE_FORTRAN
#define REMOTE_IDX_BASE 1
#else
#define REMOTE_IDX_BASE 0
#endif

namespace atlas {
namespace functionspace {
namespace detail {

namespace {

template <typename T>
array::LocalView<T,3> make_leveled_view(const field::Field &field)
{
  if( field.has_levels() )
    return array::LocalView<T,3> ( array::make_storageview<T>(field).data(),
                                   array::make_shape(field.shape(0),field.shape(1),field.stride(1)) );
  else
    return array::LocalView<T,3> ( array::make_storageview<T>(field).data(),
                                   array::make_shape(field.shape(0),1,field.stride(0)) );
}

template <typename T>
array::LocalView<T,2> surface_view(const field::Field &field)
{
  return array::LocalView<T,2> ( array::make_storageview<T>(field).data(),
                                 array::make_shape(field.shape(0),field.stride(0)) );
}

template <typename T>
array::LocalView<T,2> make_leveled_scalar_view(const field::Field &field)
{
  if( field.has_levels() )
    return array::LocalView<T,2> ( array::make_storageview<T>(field).data(),
                                   array::make_shape(field.shape(0),field.shape(1)) );
  else
    return array::LocalView<T,2> ( array::make_storageview<T>(field).data(),
                                   array::make_shape(field.shape(0),1) );
}

template <typename T>
array::LocalView<T,1> make_surface_scalar_view(const field::Field &field)
{
  return array::LocalView<T,1> ( array::make_storageview<T>(field).data(),
                                 array::make_shape(field.size()) );
}


}

NodeColumns::NodeColumns( mesh::Mesh& mesh ) :
    mesh_(mesh),
    nodes_(mesh_.nodes()),
    halo_(0),
    nb_nodes_(nodes_.size()),
    nb_nodes_global_(0) {
    constructor();
}

NodeColumns::NodeColumns( mesh::Mesh& mesh, const mesh::Halo &halo, const eckit::Parametrisation &params ) :
    mesh_(mesh),
    nodes_(mesh_.nodes()),
    halo_(halo),
    nb_nodes_(nodes_.size()),
    nb_nodes_global_(0) {
    constructor();
}

NodeColumns::NodeColumns(mesh::Mesh& mesh, const mesh::Halo &halo) :
    mesh_(mesh),
    nodes_(mesh_.nodes()),
    halo_(halo),
    nb_nodes_(nodes_.size()),
    nb_nodes_global_(0) {
    constructor();
}


void NodeColumns::constructor()
{
  mesh::actions::build_nodes_parallel_fields( mesh_.nodes() );
  mesh::actions::build_periodic_boundaries(mesh_);

  gather_scatter_.reset(new parallel::GatherScatter());
  halo_exchange_.reset(new parallel::HaloExchange());
  checksum_.reset(new parallel::Checksum());

  if( halo_.size() > 0)
  {
    mesh::actions::build_halo(mesh_,halo_.size());
    mesh::actions::renumber_nodes_glb_idx(mesh_.nodes());
    std::stringstream ss;
    ss << "nb_nodes_including_halo["<<halo_.size()<<"]";
    mesh_.metadata().get(ss.str(),nb_nodes_);

  }
  if( !nb_nodes_ ) {
    std::stringstream ss;
    ss << "nb_nodes_including_halo["<<halo_.size()<<"]";
    if( ! mesh_.metadata().get(ss.str(),nb_nodes_) ) {
      nb_nodes_ = mesh_.nodes().size();
    }
  }

  field::Field& part = mesh_.nodes().partition();
  field::Field& ridx = mesh_.nodes().remote_index();
  halo_exchange_->setup( array::make_view<int,1>(part).data(),
                         array::make_view<int,1>(ridx).data(),
                         REMOTE_IDX_BASE,nb_nodes_);

  {
    // Create new gather_scatter
    gather_scatter_.reset( new parallel::GatherScatter() );

    field::Field& ridx = mesh_.nodes().remote_index();
    field::Field& part = mesh_.nodes().partition();
    field::Field& gidx = mesh_.nodes().global_index();

    mesh::IsGhostNode is_ghost(mesh_.nodes());
    std::vector<int> mask(mesh_.nodes().size());
    const size_t npts = mask.size();
    atlas_omp_parallel_for( size_t n=0; n<npts; ++n ) {
      mask[n] = is_ghost(n) ? 1 : 0;

      // --> This would add periodic west-bc to the gather, but means that global-sums, means, etc are computed wrong
      //if( mask[j] == 1 && internals::Topology::check(flags(j),internals::Topology::BC) ) {
      //  mask[j] = 0;
      //}
    }
    gather_scatter_->setup(array::make_view<int,1>(part).data(),
                           array::make_view<int,1>(ridx).data(),
                           REMOTE_IDX_BASE,
                           array::make_view<gidx_t,1>(gidx).data(),
                           mask.data(),nb_nodes_);
  }

  {
    // Create new checksum
    checksum_.reset( new parallel::Checksum() );

    field::Field& ridx = mesh_.nodes().remote_index();
    field::Field& part = mesh_.nodes().partition();
    field::Field& gidx = mesh_.nodes().global_index();
    mesh::IsGhostNode is_ghost(mesh_.nodes());
    std::vector<int> mask(mesh_.nodes().size());
    const size_t npts = mask.size();
    atlas_omp_parallel_for( size_t n=0; n<npts; ++n ) {
      mask[n] = is_ghost(n) ? 1 : 0;
    }
    checksum_->setup(array::make_view<int,1>(part).data(),
                           array::make_view<int,1>(ridx).data(),
                           REMOTE_IDX_BASE,
                           array::make_view<gidx_t,1>(gidx).data(),
                           mask.data(),nb_nodes_);
  }

  nb_nodes_global_ = gather_scatter_->glb_dof();
}

NodeColumns::~NodeColumns() {}

size_t NodeColumns::footprint() const {
  size_t size = sizeof(*this);
  // TODO
  return size;
}

size_t NodeColumns::nb_nodes() const
{
  return nb_nodes_;
}

size_t NodeColumns::nb_nodes_global() const
{
  return nb_nodes_global_;
}

std::string NodeColumns::halo_name() const
{
  std::stringstream ss; ss << "nodes_" << halo_.size();
  return ss.str();
}

std::string NodeColumns::gather_scatter_name() const
{
  return "nodes_gather_scatter";
}

std::string NodeColumns::checksum_name() const
{
  return "nodes_checksum";
}


size_t NodeColumns::config_nb_nodes(const eckit::Parametrisation& config) const
{
  size_t size = nb_nodes();
  bool global(false);
  if( config.get("global",global) )
  {
    if( global )
    {
      size_t owner(0);
      config.get("owner",owner);
      size = (parallel::mpi::comm().rank() == owner ? nb_nodes_global() : 0);
    }
  }
  return size;
}

namespace {
void set_field_metadata(const eckit::Parametrisation& config, field::Field& field)
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


array::DataType NodeColumns::config_datatype(const eckit::Parametrisation& config) const
{
  array::DataType::kind_t kind;
  if( ! config.get("datatype",kind) ) throw eckit::AssertionFailed("datatype missing",Here());
  return array::DataType(kind);
}

std::string NodeColumns::config_name(const eckit::Parametrisation& config) const
{
  std::string name;
  config.get("name",name);
  return name;
}

size_t NodeColumns::config_levels(const eckit::Parametrisation& config) const
{
  size_t levels(0);
  config.get("levels",levels);
  return levels;
}


field::Field NodeColumns::createField(
    const std::string& name,
    array::DataType datatype,
    const eckit::Parametrisation& config ) const
{
  size_t nb_nodes = config_nb_nodes(config);
  field::Field field = field::Field(name,datatype,array::make_shape(nb_nodes));
  field.set_functionspace(this);
  set_field_metadata(config,field);
  return field;
}

field::Field NodeColumns::createField(
    const std::string& name,
    array::DataType datatype,
    size_t levels,
    const eckit::Parametrisation& config) const
{
  size_t nb_nodes = config_nb_nodes(config);
  field::Field field = field::Field(name,datatype,array::make_shape(nb_nodes,levels));
  field.set_levels(levels);
  field.set_functionspace(this);
  set_field_metadata(config,field);
  return field;
}

field::Field NodeColumns::createField(
    const std::string& name,
    array::DataType datatype,
    const std::vector<size_t>& variables,
    const eckit::Parametrisation& config ) const
{
  size_t nb_nodes = config_nb_nodes(config);
  std::vector<size_t> shape(1,nb_nodes);
  for( size_t i=0; i<variables.size(); ++i ) shape.push_back(variables[i]);
  field::Field field = field::Field(name,datatype,shape);
  field.set_functionspace(this);
  set_field_metadata(config,field);
  return field;
}

field::Field NodeColumns::createField(
    const std::string& name,
    array::DataType datatype,
    size_t levels,
    const std::vector<size_t>& variables,
    const eckit::Parametrisation& config ) const
{
  size_t nb_nodes = config_nb_nodes(config);
  std::vector<size_t> shape(1,nb_nodes); shape.push_back(levels);
  for( size_t i=0; i<variables.size(); ++i ) shape.push_back(variables[i]);
  field::Field field = field::Field(name,datatype,shape);
  field.set_levels(levels);
  field.set_functionspace(this);
  set_field_metadata(config,field);
  return field;
}

field::Field NodeColumns::createField(
    const std::string& name,
    const field::Field& other,
    const eckit::Parametrisation& config ) const
{
  array::ArrayShape shape = other.shape();
  shape[0] = config_nb_nodes(config);
  field::Field field = field::Field(name,other.datatype(),shape);
  if( other.has_levels() )
    field.set_levels(field.shape(1));
  field.set_functionspace(this);
  set_field_metadata(config,field);
  return field;
}

field::Field NodeColumns::createField(
    const eckit::Parametrisation& config ) const
{
  size_t nb_nodes = config_nb_nodes(config);
  field::Field field = field::Field("",array::DataType::create<double>(),array::make_shape(nb_nodes));
  field.set_functionspace(this);
  set_field_metadata(config,field);
  return field;
}

namespace {
template<int RANK>
void dispatch_haloExchange( const field::Field& field, const parallel::HaloExchange& halo_exchange )
{
  if     ( field.datatype() == array::DataType::kind<int>() ) {
    array::ArrayView<int,RANK> view = array::make_view<int,RANK>(field);
    halo_exchange.execute( view );
  }
  else if( field.datatype() == array::DataType::kind<long>() ) {
    array::ArrayView<long,RANK> view = array::make_view<long,RANK>(field);
    halo_exchange.execute( view );
  }
  else if( field.datatype() == array::DataType::kind<float>() ) {
    array::ArrayView<float,RANK> view = array::make_view<float,RANK>(field);
    halo_exchange.execute( view );
  }
  else if( field.datatype() == array::DataType::kind<double>() ) {
    array::ArrayView<double,RANK> view = array::make_view<double,RANK>(field);
    halo_exchange.execute( view );
  }
  else throw eckit::Exception("datatype not supported",Here());
}
}

void NodeColumns::haloExchange( field::FieldSet& fieldset ) const
{
  for( size_t f=0; f<fieldset.size(); ++f ) {
    const field::Field& field = fieldset[f];
    switch( field.rank() ) {
      case 1:
        dispatch_haloExchange<1>(field,halo_exchange());
        break;
      case 2:
        dispatch_haloExchange<2>(field,halo_exchange());
        break;
      case 3:
        dispatch_haloExchange<3>(field,halo_exchange());
        break;
      case 4:
        dispatch_haloExchange<4>(field,halo_exchange());
        break;
      default:
        throw eckit::Exception("Rank not supported", Here());
        break;
    }
  }
}

void NodeColumns::haloExchange( field::Field& field ) const
{
  field::FieldSet fieldset;
  fieldset.add(field);
  haloExchange(fieldset);
}
const parallel::HaloExchange& NodeColumns::halo_exchange() const
{
  return *halo_exchange_;;
}


void NodeColumns::gather( const field::FieldSet& local_fieldset, field::FieldSet& global_fieldset ) const
{
  ASSERT(local_fieldset.size() == global_fieldset.size());

  for( size_t f=0; f<local_fieldset.size(); ++f ) {

    const field::Field& loc = local_fieldset[f];
    field::Field& glb = global_fieldset[f];
    const size_t nb_fields = 1;
    size_t root(0);
    glb.metadata().get("owner",root);

    if     ( loc.datatype() == array::DataType::kind<int>() ) {
      parallel::Field<int const> loc_field( array::make_storageview<int>(loc).data(),loc.stride(0));
      parallel::Field<int      > glb_field( array::make_storageview<int>(glb).data(),glb.stride(0));
      gather().gather( &loc_field, &glb_field, nb_fields, root );
    }
    else if( loc.datatype() == array::DataType::kind<long>() ) {
      parallel::Field<long const> loc_field( array::make_storageview<long>(loc).data(),loc.stride(0));
      parallel::Field<long      > glb_field( array::make_storageview<long>(glb).data(),glb.stride(0));
      gather().gather( &loc_field, &glb_field, nb_fields, root );
    }
    else if( loc.datatype() == array::DataType::kind<float>() ) {
      parallel::Field<float const> loc_field( array::make_storageview<float>(loc).data(),loc.stride(0));
      parallel::Field<float      > glb_field( array::make_storageview<float>(glb).data(),glb.stride(0));
      gather().gather( &loc_field, &glb_field, nb_fields, root );
    }
    else if( loc.datatype() == array::DataType::kind<double>() ) {
      parallel::Field<double const> loc_field( array::make_storageview<double>(loc).data(),loc.stride(0));
      parallel::Field<double      > glb_field( array::make_storageview<double>(glb).data(),glb.stride(0));
      gather().gather( &loc_field, &glb_field, nb_fields, root );
    }
    else throw eckit::Exception("datatype not supported",Here());
  }
}
void NodeColumns::gather( const field::Field& local, field::Field& global ) const
{
  field::FieldSet local_fields;
  field::FieldSet global_fields;
  local_fields.add(local);
  global_fields.add(global);
  gather(local_fields,global_fields);
}
const parallel::GatherScatter& NodeColumns::gather() const
{
  return *gather_scatter_;
}
const parallel::GatherScatter& NodeColumns::scatter() const
{
  return *gather_scatter_;
}


void NodeColumns::scatter( const field::FieldSet& global_fieldset, field::FieldSet& local_fieldset ) const
{
  ASSERT(local_fieldset.size() == global_fieldset.size());

  for( size_t f=0; f<local_fieldset.size(); ++f ) {

    const field::Field& glb = global_fieldset[f];
    field::Field& loc = local_fieldset[f];
    const size_t nb_fields = 1;
    size_t root(0);
    glb.metadata().get("owner",root);

    if     ( loc.datatype() == array::DataType::kind<int>() ) {
      parallel::Field<int const> glb_field( array::make_storageview<int>(glb).data(),glb.stride(0));
      parallel::Field<int      > loc_field( array::make_storageview<int>(loc).data(),loc.stride(0));
      scatter().scatter( &glb_field, &loc_field, nb_fields, root );
    }
    else if( loc.datatype() == array::DataType::kind<long>() ) {
      parallel::Field<long const> glb_field( array::make_storageview<long>(glb).data(),glb.stride(0));
      parallel::Field<long      > loc_field( array::make_storageview<long>(loc).data(),loc.stride(0));
      scatter().scatter( &glb_field, &loc_field, nb_fields, root );
    }
    else if( loc.datatype() == array::DataType::kind<float>() ) {
      parallel::Field<float const> glb_field( array::make_storageview<float>(glb).data(),glb.stride(0));
      parallel::Field<float      > loc_field( array::make_storageview<float>(loc).data(),loc.stride(0));
      scatter().scatter( &glb_field, &loc_field, nb_fields, root );
    }
    else if( loc.datatype() == array::DataType::kind<double>() ) {
      parallel::Field<double const> glb_field( array::make_storageview<double>(glb).data(),glb.stride(0));
      parallel::Field<double      > loc_field( array::make_storageview<double>(loc).data(),loc.stride(0));
      scatter().scatter( &glb_field, &loc_field, nb_fields, root );
    }
    else throw eckit::Exception("datatype not supported",Here());

    glb.metadata().broadcast(loc.metadata(),root);
    loc.metadata().set("global",false);
  }
}

void NodeColumns::scatter( const field::Field& global, field::Field& local ) const
{
  field::FieldSet global_fields;
  field::FieldSet local_fields;
  global_fields.add(global);
  local_fields.add(local);
  scatter(global_fields,local_fields);
}

namespace {
template <typename T>
std::string checksum_3d_field(const parallel::Checksum& checksum, const field::Field& field )
{
  array::LocalView<T,3> values = make_leveled_view<T>(field);
  array::ArrayT<T> surface_field( values.shape(0), values.shape(2) );
  array::ArrayView<T,2> surface = array::make_view<T,2>(surface_field);
  const size_t npts = values.shape(0);
  atlas_omp_for( size_t n=0; n<npts; ++n ) {
    for( size_t j=0; j<surface.shape(1); ++j )
    {
      surface(n,j) = 0.;
      for( size_t l=0; l<values.shape(1);++l )
        surface(n,j) += values(n,l,j);
    }
  }
  return checksum.execute( surface.data(), surface_field.stride(0) );
}
}

std::string NodeColumns::checksum( const field::FieldSet& fieldset ) const {
  eckit::MD5 md5;
  for( size_t f=0; f<fieldset.size(); ++f ) {
    const field::Field& field=fieldset[f];
    if     ( field.datatype() == array::DataType::kind<int>() )
      md5 << checksum_3d_field<int>(checksum(),field);
    else if( field.datatype() == array::DataType::kind<long>() )
      md5 << checksum_3d_field<long>(checksum(),field);
    else if( field.datatype() == array::DataType::kind<float>() )
      md5 << checksum_3d_field<float>(checksum(),field);
    else if( field.datatype() == array::DataType::kind<double>() )
      md5 << checksum_3d_field<double>(checksum(),field);
    else throw eckit::Exception("datatype not supported",Here());
  }
  return md5;
}
std::string NodeColumns::checksum( const field::Field& field ) const {
  field::FieldSet fieldset;
  fieldset.add(field);
  return checksum(fieldset);
}

const parallel::Checksum& NodeColumns::checksum() const
{
  return *checksum_;
}



//std::string NodesFunctionSpace::checksum( const field::FieldSet& fieldset ) const {
//  const parallel::Checksum& checksum = mesh_.checksum().get(checksum_name());

//  eckit::MD5 md5;
//  for( size_t f=0; f<fieldset.size(); ++f ) {
//    const field::Field& field=fieldset[f];
//    if     ( field.datatype() == array::DataType::kind<int>() )
//      md5 << checksum.execute( field.data<int>(), field.stride(0) );
//    else if( field.datatype() == array::DataType::kind<long>() )
//      md5 << checksum.execute( field.data<long>(), field.stride(0) );
//    else if( field.datatype() == array::DataType::kind<float>() )
//      md5 << checksum.execute( field.data<float>(), field.stride(0) );
//    else if( field.datatype() == array::DataType::kind<double>() )
//      md5 << checksum.execute( field.data<double>(), field.stride(0) );
//    else throw eckit::Exception("datatype not supported",Here());
//  }
//  return md5;
//}
//std::string NodesFunctionSpace::checksum( const field::Field& field ) const {
//  field::FieldSet fieldset;
//  fieldset.add(field);
//  return checksum(fieldset);
//}

namespace { inline double sqr(const double& val) { return val*val; } }

namespace detail { // Collectives implementation



template< typename T >
void dispatch_sum( const NodeColumns& fs, const field::Field& field, T& result, size_t& N )
{
  const mesh::IsGhostNode is_ghost(fs.nodes());
  const array::LocalView<T,2> arr = make_leveled_scalar_view<T>( field );
  T local_sum = 0;
  const size_t npts = std::min(arr.shape(0),fs.nb_nodes());
  atlas_omp_pragma( omp parallel for default(shared) reduction(+:local_sum) )
  for( size_t n=0; n<npts; ++n ) {
    if( ! is_ghost(n) ) {
      for( size_t l=0; l<arr.shape(1); ++l )
        local_sum += arr(n,l);
    }
  }
  parallel::mpi::comm().allReduce(local_sum, result, eckit::mpi::sum());
  N = fs.nb_nodes_global() * arr.shape(1);
}

template< typename T >
void sum( const NodeColumns& fs , const field::Field& field, T& result, size_t& N )
{
  if( field.datatype() == array::DataType::kind<T>() ) {
    return dispatch_sum(fs,field,result,N);
  }
  else
  {
    switch( field.datatype().kind() )
    {
      case array::DataType::KIND_INT32 : {
        int tmp;
        dispatch_sum(fs,field,tmp,N);
        result = tmp;
        return;
      }
      case array::DataType::KIND_INT64 : {
        long tmp;
        dispatch_sum(fs,field,tmp,N);
        result = tmp;
        return;
      }
      case array::DataType::KIND_REAL32 : {
        float tmp;
        dispatch_sum(fs,field,tmp,N);
        result = tmp;
        return;
      }
      case array::DataType::KIND_REAL64 : {
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
void dispatch_sum( const NodeColumns& fs, const field::Field& field, std::vector<T>& result, size_t& N )
{
  const array::LocalView<T,3> arr = make_leveled_view<T>(field);
  const mesh::IsGhostNode is_ghost(fs.nodes());
  const size_t nvar = arr.shape(2);
  std::vector<T> local_sum(nvar,0);
  result.resize(nvar);

  atlas_omp_parallel
  {
    std::vector<T> local_sum_private(nvar,0);
    const size_t npts = arr.shape(0);
    atlas_omp_for( size_t n=0; n<npts; ++n )
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

  parallel::mpi::comm().allReduce(local_sum, result, eckit::mpi::sum());

  N = fs.nb_nodes_global() * arr.shape(1);
}

template< typename T >
void sum( const NodeColumns& fs , const field::Field& field, std::vector<T>& result, size_t& N )
{
  if( field.datatype() == array::DataType::kind<T>() ) {
    return dispatch_sum(fs,field,result,N);
  }
  else
  {
    switch( field.datatype().kind() )
    {
      case array::DataType::KIND_INT32 : {
        std::vector<int> tmp;
        dispatch_sum(fs,field,tmp,N);
        result.assign(tmp.begin(),tmp.end());
        return;
      }
      case array::DataType::KIND_INT64 : {
        std::vector<long> tmp;
        dispatch_sum(fs,field,tmp,N);
        result.assign(tmp.begin(),tmp.end());
        return;
      }
      case array::DataType::KIND_REAL32 : {
        std::vector<float> tmp;
        dispatch_sum(fs,field,tmp,N);
        result.assign(tmp.begin(),tmp.end());
        return;
      }
      case array::DataType::KIND_REAL64 : {
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
void dispatch_sum_per_level( const NodeColumns& fs, const field::Field& field, field::Field& sum, size_t& N )
{
  mesh::IsGhostNode is_ghost(fs.nodes());

  array::ArrayShape shape;
  shape.reserve(field.rank()-1);
  for( size_t j=1; j<field.rank(); ++j )
    shape.push_back(field.shape(j));
  sum.resize(shape);

  const array::LocalView<T,3> arr = make_leveled_view<T>(field);

  array::LocalView<T,2> sum_per_level(
    array::make_storageview<T>(sum).data(),
    array::make_shape(sum.shape(0),sum.stride(0)) );

  for( size_t l=0; l<sum_per_level.shape(0); ++l ) {
    for( size_t j=0; j<sum_per_level.shape(1); ++j ) {
      sum_per_level(l,j) = 0;
    }
  }


  atlas_omp_parallel
  {
    array::ArrayT<T> sum_per_level_private(sum_per_level.shape(0),sum_per_level.shape(1));
    array::ArrayView<T,2> sum_per_level_private_view = array::make_view<T,2>(sum_per_level_private);

    for( size_t l=0; l<sum_per_level_private_view.shape(0); ++l ) {
      for( size_t j=0; j<sum_per_level_private_view.shape(1); ++j ) {
        sum_per_level_private_view(l,j) = 0;
      }
    }

    const size_t npts = arr.shape(0);
    atlas_omp_for( size_t n=0; n<npts; ++n )
    {
      if( ! is_ghost(n) ) {
        for( size_t l=0; l<arr.shape(1); ++l ) {
          for( size_t j=0; j<arr.shape(2); ++j ) {
            sum_per_level_private_view(l,j) += arr(n,l,j);
          }
        }
      }
    }
    atlas_omp_critical
    {
      for( size_t l=0; l<sum_per_level_private.shape(0); ++l ) {
        for( size_t j=0; j<sum_per_level_private.shape(1); ++j ) {
          sum_per_level(l,j) += sum_per_level_private_view(l,j);
        }
      }
    }
  }
  parallel::mpi::comm().allReduceInPlace(sum_per_level.data(), sum.size(), eckit::mpi::sum());
  N = fs.nb_nodes_global();
}

void sum_per_level( const NodeColumns& fs, const field::Field& field, field::Field& sum, size_t& N )
{
  if( field.datatype() != sum.datatype() ) {
    throw eckit::Exception("field::Field and sum are not of same datatype.",Here());
  }
  switch( field.datatype().kind() )
  {
    case array::DataType::KIND_INT32 :
      return dispatch_sum_per_level<int>(fs,field,sum,N);
    case array::DataType::KIND_INT64 :
      return dispatch_sum_per_level<long>(fs,field,sum,N);
    case array::DataType::KIND_REAL32 :
      return dispatch_sum_per_level<float>(fs,field,sum,N);
    case array::DataType::KIND_REAL64 :
      return dispatch_sum_per_level<double>(fs,field,sum,N);
    default: throw eckit::Exception("datatype not supported",Here());
  }
}


template< typename DATATYPE >
void dispatch_order_independent_sum_2d( const NodeColumns& fs , const field::Field& field, DATATYPE& result, size_t& N )
{
  size_t root = 0;
  field::Field global = fs.createField("global",field, field::global() );
  fs.gather(field,global);
  result = std::accumulate(array::make_storageview<DATATYPE>(global).data(),
                           array::make_storageview<DATATYPE>(global).data()+global.size(),0.);
  parallel::mpi::comm().broadcast(&result, 1, root);
  N = fs.nb_nodes_global();
}

template< typename T >
void dispatch_order_independent_sum( const NodeColumns& fs , const field::Field& field, T& result, size_t& N )
{
  if( field.has_levels() )
  {
    const array::LocalView<T,2> arr = make_leveled_scalar_view<T>(field);

    field::Field surface_field = fs.createField<T>("surface");
    array::LocalView<T,1> surface = make_surface_scalar_view<T>( surface_field );

    for( size_t n=0; n<arr.shape(0); ++n ) {
      for( size_t l=0; l<arr.shape(1); ++l ) {
          surface(n) += arr(n,l);
      }
    }
    dispatch_order_independent_sum_2d( fs, surface_field, result, N );
    N *= arr.shape(1);
  }
  else
  {
    dispatch_order_independent_sum_2d( fs, field, result, N );
  }
}

template< typename T >
void order_independent_sum( const NodeColumns& fs , const field::Field& field, T& result, size_t& N )
{
  if( field.datatype() == array::DataType::kind<T>() ) {
    return dispatch_order_independent_sum(fs,field,result,N);
  }
  else
  {
    switch( field.datatype().kind() )
    {
      case array::DataType::KIND_INT32 : {
        int tmp;
        dispatch_order_independent_sum(fs,field,tmp,N);
        result = tmp;
        return;
      }
      case array::DataType::KIND_INT64 : {
        long tmp;
        dispatch_order_independent_sum(fs,field,tmp,N);
        result = tmp;
        return;
      }
      case array::DataType::KIND_REAL32 : {
        float tmp;
        dispatch_order_independent_sum(fs,field,tmp,N);
        result = tmp;
        return;
      }
      case array::DataType::KIND_REAL64 : {
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
void dispatch_order_independent_sum_2d( const NodeColumns& fs, const field::Field& field, std::vector<DATATYPE>& result, size_t& N )
{
  size_t nvar = field.stride(0);
  result.resize(nvar);
  for( size_t j=0; j<nvar; ++j ) result[j] = 0.;
  field::Field global = fs.createField("global",field, field::global() );
  fs.gather(field,global);
  if( parallel::mpi::comm().rank() == 0 ) {
    const array::LocalView<DATATYPE,2> glb( array::make_storageview<DATATYPE>(global).data(),
                                            array::make_shape(global.shape(0),global.stride(0)) );
    for( size_t n=0; n<fs.nb_nodes_global(); ++n ) {
      for( size_t j=0; j<nvar; ++j ) {
        result[j] += glb(n,j);
      }
    }
  }
  size_t root = global.metadata().get<size_t>("owner");
  parallel::mpi::comm().broadcast(result,root);
  N = fs.nb_nodes_global();
}

template< typename T >
void dispatch_order_independent_sum( const NodeColumns& fs, const field::Field& field, std::vector<T>& result, size_t& N )
{
  if( field.has_levels() )
  {
    const size_t nvar = field.stride(1);
    const array::LocalView<T,3> arr = make_leveled_view<T>(field);

    field::Field surface_field = fs.createField<T>("surface",array::make_shape(nvar));
    array::LocalView<T,2> surface = surface_view<T>( surface_field );

    for( size_t n=0; n<arr.shape(0); ++n ) {
      for( size_t l=0; l<arr.shape(1); ++l ) {
        for( size_t j=0; j<arr.shape(2); ++j ) {
          surface(n,j) += arr(n,l,j);
        }
      }
    }

    dispatch_order_independent_sum_2d( fs, surface_field, result, N );
    N *= arr.shape(1);
  }
  else
  {
    dispatch_order_independent_sum_2d( fs, field, result, N );
  }
}

template< typename T >
void order_independent_sum( const NodeColumns& fs, const field::Field& field, std::vector<T>& result, size_t& N )
{
  if( field.datatype() == array::DataType::kind<T>() ) {
    return dispatch_order_independent_sum(fs,field,result,N);
  }
  else
  {
    switch( field.datatype().kind() )
    {
      case array::DataType::KIND_INT32 : {
          std::vector<int> tmp;
          dispatch_order_independent_sum(fs,field,tmp,N);
          result.assign(tmp.begin(),tmp.end());
          return;
      }
      case array::DataType::KIND_INT64 : {
        std::vector<long> tmp;
        dispatch_order_independent_sum(fs,field,tmp,N);
        result.assign(tmp.begin(),tmp.end());
        return;
      }
      case array::DataType::KIND_REAL32 : {
        std::vector<float> tmp;
        dispatch_order_independent_sum(fs,field,tmp,N);
        result.assign(tmp.begin(),tmp.end());
        return;
      }
      case array::DataType::KIND_REAL64 : {
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
void dispatch_order_independent_sum_per_level( const NodeColumns& fs, const field::Field& field, field::Field& sumfield, size_t& N )
{
  array::ArrayShape shape;
  shape.reserve(field.rank()-1);
  for( size_t j=1; j<field.rank(); ++j )
    shape.push_back(field.shape(j));
  sumfield.resize(shape);

  array::LocalView<T,2> sum ( array::make_storageview<T>(sumfield).data(),
                              array::make_shape(sumfield.shape(0),sumfield.stride(0)) );
  for( size_t l=0; l<sum.shape(0); ++l ) {
    for( size_t j=0; j<sum.shape(1); ++j ) {
      sum(l,j) = 0.;
    }
  }
  Log::info() << field << std::endl;
  Log::info() << sumfield << std::endl;

  size_t root = 0;
  field::Field global = fs.createField("global",field,field::global());

  Log::info() << global << std::endl;

  fs.gather(field,global);
  if( parallel::mpi::comm().rank() == 0 ) {
    const array::LocalView<T,3> glb = make_leveled_view<T>(global);

    for( size_t n=0; n<glb.shape(0); ++n ) {
      for( size_t l=0; l<glb.shape(1); ++l ) {
        for( size_t j=0; j<glb.shape(2); ++j ) {
          sum(l,j) += glb(n,l,j);
        }
      }
    }
  }
  parallel::mpi::comm().broadcast( array::make_storageview<T>(sumfield).data(),sumfield.size(),root);
  N = fs.nb_nodes_global();
}

void order_independent_sum_per_level( const NodeColumns& fs, const field::Field& field, field::Field& sum, size_t& N )
{
  if( field.datatype() != sum.datatype() ) {
    throw eckit::Exception("field::Field and sum are not of same datatype.",Here());
  }
  switch( field.datatype().kind() )
  {
    case array::DataType::KIND_INT32 :
      return dispatch_order_independent_sum_per_level<int>(fs,field,sum,N);
    case array::DataType::KIND_INT64 :
      return dispatch_order_independent_sum_per_level<long>(fs,field,sum,N);
    case array::DataType::KIND_REAL32 :
      return dispatch_order_independent_sum_per_level<float>(fs,field,sum,N);
    case array::DataType::KIND_REAL64 :
      return dispatch_order_independent_sum_per_level<double>(fs,field,sum,N);
    default: throw eckit::Exception("datatype not supported",Here());
  }
}

template< typename T >
void dispatch_minimum( const NodeColumns& fs, const field::Field& field, std::vector<T>& min )
{
  const array::LocalView<T,3> arr = make_leveled_view<T>(field);
  const size_t nvar = arr.shape(2);
  min.resize(nvar);
  std::vector<T> local_minimum(nvar,std::numeric_limits<T>::max());
  atlas_omp_parallel
  {
    std::vector<T> local_minimum_private(nvar,std::numeric_limits<T>::max());
    const size_t npts = arr.shape(0);
    atlas_omp_for( size_t n=0; n<npts; ++n ) {
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

  parallel::mpi::comm().allReduce(local_minimum, min, eckit::mpi::min());
}

template< typename T >
void minimum( const NodeColumns& fs, const field::Field& field, std::vector<T>& min )
{
  if( field.datatype() == array::DataType::kind<T>() ) {
    return dispatch_minimum(fs,field,min);
  }
  else
  {
    switch( field.datatype().kind() )
    {
      case array::DataType::KIND_INT32 : {
        std::vector<int> tmp;
        dispatch_minimum(fs,field,tmp);
        min.assign(tmp.begin(),tmp.end());
        return;
      }
      case array::DataType::KIND_INT64 : {
        std::vector<long> tmp;
        dispatch_minimum(fs,field,tmp);
        min.assign(tmp.begin(),tmp.end());
        return;
      }
      case array::DataType::KIND_REAL32 : {
        std::vector<float> tmp;
        dispatch_minimum(fs,field,tmp);
        min.assign(tmp.begin(),tmp.end());
        return;
      }
      case array::DataType::KIND_REAL64 : {
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
void dispatch_maximum( const NodeColumns& fs, const field::Field& field, std::vector<T>& max )
{
  const array::LocalView<T,3> arr = make_leveled_view<T>(field);
  const size_t nvar = arr.shape(2);
  max.resize(nvar);
  std::vector<T> local_maximum(nvar,-std::numeric_limits<T>::max());
  atlas_omp_parallel
  {
    std::vector<T> local_maximum_private(nvar,-std::numeric_limits<T>::max());
    const size_t npts = arr.shape(0);
    atlas_omp_for( size_t n=0; n<npts; ++n ) {
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
  parallel::mpi::comm().allReduce(local_maximum, max, eckit::mpi::max());
}

template< typename T >
void maximum( const NodeColumns& fs, const field::Field& field, std::vector<T>& max )
{
  if( field.datatype() == array::DataType::kind<T>() ) {
    return dispatch_maximum(fs,field,max);
  }
  else
  {
    switch( field.datatype().kind() )
    {
      case array::DataType::KIND_INT32 : {
        std::vector<int> tmp;
        dispatch_maximum(fs,field,tmp);
        max.assign(tmp.begin(),tmp.end());
        return;
      }
      case array::DataType::KIND_INT64 : {
        std::vector<long> tmp;
        dispatch_maximum(fs,field,tmp);
        max.assign(tmp.begin(),tmp.end());
        return;
      }
      case array::DataType::KIND_REAL32 : {
        std::vector<float> tmp;
        dispatch_maximum(fs,field,tmp);
        max.assign(tmp.begin(),tmp.end());
        return;
      }
      case array::DataType::KIND_REAL64 : {
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
void minimum( const NodeColumns& fs, const field::Field& field, T& min )
{
  std::vector<T> v;
  minimum(fs,field,v);
  min = v[0];
}

template< typename T >
void maximum( const NodeColumns& fs, const field::Field& field, T& max )
{
  std::vector<T> v;
  maximum(fs,field,v);
  max = v[0];
}

template< typename T >
void dispatch_minimum_per_level( const NodeColumns& fs, const field::Field& field, field::Field& min_field )
{
  array::ArrayShape shape;
  shape.reserve(field.rank()-1);
  for( size_t j=1; j<field.rank(); ++j )
    shape.push_back(field.shape(j));
  min_field.resize(shape);
  const size_t nvar = field.stride(1);
  array::LocalView<T,2> min(
      array::make_storageview<T>(min_field).data(),
      array::make_shape(min_field.shape(0),min_field.stride(0)) );

  for( size_t l=0; l<min.shape(0); ++l ) {
    for( size_t j=0; j<min.shape(1); ++j ) {
      min(l,j) = std::numeric_limits<T>::max();
    }
  }

  const array::LocalView<T,3> arr = make_leveled_view<T>(field);
  atlas_omp_parallel
  {
    array::ArrayT<T> min_private( min.shape(0),min.shape(1) );
    array::ArrayView<T,2> min_private_view = array::make_view<T,2>(min_private);
    for( size_t l=0; l<min.shape(0); ++l ) {
      for( size_t j=0; j<min.shape(1); ++j ) {
        min_private_view(l,j) = std::numeric_limits<T>::max();
      }
    }

    const size_t npts = arr.shape(0);
    atlas_omp_for( size_t n=0; n<npts; ++n ) {
      for( size_t l=0; l<arr.shape(1); ++l ) {
        for( size_t j=0; j<nvar; ++j ) {
          min_private_view(l,j) = std::min(arr(n,l,j),min_private_view(l,j));
        }
      }
    }
    atlas_omp_critical
    {
      for( size_t l=0; l<arr.shape(1); ++l ) {
        for( size_t j=0; j<nvar; ++j ) {
          min(l,j) = std::min(min_private_view(l,j),min(l,j));
        }
      }
    }
  }
  parallel::mpi::comm().allReduceInPlace(min.data(),min_field.size(),eckit::mpi::min());
}

void minimum_per_level( const NodeColumns& fs, const field::Field& field, field::Field& min )
{
  if( field.datatype() != min.datatype() ) {
    throw eckit::Exception("field::Field and min are not of same datatype.",Here());
  }
  switch( field.datatype().kind() )
  {
    case array::DataType::KIND_INT32 :
      return dispatch_minimum_per_level<int>(fs,field,min);
    case array::DataType::KIND_INT64 :
      return dispatch_minimum_per_level<long>(fs,field,min);
    case array::DataType::KIND_REAL32 :
      return dispatch_minimum_per_level<float>(fs,field,min);
    case array::DataType::KIND_REAL64 :
      return dispatch_minimum_per_level<double>(fs,field,min);
    default: throw eckit::Exception("datatype not supported",Here());
  }
}

template< typename T >
void dispatch_maximum_per_level( const NodeColumns& fs, const field::Field& field, field::Field& max_field )
{
  array::ArrayShape shape;
  shape.reserve(field.rank()-1);
  for( size_t j=1; j<field.rank(); ++j )
    shape.push_back(field.shape(j));
  max_field.resize(shape);
  const size_t nvar = field.stride(1);
  array::LocalView<T,2> max (
      array::make_storageview<T>(max_field).data(),
      array::make_shape(max_field.shape(0),max_field.stride(0)) );

  for( size_t l=0; l<max.shape(0); ++l ) {
    for( size_t j=0; j<max.shape(1); ++j ) {
      max(l,j) = -std::numeric_limits<T>::max();
    }
  }

  const array::LocalView<T,3> arr = make_leveled_view<T>(field);
  atlas_omp_parallel
  {
    array::ArrayT<T> max_private(max.shape(0),max.shape(1));
    array::ArrayView<T,2> max_private_view = array::make_view<T,2>(max_private);

    for( size_t l=0; l<max_private_view.shape(0); ++l ) {
      for( size_t j=0; j<max_private_view.shape(1); ++j ) {
        max_private_view(l,j) = -std::numeric_limits<T>::max();
      }
    }

    const size_t npts = arr.shape(0);
    atlas_omp_for( size_t n=0; n<npts; ++n ) {
      for( size_t l=0; l<arr.shape(1); ++l ) {
        for( size_t j=0; j<nvar; ++j ) {
          max_private_view(l,j) = std::max(arr(n,l,j),max_private_view(l,j));
        }
      }
    }
    atlas_omp_critical
    {
      for( size_t l=0; l<arr.shape(1); ++l ) {
        for( size_t j=0; j<nvar; ++j ) {
          max(l,j) = std::max(max_private_view(l,j),max(l,j));
        }
      }
    }
  }
  parallel::mpi::comm().allReduceInPlace(max.data(),max_field.size(),eckit::mpi::max());
}

void maximum_per_level( const NodeColumns& fs, const field::Field& field, field::Field& max )
{
  if( field.datatype() != max.datatype() ) {
    throw eckit::Exception("field::Field and max are not of same datatype.",Here());
  }
  switch( field.datatype().kind() )
  {
    case array::DataType::KIND_INT32 :
      return dispatch_maximum_per_level<int>(fs,field,max);
    case array::DataType::KIND_INT64 :
      return dispatch_maximum_per_level<long>(fs,field,max);
    case array::DataType::KIND_REAL32 :
      return dispatch_maximum_per_level<float>(fs,field,max);
    case array::DataType::KIND_REAL64 :
      return dispatch_maximum_per_level<double>(fs,field,max);
    default: throw eckit::Exception("datatype not supported",Here());
  }
}


template< typename T >
void dispatch_minimum_and_location( const NodeColumns& fs, const field::Field& field, std::vector<T>& min, std::vector<gidx_t>& glb_idx, std::vector<size_t>& level )
{
  array::LocalView<T,3> arr = make_leveled_view<T>(field);
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
    const size_t npts = arr.shape(0);
    atlas_omp_for( size_t n=0; n<npts; ++n ) {
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
  const array::ArrayView<gidx_t,1> global_index = array::make_view<gidx_t,1>(fs.nodes().global_index());
  for( size_t j=0; j<nvar; ++j ) {
    gidx_t glb_idx = global_index(loc_node[j]);
    ASSERT( glb_idx < std::numeric_limits<int>::max() ); // pairs with 64bit integer for second not implemented
    min_and_gidx_loc[j] = std::make_pair(local_minimum[j],glb_idx);
    min_and_level_loc[j] = std::make_pair(local_minimum[j],loc_level[j]);
  }

  parallel::mpi::comm().allReduce(min_and_gidx_loc, min_and_gidx_glb, eckit::mpi::minloc());
  parallel::mpi::comm().allReduce(min_and_level_loc,min_and_level_glb,eckit::mpi::minloc());

  for( size_t j=0; j<nvar; ++j ) {
    min[j]     = min_and_gidx_glb[j].first;
    glb_idx[j] = min_and_gidx_glb[j].second;
    level[j]   = min_and_level_glb[j].second;
  }

}

template< typename T >
void minimum_and_location( const NodeColumns& fs, const field::Field& field, std::vector<T>& min, std::vector<gidx_t>& glb_idx, std::vector<size_t>& level )
{
  if( field.datatype() == array::DataType::kind<T>() ) {
    return dispatch_minimum_and_location(fs,field,min,glb_idx,level);
  }
  else
  {
    switch( field.datatype().kind() )
    {
      case array::DataType::KIND_INT32 : {
        std::vector<int> tmp;
        dispatch_minimum_and_location(fs,field,tmp,glb_idx,level);
        min.assign(tmp.begin(),tmp.end());
        return;
      }
      case array::DataType::KIND_INT64 : {
        std::vector<long> tmp;
        dispatch_minimum_and_location(fs,field,tmp,glb_idx,level);
        min.assign(tmp.begin(),tmp.end());
        return;
      }
      case array::DataType::KIND_REAL32 : {
        std::vector<float> tmp;
        dispatch_minimum_and_location(fs,field,tmp,glb_idx,level);
        min.assign(tmp.begin(),tmp.end());
        return;
      }
      case array::DataType::KIND_REAL64 : {
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
void dispatch_maximum_and_location( const NodeColumns& fs, const field::Field& field, std::vector<T>& max, std::vector<gidx_t>& glb_idx, std::vector<size_t>& level )
{
  array::LocalView<T,3> arr = make_leveled_view<T>(field);
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
    const size_t npts = arr.shape(0);
    atlas_omp_for( size_t n=0; n<npts; ++n ) {
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
  const array::ArrayView<gidx_t,1> global_index = array::make_view<gidx_t,1>( fs.nodes().global_index() );
  for( size_t j=0; j<nvar; ++j ) {
    gidx_t glb_idx = global_index(loc_node[j]);
    ASSERT( glb_idx < std::numeric_limits<int>::max() ); // pairs with 64bit integer for second not implemented
    max_and_gidx_loc[j] = std::make_pair(local_maximum[j],glb_idx);
    max_and_level_loc[j] = std::make_pair(local_maximum[j],loc_level[j]);
  }

  parallel::mpi::comm().allReduce(max_and_gidx_loc, max_and_gidx_glb, eckit::mpi::maxloc());
  parallel::mpi::comm().allReduce(max_and_level_loc,max_and_level_glb,eckit::mpi::maxloc());

  for( size_t j=0; j<nvar; ++j ) {
    max[j]     = max_and_gidx_glb[j].first;
    glb_idx[j] = max_and_gidx_glb[j].second;
    level[j]   = max_and_level_glb[j].second;
  }
}

template< typename T >
void maximum_and_location( const NodeColumns& fs, const field::Field& field, std::vector<T>& max, std::vector<gidx_t>& glb_idx, std::vector<size_t>& level )
{
  if( field.datatype() == array::DataType::kind<T>() ) {
    return dispatch_maximum_and_location(fs,field,max,glb_idx,level);
  }
  else
  {
    switch( field.datatype().kind() )
    {
      case array::DataType::KIND_INT32 : {
        std::vector<int> tmp;
        dispatch_maximum_and_location(fs,field,tmp,glb_idx,level);
        max.assign(tmp.begin(),tmp.end());
        return;
      }
      case array::DataType::KIND_INT64 : {
        std::vector<long> tmp;
        dispatch_maximum_and_location(fs,field,tmp,glb_idx,level);
        max.assign(tmp.begin(),tmp.end());
        return;
      }
      case array::DataType::KIND_REAL32 : {
        std::vector<float> tmp;
        dispatch_maximum_and_location(fs,field,tmp,glb_idx,level);
        max.assign(tmp.begin(),tmp.end());
        return;
      }
      case array::DataType::KIND_REAL64 : {
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
void minimum_and_location( const NodeColumns& fs, const field::Field& field, std::vector<T>& min, std::vector<gidx_t>& glb_idx)
{
  std::vector<size_t> level;
  minimum_and_location(fs,field,min,glb_idx,level);
}

template< typename T >
void maximum_and_location( const NodeColumns& fs, const field::Field& field, std::vector<T>& max, std::vector<gidx_t>& glb_idx)
{
  std::vector<size_t> level;
  maximum_and_location(fs,field,max,glb_idx,level);
}

template< typename T >
void minimum_and_location( const NodeColumns& fs, const field::Field& field, T& min, gidx_t& glb_idx, size_t& level)
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
void maximum_and_location( const NodeColumns& fs, const field::Field& field, T& max, gidx_t& glb_idx, size_t& level)
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
void minimum_and_location( const NodeColumns& fs, const field::Field& field, T& min, gidx_t& glb_idx)
{
  size_t level;
  minimum_and_location(fs,field,min,glb_idx,level);
}

template< typename T >
void maximum_and_location( const NodeColumns& fs, const field::Field& field, T& max, gidx_t& glb_idx)
{
  size_t level;
  maximum_and_location(fs,field,max,glb_idx,level);
}

template< typename T >
void dispatch_minimum_and_location_per_level( const NodeColumns& fs, const field::Field& field, field::Field& min_field, field::Field& glb_idx_field )
{
  const array::LocalView<T,3> arr = make_leveled_view<T>(field);
  array::ArrayShape shape;
  shape.reserve(field.rank()-1);
  for( size_t j=1; j<field.rank(); ++j )
    shape.push_back(field.shape(j));
  min_field.resize(shape);
  glb_idx_field.resize(shape);
  const size_t nvar = arr.shape(2);
  array::LocalView<T,2> min(
      array::make_storageview<T>(min_field).data(),
      array::make_shape(min_field.shape(0),min_field.stride(0)) );

  for( size_t l=0; l<min.shape(0); ++l ) {
    for( size_t j=0; j<min.shape(1); ++j ) {
      min(l,j) = std::numeric_limits<T>::max();
    }
  }

  array::LocalView<gidx_t,2> glb_idx(
    array::make_storageview<gidx_t>(glb_idx_field).data(),
    array::make_shape(glb_idx_field.shape(0),glb_idx_field.stride(0)) );


  atlas_omp_parallel
  {
    array::ArrayT<T> min_private(min.shape(0),min.shape(1));
    array::ArrayView<T,2> min_private_view = array::make_view<T,2>(min_private);

    for( size_t l=0; l<min_private_view.shape(0); ++l ) {
      for( size_t j=0; j<min_private_view.shape(1); ++j ) {
        min_private_view(l,j) = std::numeric_limits<T>::max();
      }
    }

    array::ArrayT<T> glb_idx_private(glb_idx.shape(0),glb_idx.shape(1));
    array::ArrayView<gidx_t,2> glb_idx_private_view = array::make_view<gidx_t,2>(glb_idx_private);
    const size_t npts = arr.shape(0);
    atlas_omp_for( size_t n=0; n<npts; ++n ) {
      for( size_t l=0; l<arr.shape(1); ++l ) {
        for( size_t j=0; j<nvar; ++j ) {
          if( arr(n,l,j) < min(l,j) ) {
            min_private_view(l,j) = arr(n,l,j);
            glb_idx_private_view(l,j) = n;
          }
        }
      }
    }
    atlas_omp_critical_ordered
    {
      for( size_t l=0; l<arr.shape(1); ++l ) {
        for( size_t j=0; j<nvar; ++j ) {
          if( min_private_view(l,j) < min(l,j) ) {
            min(l,j) = min_private_view(l,j);
            glb_idx(l,j) = glb_idx_private_view(l,j);
          }
        }
      }
    }
  }
  const size_t nlev = arr.shape(1);
  std::vector< std::pair<T,int> > min_and_gidx_loc(nlev*nvar);
  std::vector< std::pair<T,int> > min_and_gidx_glb(nlev*nvar);
  const array::ArrayView<gidx_t,1> global_index = array::make_view<gidx_t,1>( fs.nodes().global_index() );
  atlas_omp_parallel_for( size_t l=0; l<nlev; ++l ) {
    for( size_t j=0; j<nvar; ++j ) {
      gidx_t gidx = global_index(glb_idx(l,j));
      ASSERT( gidx < std::numeric_limits<int>::max() ); // pairs with 64bit integer for second not implemented
      min_and_gidx_loc[j+nvar*l] = std::make_pair(min(l,j),gidx);
    }
  }

  parallel::mpi::comm().allReduce(min_and_gidx_loc,min_and_gidx_glb,eckit::mpi::minloc());

  atlas_omp_parallel_for( size_t l=0; l<nlev; ++l ) {
    for( size_t j=0; j<nvar; ++j ) {
      min(l,j)     = min_and_gidx_glb[j+l*nvar].first;
      glb_idx(l,j) = min_and_gidx_glb[j+l*nvar].second;
    }
  }
}

void minimum_and_location_per_level( const NodeColumns& fs, const field::Field& field, field::Field& min, field::Field& glb_idx )
{
  if( field.datatype() != min.datatype() ) {
    throw eckit::Exception("field::Field and min are not of same datatype.",Here());
  }
  if( glb_idx.datatype() != array::DataType::kind<gidx_t>() ) {
    throw eckit::Exception("glb_idx field::Field is not of correct datatype",Here());
  }
  switch( field.datatype().kind() )
  {
    case array::DataType::KIND_INT32 :
      return dispatch_minimum_and_location_per_level<int>(fs,field,min,glb_idx);
    case array::DataType::KIND_INT64 :
      return dispatch_minimum_and_location_per_level<long>(fs,field,min,glb_idx);
    case array::DataType::KIND_REAL32 :
      return dispatch_minimum_and_location_per_level<float>(fs,field,min,glb_idx);
    case array::DataType::KIND_REAL64 :
      return dispatch_minimum_and_location_per_level<double>(fs,field,min,glb_idx);
    default: throw eckit::Exception("datatype not supported",Here());
  }
}


template< typename T >
void dispatch_maximum_and_location_per_level( const NodeColumns& fs, const field::Field& field, field::Field& max_field, field::Field& glb_idx_field )
{
  const array::LocalView<T,3> arr = make_leveled_view<T>(field);
  array::ArrayShape shape;
  shape.reserve(field.rank()-1);
  for( size_t j=1; j<field.rank(); ++j )
    shape.push_back(field.shape(j));
  max_field.resize(shape);
  glb_idx_field.resize(shape);
  const size_t nvar = arr.shape(2);
  array::LocalView<T,2> max(
      array::make_storageview<T>(max_field).data(),
      array::make_shape(max_field.shape(0),max_field.stride(0)) );

  for( size_t l=0; l<max.shape(0); ++l ) {
    for( size_t j=0; j<max.shape(1); ++j ) {
      max(l,j) = -std::numeric_limits<T>::max();
    }
  }

  array::LocalView<gidx_t,2> glb_idx(
      array::make_storageview<gidx_t>(glb_idx_field).data(),
      array::make_shape(glb_idx_field.shape(0),glb_idx_field.stride(0)) );

  atlas_omp_parallel
  {
    array::ArrayT<T> max_private(max.shape(0),max.shape(1));
    array::ArrayView<T,2> max_private_view = array::make_view<T,2>(max_private);

    for( size_t l=0; l<max_private_view.shape(0); ++l ) {
      for( size_t j=0; j<max_private_view.shape(1); ++j ) {
        max_private_view(l,j) = -std::numeric_limits<T>::max();
      }
    }

    array::ArrayT<T> glb_idx_private(glb_idx.shape(0),glb_idx.shape(1));
    array::ArrayView<gidx_t,2> glb_idx_private_view = array::make_view<gidx_t,2>(glb_idx_private);
    const size_t npts = arr.shape(0);
    atlas_omp_for( size_t n=0; n<npts; ++n ) {
      for( size_t l=0; l<arr.shape(1); ++l ) {
        for( size_t j=0; j<nvar; ++j ) {
          if( arr(n,l,j) > max(l,j) ) {
            max_private_view(l,j) = arr(n,l,j);
            glb_idx_private_view(l,j) = n;
          }
        }
      }
    }
    atlas_omp_critical_ordered
    {
      for( size_t l=0; l<arr.shape(1); ++l ) {
        for( size_t j=0; j<nvar; ++j ) {
          if( max_private_view(l,j) > max(l,j) ) {
            max(l,j) = max_private_view(l,j);
            glb_idx(l,j) = glb_idx_private_view(l,j);
          }
        }
      }
    }
  }

  const size_t nlev = arr.shape(1);
  std::vector< std::pair<T,int> > max_and_gidx_loc(nlev*nvar);
  std::vector< std::pair<T,int> > max_and_gidx_glb(nlev*nvar);
  const array::ArrayView<gidx_t,1> global_index = array::make_view<gidx_t,1>( fs.nodes().global_index() );
  atlas_omp_parallel_for( size_t l=0; l<nlev; ++l ) {
    for( size_t j=0; j<nvar; ++j ) {
      gidx_t gidx = global_index(glb_idx(l,j));
      ASSERT( gidx < std::numeric_limits<int>::max() ); // pairs with 64bit integer for second not implemented
      max_and_gidx_loc[j+nvar*l] = std::make_pair(max(l,j),gidx);
    }
  }

  parallel::mpi::comm().allReduce(max_and_gidx_loc,max_and_gidx_glb,eckit::mpi::maxloc());

  atlas_omp_parallel_for( size_t l=0; l<nlev; ++l ) {
    for( size_t j=0; j<nvar; ++j ) {
      max(l,j)     = max_and_gidx_glb[j+l*nvar].first;
      glb_idx(l,j) = max_and_gidx_glb[j+l*nvar].second;
    }
  }
}

void maximum_and_location_per_level( const NodeColumns& fs, const field::Field& field, field::Field& max, field::Field& glb_idx )
{
  if( field.datatype() != max.datatype() ) {
    throw eckit::Exception("field::Field and max are not of same datatype.",Here());
  }
  if( glb_idx.datatype() != array::DataType::kind<gidx_t>() ) {
    throw eckit::Exception("glb_idx field::Field is not of correct datatype",Here());
  }
  switch( field.datatype().kind() )
  {
    case array::DataType::KIND_INT32 :
      return dispatch_maximum_and_location_per_level<int>(fs,field,max,glb_idx);
    case array::DataType::KIND_INT64 :
      return dispatch_maximum_and_location_per_level<long>(fs,field,max,glb_idx);
    case array::DataType::KIND_REAL32 :
      return dispatch_maximum_and_location_per_level<float>(fs,field,max,glb_idx);
    case array::DataType::KIND_REAL64 :
      return dispatch_maximum_and_location_per_level<double>(fs,field,max,glb_idx);
    default: throw eckit::Exception("datatype not supported",Here());
  }
}


template< typename T >
void mean( const NodeColumns& fs, const field::Field& field, T& result, size_t& N )
{
  sum(fs,field,result,N);
  result /= static_cast<double>(N);
}

template< typename T >
void mean( const NodeColumns& fs, const field::Field& field, std::vector<T>& result, size_t& N )
{
  sum(fs,field,result,N);
  for( size_t j=0; j<result.size(); ++j ) {
    result[j] /= static_cast<double>(N);
  }
}

template< typename T >
void dispatch_mean_per_level( const NodeColumns& fs, const field::Field& field, field::Field& mean, size_t& N )
{
  dispatch_sum_per_level<T>(fs,field,mean,N);
  T* rawdata = array::make_storageview<T>(mean).data();
  for( size_t j=0; j<mean.size(); ++j ) {
    rawdata[j] /= static_cast<double>(N);
  }
}


void mean_per_level( const NodeColumns& fs, const field::Field& field, field::Field& mean, size_t& N )
{
  if( field.datatype() != mean.datatype() ) {
    throw eckit::Exception("field::Field and sum are not of same datatype.",Here());
  }
  switch( field.datatype().kind() )
  {
    case array::DataType::KIND_INT32 :
      return dispatch_mean_per_level<int>(fs,field,mean,N);
    case array::DataType::KIND_INT64 :
      return dispatch_mean_per_level<long>(fs,field,mean,N);
    case array::DataType::KIND_REAL32 :
      return dispatch_mean_per_level<float>(fs,field,mean,N);
    case array::DataType::KIND_REAL64 :
      return dispatch_mean_per_level<double>(fs,field,mean,N);
    default: throw eckit::Exception("datatype not supported",Here());
  }
}

template< typename T >
void mean_and_standard_deviation( const NodeColumns& fs, const field::Field& field, T& mu, T& sigma, size_t& N )
{
  mean(fs,field,mu,N);
  field::Field squared_diff_field = fs.createField("sqr_diff",field);
  array::LocalView<T,2> squared_diff = make_leveled_scalar_view<T>( squared_diff_field );
  array::LocalView<T,2> values = make_leveled_scalar_view<T>( field );

  const size_t npts = std::min(values.shape(0),fs.nb_nodes());
  atlas_omp_parallel_for( size_t n=0; n<npts; ++n ) {
    for( size_t l=0; l<values.shape(1); ++l ) {
      squared_diff(n,l) = sqr( values(n,l) - mu );
    }
  }
  mean(fs,squared_diff_field,sigma,N);
  sigma = std::sqrt(sigma);
}

template< typename T >
void mean_and_standard_deviation( const NodeColumns& fs, const field::Field& field, std::vector<T>& mu, std::vector<T>& sigma, size_t& N )
{
  mean(fs,field,mu,N);
  field::Field squared_diff_field = fs.createField("sqr_diff",field);
  array::LocalView<T,3> squared_diff = make_leveled_view<T>( squared_diff_field );
  array::LocalView<T,3> values = make_leveled_view<T>( field );

  const size_t npts = values.shape(0);
  atlas_omp_parallel_for( size_t n=0; n<npts; ++n ) {
    for( size_t l=0; l<values.shape(1); ++l ) {
      for( size_t j=0; j<values.shape(2); ++j ) {
        squared_diff(n,l,j) = sqr( values(n,l,j) - mu[j] );
      }
    }
  }
  mean(fs,squared_diff_field,sigma,N);
  for( size_t j=0; j<sigma.size(); ++j ) {
    sigma[j] = std::sqrt(sigma[j]);
  }
}

template< typename T >
void dispatch_mean_and_standard_deviation_per_level( const NodeColumns& fs, const field::Field& field, field::Field& mean, field::Field& stddev, size_t& N )
{
  dispatch_mean_per_level<T>(fs,field,mean,N);
  field::Field squared_diff_field = fs.createField("sqr_diff",field);
  array::LocalView<T,3> squared_diff = make_leveled_view<T>( squared_diff_field );
  array::LocalView<T,3> values = make_leveled_view<T>( field );
  array::LocalView<T,2> mu( array::make_storageview<T>(mean).data(), array::make_shape(values.shape(1),values.shape(2)) );

  const size_t npts = values.shape(0);
  atlas_omp_parallel_for( size_t n=0; n<npts; ++n ) {
    for( size_t l=0; l<values.shape(1); ++l ) {
      for( size_t j=0; j<values.shape(2); ++j ) {
        squared_diff(n,l,j) = sqr( values(n,l,j) - mu(l,j) );
      }
    }
  }
  dispatch_mean_per_level<T>(fs,squared_diff_field,stddev,N);
  T* sigma = array::make_storageview<T>(stddev).data();
  const size_t size = stddev.size();
  atlas_omp_parallel_for( size_t j=0; j<size; ++j ) {
    sigma[j] = std::sqrt(sigma[j]);
  }
}


void mean_and_standard_deviation_per_level( const NodeColumns& fs, const field::Field& field, field::Field& mean, field::Field& stddev, size_t& N )
{
  if( field.datatype() != mean.datatype() ) {
    throw eckit::Exception("field::Field and mean are not of same datatype.",Here());
  }
  if( field.datatype() != stddev.datatype() ) {
    throw eckit::Exception("field::Field and stddev are not of same datatype.",Here());
  }
  switch( field.datatype().kind() )
  {
    case array::DataType::KIND_INT32 :
      return dispatch_mean_and_standard_deviation_per_level<int>(fs,field,mean,stddev,N);
    case array::DataType::KIND_INT64 :
      return dispatch_mean_and_standard_deviation_per_level<long>(fs,field,mean,stddev,N);
    case array::DataType::KIND_REAL32 :
      return dispatch_mean_and_standard_deviation_per_level<float>(fs,field,mean,stddev,N);
    case array::DataType::KIND_REAL64 :
      return dispatch_mean_and_standard_deviation_per_level<double>(fs,field,mean,stddev,N);
    default: throw eckit::Exception("datatype not supported",Here());
  }
}


} // end collectives implementation




template< typename Value >
NodeColumns::FieldStatisticsT<Value>::FieldStatisticsT(const NodeColumns* f) :
  functionspace(*f) {
}

template< typename Vector >
NodeColumns::FieldStatisticsVectorT<Vector>::FieldStatisticsVectorT(const NodeColumns* f) :
  functionspace(*f) {
}

NodeColumns::FieldStatistics::FieldStatistics(const NodeColumns* f) :
  functionspace(*f) {
}

template< typename Value >
void NodeColumns::FieldStatisticsT<Value>::sum( const field::Field& field, Value& result, size_t& N ) const {
  detail::sum(functionspace,field,result,N);
}

template< typename Vector >
void NodeColumns::FieldStatisticsVectorT<Vector>::sum( const field::Field& field, Vector& result, size_t& N ) const {
  detail::sum(functionspace,field,result,N);
}

void NodeColumns::FieldStatistics::sumPerLevel( const field::Field& field, field::Field& result, size_t& N ) const {
  detail::sum_per_level(functionspace,field,result,N);
}

template< typename Value >
void NodeColumns::FieldStatisticsT<Value>::orderIndependentSum( const field::Field& field, Value& result, size_t& N ) const {
  detail::order_independent_sum(functionspace,field,result,N);
}

template< typename Vector >
void NodeColumns::FieldStatisticsVectorT<Vector>::orderIndependentSum( const field::Field& field, Vector& result, size_t& N ) const {
  detail::order_independent_sum(functionspace,field,result,N);
}

void NodeColumns::FieldStatistics::orderIndependentSumPerLevel( const field::Field& field, field::Field& result, size_t& N ) const {
  detail::order_independent_sum_per_level(functionspace,field,result,N);
}

template< typename Value >
void NodeColumns::FieldStatisticsT<Value>::minimum( const field::Field& field, Value& minimum ) const {
  detail::minimum(functionspace,field,minimum);
}

template< typename Value >
void NodeColumns::FieldStatisticsT<Value>::maximum( const field::Field& field, Value& maximum ) const {
  detail::maximum(functionspace,field,maximum);

}

template< typename Vector >
void NodeColumns::FieldStatisticsVectorT<Vector>::minimum( const field::Field& field, Vector& minimum ) const {
  detail::minimum(functionspace,field,minimum);
}

template< typename Vector >
void NodeColumns::FieldStatisticsVectorT<Vector>::maximum( const field::Field& field, Vector& maximum) const {
  detail::maximum(functionspace,field,maximum);
}

void NodeColumns::FieldStatistics::minimumPerLevel( const field::Field& field, field::Field& minimum ) const {
  detail::minimum_per_level(functionspace,field,minimum);
}

void NodeColumns::FieldStatistics::maximumPerLevel( const field::Field& field, field::Field& maximum ) const {
  detail::maximum_per_level(functionspace,field,maximum);
}

template< typename Value >
void NodeColumns::FieldStatisticsT<Value>::minimumAndLocation( const field::Field& field, Value& minimum, gidx_t& glb_idx ) const {
  detail::minimum_and_location(functionspace,field,minimum,glb_idx);

}

template< typename Value >
void NodeColumns::FieldStatisticsT<Value>::maximumAndLocation( const field::Field& field, Value& maximum, gidx_t& glb_idx ) const {
  detail::maximum_and_location(functionspace,field,maximum,glb_idx);
}

template< typename Value >
void NodeColumns::FieldStatisticsT<Value>::minimumAndLocation( const field::Field& field, Value& minimum, gidx_t& glb_idx, size_t& level ) const {
  detail::minimum_and_location(functionspace,field,minimum,glb_idx,level);
}

template< typename Value >
void NodeColumns::FieldStatisticsT<Value>::maximumAndLocation( const field::Field& field, Value& maximum, gidx_t& glb_idx, size_t& level ) const {
  detail::maximum_and_location(functionspace,field,maximum,glb_idx,level);
}

template< typename Vector >
void NodeColumns::FieldStatisticsVectorT<Vector>::minimumAndLocation( const field::Field& field, Vector& minimum, std::vector<gidx_t>& glb_idx ) const {
  detail::minimum_and_location(functionspace,field,minimum,glb_idx);
}

template< typename Vector >
void NodeColumns::FieldStatisticsVectorT<Vector>::maximumAndLocation( const field::Field& field, Vector& maximum, std::vector<gidx_t>& glb_idx ) const {
  detail::maximum_and_location(functionspace,field,maximum,glb_idx);
}

template< typename Vector >
void NodeColumns::FieldStatisticsVectorT<Vector>::minimumAndLocation( const field::Field& field, Vector& minimum, std::vector<gidx_t>& glb_idx, std::vector<size_t>& level ) const {
  detail::minimum_and_location(functionspace,field,minimum,glb_idx,level);
}

template< typename Vector >
void NodeColumns::FieldStatisticsVectorT<Vector>::maximumAndLocation( const field::Field& field, Vector& maximum, std::vector<gidx_t>& glb_idx, std::vector<size_t>& level ) const {
  detail::maximum_and_location(functionspace,field,maximum,glb_idx,level);
}

void NodeColumns::FieldStatistics::minimumAndLocationPerLevel( const field::Field& field, field::Field& column, field::Field& glb_idx ) const {
  detail::minimum_and_location_per_level(functionspace,field,column,glb_idx);
}

void NodeColumns::FieldStatistics::maximumAndLocationPerLevel( const field::Field& field, field::Field& column, field::Field& glb_idx ) const {
  detail::maximum_and_location_per_level(functionspace,field,column,glb_idx);
}

template< typename Value >
void NodeColumns::FieldStatisticsT<Value>::mean( const field::Field& field, Value& mean, size_t& N ) const {
  detail::mean(functionspace,field,mean,N);
}

template< typename Vector >
void NodeColumns::FieldStatisticsVectorT<Vector>::mean( const field::Field& field, Vector& mean, size_t& N ) const {
  detail::mean(functionspace,field,mean,N);
}

void NodeColumns::FieldStatistics::meanPerLevel(const field::Field &field, field::Field &mean, size_t &N) const {
  detail::mean_per_level(functionspace,field,mean,N);
}

template< typename Value >
void NodeColumns::FieldStatisticsT<Value>::meanAndStandardDeviation( const field::Field& field, Value& mean, Value& stddev, size_t& N ) const {
  detail::mean_and_standard_deviation(functionspace,field,mean,stddev,N);
}

template< typename Vector >
void NodeColumns::FieldStatisticsVectorT<Vector>::meanAndStandardDeviation( const field::Field& field, Vector& mean, Vector& stddev, size_t& N ) const {
  detail::mean_and_standard_deviation(functionspace,field,mean,stddev,N);
}

void NodeColumns::FieldStatistics::meanAndStandardDeviationPerLevel( const field::Field& field, field::Field& mean, field::Field& stddev, size_t& N ) const {
  detail::mean_and_standard_deviation_per_level(functionspace,field,mean,stddev,N);
}

template class NodeColumns::FieldStatisticsT<int>;
template class NodeColumns::FieldStatisticsT<long>;
template class NodeColumns::FieldStatisticsT<float>;
template class NodeColumns::FieldStatisticsT<double>;
//template class NodeColumns::FieldStatisticsT<unsigned long>;
template class NodeColumns::FieldStatisticsVectorT< std::vector<int> >;
template class NodeColumns::FieldStatisticsVectorT< std::vector<long> >;
template class NodeColumns::FieldStatisticsVectorT< std::vector<float> >;
template class NodeColumns::FieldStatisticsVectorT< std::vector<double> >;
//template class NodeColumns::FieldStatisticsVectorT< std::vector<unsigned long> >;

} // namespace detail


NodeColumns::NodeColumns() :
  FunctionSpace(),
  functionspace_(nullptr) {
}

NodeColumns::NodeColumns( const FunctionSpace& functionspace ) :
  FunctionSpace( functionspace ),
  functionspace_( dynamic_cast< const detail::NodeColumns* >( get() ) ) {
}

NodeColumns::NodeColumns( mesh::Mesh& mesh, const mesh::Halo& halo, const eckit::Parametrisation& config ) :
  FunctionSpace( new detail::NodeColumns(mesh,halo,config) ),
  functionspace_( dynamic_cast< const detail::NodeColumns* >( get() ) ) {
}

NodeColumns::NodeColumns( mesh::Mesh& mesh, const mesh::Halo& halo ) :
  FunctionSpace( new detail::NodeColumns(mesh,halo) ),
  functionspace_( dynamic_cast< const detail::NodeColumns* >( get() ) ) {
}

NodeColumns::NodeColumns( mesh::Mesh& mesh ) :
  FunctionSpace( new detail::NodeColumns(mesh) ),
  functionspace_( dynamic_cast< const detail::NodeColumns* >( get() ) ) {
}

size_t NodeColumns::nb_nodes() const {
  return functionspace_->nb_nodes();
}

size_t NodeColumns::nb_nodes_global() const { // All MPI ranks will have same output
  return functionspace_->nb_nodes_global();
}

const mesh::Mesh& NodeColumns::mesh() const {
  return functionspace_->mesh();
}

mesh::Nodes& NodeColumns::nodes() const{
  return functionspace_->nodes();
}

field::Field NodeColumns::createField(
    const std::string& name,
    array::DataType datatype,
    const eckit::Parametrisation& config ) const {
  return functionspace_->createField(name,datatype,config);
}

field::Field NodeColumns::createField(
    const std::string& name,
    array::DataType, size_t levels,
    const eckit::Parametrisation& config) const {
  return functionspace_->createField(name,levels,config);
}

field::Field NodeColumns::createField(
    const std::string& name,
    array::DataType datatype,
    const std::vector<size_t>& variables,
    const eckit::Parametrisation& config) const {
  return functionspace_->createField(name,datatype,variables,config);
}

field::Field NodeColumns::createField(
    const std::string& name,
    array::DataType, size_t levels,
    const std::vector<size_t>& variables,
    const eckit::Parametrisation& config)  const {
  return functionspace_->createField(name,levels,variables,config);
}

field::Field NodeColumns::createField(
    const std::string& name,
    const field::Field& field,
    const eckit::Parametrisation& config ) const {
  return functionspace_->createField(name,field,config);
}

field::Field NodeColumns::createField(const eckit::Parametrisation& config) const {
  return functionspace_->createField(config);
}

// -- Parallelisation aware methods

const mesh::Halo& NodeColumns::halo() const {
  return functionspace_->halo();
}

void NodeColumns::haloExchange( field::FieldSet& fieldset ) const {
  functionspace_->haloExchange(fieldset);
}

void NodeColumns::haloExchange( field::Field& field ) const {
  functionspace_->haloExchange(field);
}

const parallel::HaloExchange& NodeColumns::halo_exchange() const {
  return functionspace_->halo_exchange();
}

void NodeColumns::gather( const field::FieldSet& local, field::FieldSet& global ) const {
  functionspace_->gather(local,global);
}

void NodeColumns::gather( const field::Field& local, field::Field& global ) const {
  functionspace_->gather(local,global);
}

const parallel::GatherScatter& NodeColumns::gather() const {
  return functionspace_->gather();
}

void NodeColumns::scatter( const field::FieldSet& global, field::FieldSet& local ) const {
  functionspace_->scatter(global,local);
}

void NodeColumns::scatter( const field::Field& global, field::Field& local ) const {
  functionspace_->scatter(global,local);
}

const parallel::GatherScatter& NodeColumns::scatter() const {
  return functionspace_->scatter();
}

std::string NodeColumns::checksum( const field::FieldSet& fieldset ) const {
  return functionspace_->checksum(fieldset);
}

std::string NodeColumns::checksum( const field::Field& field ) const {
  return functionspace_->checksum(field);
}

const parallel::Checksum& NodeColumns::checksum() const {
  return functionspace_->checksum();
}






} // namespace functionspace
} // namespace atlas

