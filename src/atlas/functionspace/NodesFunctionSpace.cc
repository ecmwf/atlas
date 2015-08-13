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

NodesFunctionSpace::NodesFunctionSpace(const std::string& name, Mesh& mesh, const Halo& halo)
  : next::FunctionSpace(name),
    mesh_(mesh),
    nodes_(mesh_.nodes()),
    halo_(halo.size()),
    nb_nodes_(0),
    nb_nodes_global_(0),
    nb_nodes_global_broadcasted_(0)
{
  if( ! mesh_.halo_exchange().has(halo_name()) && halo_ > 0)
  {
    // Create new halo-exchange
    mpl::HaloExchange* halo_exchange = new mpl::HaloExchange( halo_name() );

    // Set it up.
    actions::build_nodes_parallel_fields( mesh_.nodes() );

    actions::build_periodic_boundaries(mesh_);

    actions::build_halo(mesh_,halo_);

    actions::renumber_nodes_glb_idx(mesh_.nodes());

    Field& ridx = mesh_.nodes().remote_index();
    Field& part = mesh_.nodes().partition();

    std::stringstream ss;
    ss << "nb_nodes_including_halo["<<halo_<<"]";
    mesh.metadata().get(ss.str(),nb_nodes_);

    halo_exchange->setup(part.data<int>(),ridx.data<int>(),REMOTE_IDX_BASE,nb_nodes_);

    // Store it in the mesh
    mesh_.halo_exchange().add(halo_exchange);
  }
  if( !nb_nodes_ ) {
    std::stringstream ss;
    ss << "nb_nodes_including_halo["<<halo_<<"]";
    if( ! mesh.metadata().get(ss.str(),nb_nodes_) ) {
      nb_nodes_ = mesh_.nodes().metadata().get<size_t>("nb_owned");
    }
  }

  if( !mesh_.gather_scatter().has(gather_scatter_name()) )
  {
    // Create new gather_scatter
    mpl::GatherScatter* gather_scatter = new mpl::GatherScatter( gather_scatter_name() );

    // Set it up.
    if( halo_ == 0 )
      actions::build_nodes_parallel_fields( mesh_.nodes() );

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

    // Set it up.
    if( halo_ == 0 )
      actions::build_nodes_parallel_fields( mesh_.nodes() );

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
  nb_nodes_global_broadcasted_ = mesh_.gather_scatter().get(gather_scatter_name()).glb_dof(0);
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

std::string NodesFunctionSpace::halo_name() const
{
  std::stringstream ss; ss << "nodes_" << halo_;
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

Field* NodesFunctionSpace::createField(const std::string& datatype) const {
  return Field::create(make_shape(nb_nodes()),datatype);
}

Field* NodesFunctionSpace::createField(const std::string& name, const std::string& datatype) const {
  return Field::create(name,make_shape(nb_nodes()),datatype);
}

Field* NodesFunctionSpace::createField(const std::vector<size_t>& variables, const std::string& datatype) const {
  std::vector<size_t> shape(1,nb_nodes());
  for( size_t i=0; i<variables.size(); ++i ) shape.push_back(variables[i]);
  return Field::create(shape,datatype);
}

Field* NodesFunctionSpace::createField(const std::string& name, const std::vector<size_t>& variables, const std::string& datatype) const {
  std::vector<size_t> shape(1,nb_nodes());
  for( size_t i=0; i<variables.size(); ++i ) shape.push_back(variables[i]);
  return Field::create(name, shape,datatype);
}

Field* NodesFunctionSpace::createField(const Field& other) const {
  ArrayShape shape = other.shape();
  shape[0] = nb_nodes();
  return Field::create(shape,other.datatype());
}

Field* NodesFunctionSpace::createField(const std::string& name, const Field& other) const {
  ArrayShape shape = other.shape();
  shape[0] = nb_nodes();
  return Field::create(name,shape,other.datatype());
}

Field* NodesFunctionSpace::createGlobalField(const std::string& datatype) const {
  return Field::create(make_shape(nb_nodes_global()),datatype);
}

Field* NodesFunctionSpace::createGlobalField(const std::string& name, const std::string& datatype) const {
  return Field::create(name,make_shape(nb_nodes_global()),datatype);
}

Field* NodesFunctionSpace::createGlobalField(const std::vector<size_t>& variables, const std::string& datatype) const {
  std::vector<size_t> shape(1,nb_nodes_global());
  for( size_t i=0; i<variables.size(); ++i ) shape.push_back(variables[i]);
  return Field::create(shape,datatype);
}

Field* NodesFunctionSpace::createGlobalField(const std::string& name, const std::vector<size_t>& variables, const std::string& datatype) const {
  std::vector<size_t> shape(1,nb_nodes_global());
  for( size_t i=0; i<variables.size(); ++i ) shape.push_back(variables[i]);
  return Field::create(name, shape,datatype);
}

Field* NodesFunctionSpace::createGlobalField(const Field& other) const {
  ArrayShape shape = other.shape();
  shape[0] = nb_nodes_global();
  return Field::create(shape,other.datatype());
}

Field* NodesFunctionSpace::createGlobalField(const std::string& name, const Field& other) const {
  ArrayShape shape = other.shape();
  shape[0] = nb_nodes_global();
  return Field::create(name,shape,other.datatype());
}

void NodesFunctionSpace::haloExchange( FieldSet& fieldset ) const
{
  if( halo_ ) {
    const mpl::HaloExchange& halo_exchange = mesh_.halo_exchange().get(halo_name());
    for( size_t f=0; f<fieldset.size(); ++f ) {
      const Field& field = fieldset[f];
      ArrayStrides strides = make_strides(field.stride(0),1);
      ArrayShape   shape   = make_shape(field.shape(0),field.stride(0));
      if     ( field.kind() == DataType::kind<int>() ) {
        ArrayView<int,2> view(field.data<int>(),strides.data(),shape.data());
        halo_exchange.execute( view );
      }
      else if( field.kind() == DataType::kind<long>() ) {
        ArrayView<long,2> view(field.data<long>(),strides.data(),shape.data());
        halo_exchange.execute( view );
      }
      else if( field.kind() == DataType::kind<float>() ) {
        ArrayView<float,2> view(field.data<float>(),strides.data(),shape.data());
        halo_exchange.execute( view );
      }
      else if( field.kind() == DataType::kind<double>() ) {
        ArrayView<double,2> view(field.data<double>(),strides.data(),shape.data());
        halo_exchange.execute( view );
      }
      else throw eckit::Exception("datatype not supported",Here());
    }
  }
}
void NodesFunctionSpace::haloExchange( Field& field ) const
{
  if( halo_ ) {
    FieldSet fieldset;
    fieldset.add(field);
    haloExchange(fieldset);
  }
}

void NodesFunctionSpace::gather( const FieldSet& local_fieldset, FieldSet& global_fieldset ) const
{
  const mpl::GatherScatter& gather_scatter = mesh_.gather_scatter().get(gather_scatter_name());

  ASSERT(local_fieldset.size() == global_fieldset.size());

  for( size_t f=0; f<local_fieldset.size(); ++f ) {

    const Field& loc = local_fieldset[f];
    Field& glb = global_fieldset[f];

    if     ( loc.kind() == DataType::kind<int>() ) {
      mpl::Field<int const> loc_field(loc.data<int>(),loc.stride(0));
      mpl::Field<int      > glb_field(glb.data<int>(),glb.stride(0));
      gather_scatter.gather( &loc_field, &glb_field, 1 );
    }
    else if( loc.kind() == DataType::kind<long>() ) {
      mpl::Field<long const> loc_field(loc.data<long>(),loc.stride(0));
      mpl::Field<long      > glb_field(glb.data<long>(),glb.stride(0));
      gather_scatter.gather( &loc_field, &glb_field, 1 );
    }
    else if( loc.kind() == DataType::kind<float>() ) {
      mpl::Field<float const> loc_field(loc.data<float>(),loc.stride(0));
      mpl::Field<float      > glb_field(glb.data<float>(),glb.stride(0));
      gather_scatter.gather( &loc_field, &glb_field, 1 );
    }
    else if( loc.kind() == DataType::kind<double>() ) {
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

void NodesFunctionSpace::scatter( const FieldSet& global_fieldset, FieldSet& local_fieldset ) const
{
  const mpl::GatherScatter& gather_scatter = mesh_.gather_scatter().get(gather_scatter_name());

  ASSERT(local_fieldset.size() == global_fieldset.size());

  for( size_t f=0; f<local_fieldset.size(); ++f ) {

    const Field& glb = global_fieldset[f];
    Field& loc = local_fieldset[f];

    if     ( loc.kind() == DataType::kind<int>() ) {
      mpl::Field<int const> glb_field(glb.data<int>(),glb.stride(0));
      mpl::Field<int      > loc_field(loc.data<int>(),loc.stride(0));
      gather_scatter.scatter( &glb_field, &loc_field, 1 );
    }
    else if( loc.kind() == DataType::kind<long>() ) {
      mpl::Field<long const> glb_field(glb.data<long>(),glb.stride(0));
      mpl::Field<long      > loc_field(loc.data<long>(),loc.stride(0));
      gather_scatter.scatter( &glb_field, &loc_field, 1 );
    }
    else if( loc.kind() == DataType::kind<float>() ) {
      mpl::Field<float const> glb_field(glb.data<float>(),glb.stride(0));
      mpl::Field<float      > loc_field(loc.data<float>(),loc.stride(0));
      gather_scatter.scatter( &glb_field, &loc_field, 1 );
    }
    else if( loc.kind() == DataType::kind<double>() ) {
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

std::string NodesFunctionSpace::checksum( const FieldSet& fieldset ) const {
  const mpl::Checksum& checksum = mesh_.checksum().get(checksum_name());

  eckit::MD5 md5;
  for( size_t f=0; f<fieldset.size(); ++f ) {
    const Field& field=fieldset[f];
    if     ( field.kind() == DataType::kind<int>() )
      md5 << checksum.execute( field.data<int>(), field.stride(0) );
    else if( field.kind() == DataType::kind<long>() )
      md5 << checksum.execute( field.data<long>(), field.stride(0) );
    else if( field.kind() == DataType::kind<float>() )
      md5 << checksum.execute( field.data<float>(), field.stride(0) );
    else if( field.kind() == DataType::kind<double>() )
      md5 << checksum.execute( field.data<double>(), field.stride(0) );
    else throw eckit::Exception("datatype not supported",Here());
  }
  return md5;
}
std::string NodesFunctionSpace::checksum( const Field& field ) const {
  FieldSet fieldset;
  fieldset.add(field);
  return checksum(fieldset);
}

namespace { inline double sqr(const double& val) { return val*val; } }

namespace detail { // Collectives implementation

template< typename DATATYPE >
void dispatch_sum( const NodesFunctionSpace& fs, const Field& field, DATATYPE& result, size_t& N )
{
  if( field.size() != fs.nb_nodes() )
    throw eckit::SeriousBug("Cannot sum multi-variable-field into a single scalar value");
  util::IsGhost is_ghost(fs.nodes());
  ArrayView<DATATYPE,1> arr( field.data<DATATYPE>(), make_shape(field.size()) );
  DATATYPE local_sum = 0;
  atlas_omp_pragma( omp parallel for default(shared) reduction(+:local_sum) )
  for( size_t n=0; n<arr.shape(0); ++n ) {
    if( ! is_ghost(n) ) {
      local_sum += arr(n);
    }
  }
  eckit::mpi::all_reduce(local_sum,result,eckit::mpi::sum());
  N = fs.nb_nodes_global_broadcasted_;
}

template< typename DATATYPE >
void sum( const NodesFunctionSpace& fs , const Field& field, DATATYPE& result, size_t& N )
{
  if( field.kind() == DataType::kind<DATATYPE>() ) {
    return dispatch_sum(fs,field,result,N);
  }
  else
  {
    switch( field.kind() )
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


template< typename DATATYPE >
void dispatch_sum( const NodesFunctionSpace& fs, const Field& field, std::vector<DATATYPE>& result, size_t& N )
{
  size_t nvar = field.stride(0);
  ArrayView<DATATYPE,2> arr( field.data<DATATYPE>(), make_shape(field.shape(0),nvar) );
  util::IsGhost is_ghost(fs.nodes());
  std::vector<DATATYPE> local_sum(nvar,0);

  atlas_omp_parallel
  {
    std::vector<DATATYPE> local_sum_private(nvar,0);
    atlas_omp_for( size_t n=0; n<arr.shape(0); ++n )
    {
      if( ! is_ghost(n) ) {
        for( size_t j=0; j<nvar; ++j ) {
          local_sum_private[j] += arr(n,j);
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
  N = fs.nb_nodes_global_broadcasted_;
}

template< typename DATATYPE >
void sum( const NodesFunctionSpace& fs , const Field& field, std::vector<DATATYPE>& result, size_t& N )
{
  if( field.kind() == DataType::kind<DATATYPE>() ) {
    return dispatch_sum(fs,field,result,N);
  }
  else
  {
    switch( field.kind() )
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

template< typename DATATYPE >
void dispatch_order_independent_sum( const NodesFunctionSpace& fs , const Field& field, DATATYPE& result, size_t& N )
{
  size_t root = 0;
  Field::Ptr global( fs.createGlobalField(field) );
  fs.gather(field,*global);
  result = std::accumulate(global->data<DATATYPE>(),global->data<DATATYPE>()+global->size(),0.);
  eckit::mpi::broadcast(result,root);
  N = fs.nb_nodes_global_broadcasted_;
}

template< typename DATATYPE >
void order_independent_sum( const NodesFunctionSpace& fs , const Field& field, DATATYPE& result, size_t& N )
{
  if( field.kind() == DataType::kind<DATATYPE>() ) {
    return dispatch_order_independent_sum(fs,field,result,N);
  }
  else
  {
    switch( field.kind() )
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
void dispatch_order_independent_sum( const NodesFunctionSpace& fs, const Field& field, std::vector<DATATYPE>& result, size_t& N )
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
  N = fs.nb_nodes_global_broadcasted_;
}

template< typename DATATYPE >
void order_independent_sum( const NodesFunctionSpace& fs, const Field& field, std::vector<DATATYPE>& result, size_t& N )
{
  if( field.kind() == DataType::kind<DATATYPE>() ) {
    return dispatch_order_independent_sum(fs,field,result,N);
  }
  else
  {
    switch( field.kind() )
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


template< typename DATATYPE >
void dispatch_minimum( const NodesFunctionSpace& fs, const Field& field, DATATYPE& min )
{
  const DATATYPE* values = field.data<DATATYPE>();
  DATATYPE minimum( std::numeric_limits<DATATYPE>::max() );
  atlas_omp_pragma( omp parallel for default(shared) reduction(min: minimum) )
  for( size_t n=0; n<field.size(); ++n )
    minimum = std::min(values[n],minimum);
  eckit::mpi::all_reduce(minimum,min,eckit::mpi::min());
}

template< typename DATATYPE >
void minimum( const NodesFunctionSpace& fs, const Field& field, DATATYPE& min )
{
  if( field.kind() == DataType::kind<DATATYPE>() ) {
    return dispatch_minimum(fs,field,min);
  }
  else
  {
    switch( field.kind() )
    {
      case DataType::KIND_INT32 : {
        int tmp;
        dispatch_minimum(fs,field,tmp);
        min = tmp;
        return;
      }
      case DataType::KIND_INT64 : {
        long tmp;
        dispatch_minimum(fs,field,tmp);
        min = tmp;
        return;
      }
      case DataType::KIND_REAL32 : {
        float tmp;
        dispatch_minimum(fs,field,tmp);
        min = tmp;
        return;
      }
      case DataType::KIND_REAL64 : {
        double tmp;
        dispatch_minimum(fs,field,tmp);
        min = tmp;
        return;
      }
      default: throw eckit::Exception("datatype not supported",Here());
    }
  }
}

template< typename DATATYPE >
void dispatch_maximum( const NodesFunctionSpace& fs, const Field& field, DATATYPE& max )
{
  const DATATYPE* values = field.data<DATATYPE>();
  DATATYPE maximum( -std::numeric_limits<DATATYPE>::max() );
  atlas_omp_pragma( omp parallel for default(shared) reduction(max: maximum) )
  for( size_t n=0; n<field.size(); ++n )
     maximum = std::max(values[n],maximum);
  eckit::mpi::all_reduce(maximum,max,eckit::mpi::max());
}

template< typename DATATYPE >
void maximum( const NodesFunctionSpace& fs, const Field& field, DATATYPE& max )
{
  if( field.kind() == DataType::kind<DATATYPE>() ) {
    return dispatch_maximum(fs,field,max);
  }
  else
  {
    switch( field.kind() )
    {
      case DataType::KIND_INT32 : {
        int tmp;
        dispatch_maximum(fs,field,tmp);
        max = tmp;
        return;
      }
      case DataType::KIND_INT64 : {
        long tmp;
        dispatch_maximum(fs,field,tmp);
        max = tmp;
        return;
      }
      case DataType::KIND_REAL32 : {
        float tmp;
        dispatch_maximum(fs,field,tmp);
        max = tmp;
        return;
      }
      case DataType::KIND_REAL64 : {
        double tmp;
        dispatch_maximum(fs,field,tmp);
        max = tmp;
        return;
      }
      default: throw eckit::Exception("datatype not supported",Here());
    }
  }
}

template< typename DATATYPE >
void dispatch_minimum( const NodesFunctionSpace& fs, const Field& field, std::vector<DATATYPE>& min )
{
  const size_t nvar = field.stride(0);
  min.resize(nvar);
  std::vector<DATATYPE> local_minimum(nvar,std::numeric_limits<DATATYPE>::max());
  const ArrayView<DATATYPE,2> arr( field.data<DATATYPE>(), make_shape(field.shape(0),nvar) );
  atlas_omp_parallel
  {
    std::vector<DATATYPE> local_minimum_private(nvar,std::numeric_limits<DATATYPE>::max());
    atlas_omp_for( size_t n=0; n<arr.shape(0); ++n ) {
      for( size_t j=0; j<nvar; ++j ) {
        local_minimum[j] = std::min(arr(n,j),local_minimum[j]);
      }
    }
    atlas_omp_critical
    {
      for( size_t j=0; j<nvar; ++j ) {
        local_minimum[j] = std::min(local_minimum_private[j],local_minimum[j]);
      }
    }
  }
  eckit::mpi::all_reduce(local_minimum,min,eckit::mpi::min());
}

template< typename DATATYPE >
void minimum( const NodesFunctionSpace& fs, const Field& field, std::vector<DATATYPE>& min )
{
  if( field.kind() == DataType::kind<DATATYPE>() ) {
    return dispatch_minimum(fs,field,min);
  }
  else
  {
    switch( field.kind() )
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

template< typename DATATYPE >
void dispatch_maximum( const NodesFunctionSpace& fs, const Field& field, std::vector<DATATYPE>& max )
{
  const size_t nvar = field.stride(0);
  max.resize(nvar);
  std::vector<DATATYPE> local_maximum(nvar,-std::numeric_limits<DATATYPE>::max());
  const ArrayView<DATATYPE,2> arr( field.data<DATATYPE>(), make_shape(field.shape(0),nvar) );
  atlas_omp_parallel
  {
    std::vector<DATATYPE> local_maximum_private(nvar,-std::numeric_limits<DATATYPE>::max());
    atlas_omp_for( size_t n=0; n<arr.shape(0); ++n ) {
      for( size_t j=0; j<nvar; ++j ) {
        local_maximum_private[j] = std::max(arr(n,j),local_maximum_private[j]);
      }
    }
    atlas_omp_critical
    {
      for( size_t j=0; j<nvar; ++j ) {
        local_maximum[j] = std::max(local_maximum[j],local_maximum_private[j]);
      }
    }
  }
  eckit::mpi::all_reduce(local_maximum,max,eckit::mpi::max());
}

template< typename DATATYPE >
void maximum( const NodesFunctionSpace& fs, const Field& field, std::vector<DATATYPE>& max )
{
  if( field.kind() == DataType::kind<DATATYPE>() ) {
    return dispatch_maximum(fs,field,max);
  }
  else
  {
    switch( field.kind() )
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

template< typename DATATYPE >
void dispatch_minimum_and_location( const NodesFunctionSpace& fs, const Field& field, DATATYPE& min, gidx_t& glb_idx )
{
  DATATYPE local_minimum ( std::numeric_limits<DATATYPE>::max() );
  size_t location;

  const DATATYPE* values = field.data<DATATYPE>();
  const size_t size = field.size();
  atlas_omp_parallel
  {
    DATATYPE local_minimum_private( std::numeric_limits<DATATYPE>::max() );
    size_t location_private;
    atlas_omp_for( size_t n=0; n<size; ++n )
    {
      if( values[n] < local_minimum_private )
      {
        local_minimum_private = values[n];
        location_private = n;
      }
    }
    atlas_omp_critical_ordered
    {
      if( local_minimum_private < local_minimum )
      {
        local_minimum = local_minimum_private;
        location = location_private;
      }
    }
  }
  glb_idx = fs.nodes().global_index().data<gidx_t>()[location];
  ASSERT( glb_idx < std::numeric_limits<int>::max() ); // pairs with 64bit integer for second not implemented
  const std::pair<DATATYPE,int> min_and_gidx_loc (local_minimum,glb_idx);
  std::pair<DATATYPE,int> min_and_gidx_glb;
  eckit::mpi::all_reduce(min_and_gidx_loc,min_and_gidx_glb,eckit::mpi::minloc());
  min     = min_and_gidx_glb.first;
  glb_idx = min_and_gidx_glb.second;
}

template< typename DATATYPE >
void minimum_and_location( const NodesFunctionSpace& fs, const Field& field, DATATYPE& min, gidx_t& glb_idx )
{
  if( field.kind() == DataType::kind<DATATYPE>() ) {
    return dispatch_minimum_and_location(fs,field,min,glb_idx);
  }
  else
  {
    switch( field.kind() )
    {
      case DataType::KIND_INT32 : {
        int tmp;
        dispatch_minimum_and_location(fs,field,tmp,glb_idx);
        min = tmp;
        return;
      }
      case DataType::KIND_INT64 : {
        long tmp;
        dispatch_minimum_and_location(fs,field,tmp,glb_idx);
        min = tmp;
        return;
      }
      case DataType::KIND_REAL32 : {
        float tmp;
        dispatch_minimum_and_location(fs,field,tmp,glb_idx);
        min = tmp;
        return;
      }
      case DataType::KIND_REAL64 : {
        double tmp;
        dispatch_minimum_and_location(fs,field,tmp,glb_idx);
        min = tmp;
        return;
      }
      default: throw eckit::Exception("datatype not supported",Here());
    }
  }
}

template< typename DATATYPE >
void dispatch_maximum_and_location( const NodesFunctionSpace& fs, const Field& field, DATATYPE& max, gidx_t& glb_idx )
{
  DATATYPE local_maximum ( -std::numeric_limits<DATATYPE>::max() );
  size_t location;

  const DATATYPE* values = field.data<DATATYPE>();
  const size_t size = field.size();
  atlas_omp_parallel
  {
    DATATYPE local_maximum_private( -std::numeric_limits<DATATYPE>::max() );
    size_t location_private;
    atlas_omp_for( size_t n=0; n<size; ++n )
    {
      if( values[n] > local_maximum_private )
      {
        local_maximum_private = values[n];
        location_private = n;
      }
    }
    atlas_omp_critical_ordered
    {
      if( local_maximum_private > local_maximum )
      {
        local_maximum = local_maximum_private;
        location = location_private;
      }
    }
  }
  glb_idx = fs.nodes().global_index().data<gidx_t>()[location];
  ASSERT( glb_idx < std::numeric_limits<int>::max() ); // pairs with 64bit integer for second not implemented
  const std::pair<DATATYPE,int> max_and_gidx_loc (local_maximum,glb_idx);
  std::pair<DATATYPE,int> max_and_gidx_glb;
  eckit::mpi::all_reduce(max_and_gidx_loc,max_and_gidx_glb,eckit::mpi::maxloc());
  max     = max_and_gidx_glb.first;
  glb_idx = max_and_gidx_glb.second;
}

template< typename DATATYPE >
void maximum_and_location( const NodesFunctionSpace& fs, const Field& field, DATATYPE& max, gidx_t& glb_idx )
{
  if( field.kind() == DataType::kind<DATATYPE>() ) {
    return dispatch_maximum_and_location(fs,field,max,glb_idx);
  }
  else
  {
    switch( field.kind() )
    {
      case DataType::KIND_INT32 : {
        int tmp;
        dispatch_maximum_and_location(fs,field,tmp,glb_idx);
        max = tmp;
        return;
      }
      case DataType::KIND_INT64 : {
        long tmp;
        dispatch_maximum_and_location(fs,field,tmp,glb_idx);
        max = tmp;
        return;
      }
      case DataType::KIND_REAL32 : {
        float tmp;
        dispatch_maximum_and_location(fs,field,tmp,glb_idx);
        max = tmp;
        return;
      }
      case DataType::KIND_REAL64 : {
        double tmp;
        dispatch_maximum_and_location(fs,field,tmp,glb_idx);
        max = tmp;
        return;
      }
      default: throw eckit::Exception("datatype not supported",Here());
    }
  }
}


template< typename DATATYPE >
void dispatch_minimum_and_location( const NodesFunctionSpace& fs, const Field& field, std::vector<DATATYPE>& min, std::vector<gidx_t>& glb_idx )
{
  size_t nvar = field.stride(0);
  min.resize(nvar);
  glb_idx.resize(nvar);
  std::vector<DATATYPE> local_minimum(nvar,std::numeric_limits<DATATYPE>::max());
  std::vector<size_t> location(nvar);
  const ArrayView<DATATYPE,2> arr( field.data<DATATYPE>(), make_shape(field.shape(0),nvar) );
  atlas_omp_parallel
  {
    std::vector<DATATYPE> local_minimum_private(nvar,std::numeric_limits<DATATYPE>::max());
    std::vector<size_t> location_private(nvar);
    atlas_omp_for( size_t n=0; n<arr.shape(0); ++n ) {
      for( size_t j=0; j<arr.shape(1); ++j ) {
        if( arr(n,j) < local_minimum_private[j] ) {
          local_minimum_private[j] = arr(n,j);
          location_private[j] = n;
        }
      }
    }
    atlas_omp_critical_ordered
    {
      for( size_t j=0; j<arr.shape(1); ++j ) {
        if( local_minimum_private[j] < local_minimum[j] ) {
          local_minimum[j] = local_minimum_private[j];
          location[j] = location_private[j];
        }
      }
    }
  }
  std::vector< std::pair<DATATYPE,int> > min_and_gidx_loc(nvar);
  std::vector< std::pair<DATATYPE,int> > min_and_gidx_glb(nvar);
  const gidx_t* global_index = fs.nodes().global_index().data<gidx_t>();
  for( size_t j=0; j<nvar; ++j ) {
    gidx_t glb_idx = global_index[location[j]];
    ASSERT( glb_idx < std::numeric_limits<int>::max() ); // pairs with 64bit integer for second not implemented
    min_and_gidx_loc[j] = std::make_pair(local_minimum[j],glb_idx);
  }
  eckit::mpi::all_reduce(min_and_gidx_loc,min_and_gidx_glb,eckit::mpi::minloc());
  for( size_t j=0; j<nvar; ++j ) {
    min[j]     = min_and_gidx_glb[j].first;
    glb_idx[j] = min_and_gidx_glb[j].second;
  }
}

template< typename DATATYPE >
void minimum_and_location( const NodesFunctionSpace& fs, const Field& field, std::vector<DATATYPE>& min, std::vector<gidx_t>& glb_idx )
{
  if( field.kind() == DataType::kind<DATATYPE>() ) {
    return dispatch_minimum_and_location(fs,field,min,glb_idx);
  }
  else
  {
    switch( field.kind() )
    {
      case DataType::KIND_INT32 : {
        std::vector<int> tmp;
        dispatch_minimum_and_location(fs,field,tmp,glb_idx);
        min.assign(tmp.begin(),tmp.end());
        return;
      }
      case DataType::KIND_INT64 : {
        std::vector<long> tmp;
        dispatch_minimum_and_location(fs,field,tmp,glb_idx);
        min.assign(tmp.begin(),tmp.end());
        return;
      }
      case DataType::KIND_REAL32 : {
        std::vector<float> tmp;
        dispatch_minimum_and_location(fs,field,tmp,glb_idx);
        min.assign(tmp.begin(),tmp.end());
        return;
      }
      case DataType::KIND_REAL64 : {
        std::vector<double> tmp;
        dispatch_minimum_and_location(fs,field,tmp,glb_idx);
        min.assign(tmp.begin(),tmp.end());
        return;
      }
      default: throw eckit::Exception("datatype not supported",Here());
    }
  }
}


template< typename DATATYPE >
void dispatch_maximum_and_location( const NodesFunctionSpace& fs, const Field& field, std::vector<DATATYPE>& max, std::vector<gidx_t>& glb_idx )
{
  size_t nvar = field.stride(0);
  max.resize(nvar);
  glb_idx.resize(nvar);
  std::vector<DATATYPE> local_maximum(nvar,-std::numeric_limits<DATATYPE>::max());
  std::vector<size_t> location(nvar);
  ArrayView<DATATYPE,2> arr( field.data<DATATYPE>(), make_shape(field.shape(0),nvar) );
  atlas_omp_parallel
  {
    std::vector<DATATYPE> local_maximum_private(nvar,-std::numeric_limits<DATATYPE>::max());
    std::vector<size_t> location_private(nvar);
    atlas_omp_for( size_t n=0; n<arr.shape(0); ++n ) {
      for( size_t j=0; j<arr.shape(1); ++j ) {
        if( arr(n,j) > local_maximum_private[j] ) {
          local_maximum_private[j] = arr(n,j);
          location_private[j] = n;
        }
      }
    }
    atlas_omp_critical_ordered
    {
      for( size_t j=0; j<arr.shape(1); ++j ) {
        if( local_maximum_private[j] > local_maximum[j] ) {
          local_maximum[j] = local_maximum_private[j];
          location[j] = location_private[j];
        }
      }
    }
  }
  std::vector< std::pair<DATATYPE,int> > max_and_gidx_loc(nvar);
  std::vector< std::pair<DATATYPE,int> > max_and_gidx_glb(nvar);
  const gidx_t* global_index = fs.nodes().global_index().data<gidx_t>();
  for( size_t j=0; j<nvar; ++j ) {
    gidx_t glb_idx = global_index[location[j]];
    ASSERT( glb_idx < std::numeric_limits<int>::max() ); // pairs with 64bit integer for second not implemented
    max_and_gidx_loc[j] = std::make_pair(local_maximum[j],glb_idx);
  }
  eckit::mpi::all_reduce(max_and_gidx_loc,max_and_gidx_glb,eckit::mpi::maxloc());
  for( size_t j=0; j<nvar; ++j ) {
    max[j]     = max_and_gidx_glb[j].first;
    glb_idx[j] = max_and_gidx_glb[j].second;
  }
}

template< typename DATATYPE >
void maximum_and_location( const NodesFunctionSpace& fs, const Field& field, std::vector<DATATYPE>& max, std::vector<gidx_t>& glb_idx )
{
  if( field.kind() == DataType::kind<DATATYPE>() ) {
    return dispatch_maximum_and_location(fs,field,max,glb_idx);
  }
  else
  {
    switch( field.kind() )
    {
      case DataType::KIND_INT32 : {
        std::vector<int> tmp;
        dispatch_maximum_and_location(fs,field,tmp,glb_idx);
        max.assign(tmp.begin(),tmp.end());
        return;
      }
      case DataType::KIND_INT64 : {
        std::vector<long> tmp;
        dispatch_maximum_and_location(fs,field,tmp,glb_idx);
        max.assign(tmp.begin(),tmp.end());
        return;
      }
      case DataType::KIND_REAL32 : {
        std::vector<float> tmp;
        dispatch_maximum_and_location(fs,field,tmp,glb_idx);
        max.assign(tmp.begin(),tmp.end());
        return;
      }
      case DataType::KIND_REAL64 : {
        std::vector<double> tmp;
        dispatch_maximum_and_location(fs,field,tmp,glb_idx);
        max.assign(tmp.begin(),tmp.end());
        return;
      }
      default: throw eckit::Exception("datatype not supported",Here());
    }
  }
}



template< typename DATATYPE >
void mean( const NodesFunctionSpace& fs, const Field& field, DATATYPE& result, size_t& N )
{
  sum(fs,field,result,N);
  result /= static_cast<double>(N);
}

template< typename DATATYPE >
void mean( const NodesFunctionSpace& fs, const Field& field, std::vector<DATATYPE>& result, size_t& N )
{
  sum(fs,field,result,N);
  for( size_t j=0; j<result.size(); ++j ) {
    result[j] /= static_cast<double>(N);
  }
}

template< typename DATATYPE >
void mean_and_standard_deviation( const NodesFunctionSpace& fs, const Field& field, DATATYPE& mu, DATATYPE& sigma, size_t& N )
{
  mean(fs,field,mu,N);
  Field::Ptr squared_diff_field( fs.createField(field) );
  ArrayView<DATATYPE,1> squared_diff( *squared_diff_field );
  ArrayView<DATATYPE,1> values( field );
  atlas_omp_parallel_for( size_t n=0; n<field.size(); ++n )
    squared_diff(n) = sqr( values(n) - mu );

  mean(fs,*squared_diff_field,sigma,N);
  sigma = std::sqrt(sigma);
}

template< typename DATATYPE >
void mean_and_standard_deviation( const NodesFunctionSpace& fs, const Field& field, std::vector<DATATYPE>& mu, std::vector<DATATYPE>& sigma, size_t& N )
{
  mean(fs,field,mu,N);
  Field::Ptr squared_diff_field( fs.createField(field) );
  ArrayView<DATATYPE,2> squared_diff( squared_diff_field->data<DATATYPE>(), make_shape(squared_diff_field->shape(0),squared_diff_field->stride(0)) );
  ArrayView<DATATYPE,2> values( field.data<DATATYPE>(), make_shape(field.shape(0),field.stride(0)) );
  atlas_omp_parallel_for( size_t n=0; n<values.shape(0); ++n ) {
    for( size_t j=0; j<values.shape(1); ++j ) {
      squared_diff(n,j) = sqr( values(n,j) - mu[j] );
    }
  }
  mean(fs,*squared_diff_field,sigma,N);
  for( size_t j=0; j<sigma.size(); ++j ) {
    sigma[j] = std::sqrt(sigma[j]);
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

template<> void NodesFunctionSpace::mean( const Field& field, float&  result, size_t& N ) const { return detail::mean(*this,field,result,N); }
template<> void NodesFunctionSpace::mean( const Field& field, double& result, size_t& N ) const { return detail::mean(*this,field,result,N); }

template<> void NodesFunctionSpace::mean( const Field& field, std::vector<float>&  result, size_t& N ) const { return detail::mean(*this,field,result,N); }
template<> void NodesFunctionSpace::mean( const Field& field, std::vector<double>& result, size_t& N ) const { return detail::mean(*this,field,result,N); }

template<> void NodesFunctionSpace::meanAndStandardDeviation( const Field& field, float&  mu, float&  sigma, size_t& N ) const { return detail::mean_and_standard_deviation(*this,field,mu,sigma,N); }
template<> void NodesFunctionSpace::meanAndStandardDeviation( const Field& field, double& mu, double& sigma, size_t& N ) const { return detail::mean_and_standard_deviation(*this,field,mu,sigma,N); }

template<> void NodesFunctionSpace::meanAndStandardDeviation( const Field& field, std::vector<float>&  mu, std::vector<float>&  sigma, size_t& N ) const { return detail::mean_and_standard_deviation(*this,field,mu,sigma,N); }
template<> void NodesFunctionSpace::meanAndStandardDeviation( const Field& field, std::vector<double>& mu, std::vector<double>& sigma, size_t& N ) const { return detail::mean_and_standard_deviation(*this,field,mu,sigma,N); }


// --------------------------------------------------------------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------------------------------------------------------------

NodesColumnFunctionSpace::NodesColumnFunctionSpace(const std::string& name, Mesh& mesh, size_t nb_levels, const Halo& halo)
  : NodesFunctionSpace(name,mesh,halo),
    nb_levels_(nb_levels)
{
}

NodesColumnFunctionSpace::~NodesColumnFunctionSpace() {}

size_t NodesColumnFunctionSpace::nb_levels() const
{
  return nb_levels_;
}

Field* NodesColumnFunctionSpace::createField(const std::string& datatype) const {
  return Field::create(make_shape(nb_nodes(),nb_levels()),datatype);
}

Field* NodesColumnFunctionSpace::createField(const std::string& name, const std::string& datatype) const {
  return Field::create(name,make_shape(nb_nodes(),nb_levels()),datatype);
}

Field* NodesColumnFunctionSpace::createField(const std::vector<size_t>& variables, const std::string& datatype) const {
  std::vector<size_t> shape(1,nb_nodes()); shape.push_back(nb_levels());
  for( size_t i=0; i<variables.size(); ++i ) shape.push_back(variables[i]);
  return Field::create(shape,datatype);
}

Field* NodesColumnFunctionSpace::createField(const std::string& name, const std::vector<size_t>& variables, const std::string& datatype) const {
  std::vector<size_t> shape(1,nb_nodes()); shape.push_back(nb_levels());
  for( size_t i=0; i<variables.size(); ++i ) shape.push_back(variables[i]);
  return Field::create(name, shape,datatype);
}

Field* NodesColumnFunctionSpace::createField(const Field& other) const {
  ArrayShape shape = other.shape();
  shape[0] = nb_nodes();
  shape[1] = nb_levels();
  return Field::create(shape,other.datatype());
}

Field* NodesColumnFunctionSpace::createField(const std::string& name, const Field& other) const {
  ArrayShape shape = other.shape();
  shape[0] = nb_nodes();
  shape[1] = nb_levels();
  return Field::create(name,shape,other.datatype());
}

Field* NodesColumnFunctionSpace::createGlobalField(const std::string& datatype) const {
  return Field::create(make_shape(nb_nodes_global(),nb_levels()),datatype);
}

Field* NodesColumnFunctionSpace::createGlobalField(const std::string& name, const std::string& datatype) const {
  return Field::create(name,make_shape(nb_nodes_global(),nb_levels()),datatype);
}

Field* NodesColumnFunctionSpace::createGlobalField(const std::vector<size_t>& variables, const std::string& datatype) const {
  std::vector<size_t> shape(1,nb_nodes_global()); shape.push_back(nb_levels());
  for( size_t i=0; i<variables.size(); ++i ) shape.push_back(variables[i]);
  return Field::create(shape,datatype);
}

Field* NodesColumnFunctionSpace::createGlobalField(const std::string& name, const std::vector<size_t>& variables, const std::string& datatype) const {
  std::vector<size_t> shape(1,nb_nodes_global()); shape.push_back(nb_levels());
  for( size_t i=0; i<variables.size(); ++i ) shape.push_back(variables[i]);
  return Field::create(name, shape,datatype);
}

Field* NodesColumnFunctionSpace::createGlobalField(const Field& other) const {
  ArrayShape shape = other.shape();
  shape[0] = nb_nodes_global();
  shape[1] = nb_levels();
  return Field::create(shape,other.datatype());
}

Field* NodesColumnFunctionSpace::createGlobalField(const std::string& name, const Field& other) const {
  ArrayShape shape = other.shape();
  shape[0] = nb_nodes_global();
  shape[1] = nb_levels();
  return Field::create(name,shape,other.datatype());
}

namespace {

template <typename DATATYPE>
std::string checksum_3d_field(const mpl::Checksum& checksum, const Field& field )
{
  ArrayShape shape(1,field.shape(0));
  for( size_t j=2; j<field.rank(); ++j )
    shape.push_back(field.shape(j));
  Field::Ptr surface_field ( Field::create<DATATYPE>(shape) );
  ArrayView<DATATYPE,2> surface( surface_field->data<DATATYPE>(), make_shape(surface_field->shape(0), surface_field->stride(0)) );
  ArrayView<DATATYPE,3> values( field.data<DATATYPE>(), make_shape(field.shape(0),field.shape(1),field.stride(1) ) );
  for( size_t n=0; n<field.shape(0); ++n ) {
    for( size_t j=0; j<surface.shape(1); ++j )
    {
      surface(n,j) = 0.;
      for( size_t l=0; l<field.shape(1);++l )
        surface(n,j) += values(n,l,j);
    }
  }
  return checksum.execute( surface.data(), surface.stride(0) );
}
}

std::string NodesColumnFunctionSpace::checksum( const FieldSet& fieldset ) const {
  const mpl::Checksum& checksum = mesh().checksum().get(checksum_name());
  eckit::MD5 md5;
  for( size_t f=0; f<fieldset.size(); ++f ) {
    const Field& field=fieldset[f];
    if     ( field.kind() == DataType::kind<int>() )
      md5 << checksum_3d_field<int>(checksum,field);
    else if( field.kind() == DataType::kind<long>() )
      md5 << checksum_3d_field<long>(checksum,field);
    else if( field.kind() == DataType::kind<float>() )
      md5 << checksum_3d_field<float>(checksum,field);
    else if( field.kind() == DataType::kind<double>() )
      md5 << checksum_3d_field<double>(checksum,field);
    else throw eckit::Exception("datatype not supported",Here());
  }
  return md5;
}
std::string NodesColumnFunctionSpace::checksum( const Field& field ) const {
  FieldSet fieldset;
  fieldset.add(field);
  return checksum(fieldset);
}


namespace detail { // Collectives implementation

template< typename DATATYPE >
void dispatch_sum( const NodesColumnFunctionSpace& fs, const Field& field, DATATYPE& result, size_t& N )
{
  const util::IsGhost is_ghost(fs.nodes());
  const ArrayView<DATATYPE,2> arr( field );
  DATATYPE local_sum = 0;
  atlas_omp_pragma( omp parallel for default(shared) reduction(+:local_sum) )
  for( size_t n=0; n<arr.shape(0); ++n ) {
    if( ! is_ghost(n) ) {
      for( size_t l=0; l<arr.shape(1); ++l )
        local_sum += arr(n,l);
    }
  }
  eckit::mpi::all_reduce(local_sum,result,eckit::mpi::sum());
  N = field.shape(1) * fs.nb_nodes_global_broadcasted_;
}

template< typename DATATYPE >
void sum( const NodesColumnFunctionSpace& fs , const Field& field, DATATYPE& result, size_t& N )
{
  if( field.kind() == DataType::kind<DATATYPE>() ) {
    return dispatch_sum(fs,field,result,N);
  }
  else
  {
    switch( field.kind() )
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



template< typename DATATYPE >
void dispatch_sum( const NodesColumnFunctionSpace& fs, const Field& field, std::vector<DATATYPE>& result, size_t& N )
{
  const size_t nvar = field.stride(1);
  const ArrayView<DATATYPE,3> arr( field.data<DATATYPE>(), make_shape(field.shape(0),field.shape(1),nvar) );
  const util::IsGhost is_ghost(fs.nodes());
  std::vector<DATATYPE> local_sum(nvar,0);

  atlas_omp_parallel
  {
    std::vector<DATATYPE> local_sum_private(nvar,0);
    atlas_omp_for( size_t n=0; n<arr.shape(0); ++n )
    {
      if( ! is_ghost(n) ) {
        for( size_t l=0; l<arr.shape(1); ++l ) {
          for( size_t j=0; j<nvar; ++j ) {
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
  N = fs.nb_nodes_global_broadcasted_*field.shape(1);
}

template< typename DATATYPE >
void sum( const NodesColumnFunctionSpace& fs , const Field& field, std::vector<DATATYPE>& result, size_t& N )
{
  if( field.kind() == DataType::kind<DATATYPE>() ) {
    return dispatch_sum(fs,field,result,N);
  }
  else
  {
    switch( field.kind() )
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


template< typename DATATYPE >
void dispatch_sum_per_level( const NodesColumnFunctionSpace& fs, const Field& field, Field& sum, size_t& N )
{
  ArrayShape shape;
  shape.reserve(field.rank()-1);
  for( size_t j=1; j<field.rank(); ++j )
    shape.push_back(field.shape(j));
  sum.resize(shape);

  const size_t nvar = field.stride(1);
  const ArrayView<DATATYPE,3> arr( field.data<DATATYPE>(), make_shape(field.shape(0),field.shape(1),nvar) );
  ArrayView<DATATYPE,2> sum_per_level( sum.data<DATATYPE>(), make_shape(sum.shape(0),sum.stride(0)) );
  sum_per_level = 0;
  const util::IsGhost is_ghost(fs.nodes());

  atlas_omp_parallel
  {
    Array<DATATYPE> sum_per_level_private(sum_per_level.shape(0),sum_per_level.shape(1));
    ArrayView<DATATYPE> sum_per_level_private_view(sum_per_level_private); sum_per_level_private_view = 0.;
    atlas_omp_for( size_t n=0; n<arr.shape(0); ++n )
    {
      if( ! is_ghost(n) ) {
        for( size_t l=0; l<arr.shape(1); ++l ) {
          for( size_t j=0; j<nvar; ++j ) {
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
  N = fs.nb_nodes_global_broadcasted_;
}

void sum_per_level( const NodesColumnFunctionSpace& fs, const Field& field, Field& sum, size_t& N )
{
  if( field.kind() != sum.kind() ) {
    throw eckit::Exception("Field and sum are not of same datatype.",Here());
  }
  switch( field.kind() )
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
void dispatch_order_independent_sum( const NodesColumnFunctionSpace& fs , const Field& field, DATATYPE& result, size_t& N )
{
  const ArrayView<DATATYPE,2> arr( field.data<DATATYPE>(), make_shape(field.shape(0),field.shape(1)) );

  Field::Ptr surface_field( dynamic_cast<const NodesFunctionSpace&>(fs).createField<DATATYPE>() );
  ArrayView<DATATYPE,1> surface( *surface_field );

  for( size_t n=0; n<field.shape(0); ++n ) {
    for( size_t l=0; l<field.shape(1); ++l ) {
        surface(n) += arr(n,l);
    }
  }

  dispatch_order_independent_sum( dynamic_cast<const NodesFunctionSpace&>(fs), *surface_field, result, N );
  N *= field.shape(1);
}

template< typename DATATYPE >
void order_independent_sum( const NodesColumnFunctionSpace& fs , const Field& field, DATATYPE& result, size_t& N )
{
  if( field.kind() == DataType::kind<DATATYPE>() ) {
    return dispatch_order_independent_sum(fs,field,result,N);
  }
  else
  {
    switch( field.kind() )
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
void dispatch_order_independent_sum( const NodesColumnFunctionSpace& fs, const Field& field, std::vector<DATATYPE>& result, size_t& N )
{
  const size_t nvar = field.stride(1);
  const ArrayView<DATATYPE,3> arr( field.data<DATATYPE>(), make_shape(field.shape(0),field.shape(1),nvar) );

  Field::Ptr surface_field( dynamic_cast<const NodesFunctionSpace&>(fs).createField<DATATYPE>(make_shape(nvar)) );
  ArrayView<DATATYPE,2> surface( *surface_field );

  for( size_t n=0; n<field.shape(0); ++n ) {
    for( size_t l=0; l<field.shape(1); ++l ) {
      for( size_t j=0; j<nvar; ++j ) {
        surface(n,j) += arr(n,l,j);
      }
    }
  }

  dispatch_order_independent_sum( dynamic_cast<const NodesFunctionSpace&>(fs), *surface_field, result, N );
  N *= field.shape(1);
}

template< typename DATATYPE >
void order_independent_sum( const NodesColumnFunctionSpace& fs, const Field& field, std::vector<DATATYPE>& result, size_t& N )
{
  if( field.kind() == DataType::kind<DATATYPE>() ) {
    return dispatch_order_independent_sum(fs,field,result,N);
  }
  else
  {
    switch( field.kind() )
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

template< typename DATATYPE >
void dispatch_order_independent_sum_per_level( const NodesColumnFunctionSpace& fs, const Field& field, Field& sumfield, size_t& N )
{
  ArrayShape shape;
  shape.reserve(field.rank()-1);
  for( size_t j=1; j<field.rank(); ++j )
    shape.push_back(field.shape(j));
  sumfield.resize(shape);

  ArrayView<DATATYPE,2> sum ( sumfield.data<DATATYPE>(), make_shape(sumfield.shape(0),sumfield.stride(0)) );
  sum = 0.;

  size_t root = 0;
  Field::Ptr global( fs.createGlobalField(field) );
  fs.gather(field,*global);
  if( eckit::mpi::rank() == 0 ) {
    const ArrayView<DATATYPE,3> glb( global->data<DATATYPE>(), make_shape(global->shape(0),global->shape(1),global->stride(1)) );
    for( size_t n=0; n<glb.shape(0); ++n ) {
      for( size_t l=0; l<glb.shape(1); ++l ) {
        for( size_t j=0; j<glb.shape(2); ++j ) {
          sum(l,j) += glb(n,l,j);
        }
      }
    }
  }
  eckit::mpi::broadcast(sumfield.data<DATATYPE>(),sumfield.size(),root);
  N = fs.nb_nodes_global_broadcasted_;
}

void order_independent_sum_per_level( const NodesColumnFunctionSpace& fs, const Field& field, Field& sum, size_t& N )
{
  if( field.kind() != sum.kind() ) {
    throw eckit::Exception("Field and sum are not of same datatype.",Here());
  }
  switch( field.kind() )
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

template< typename DATATYPE >
void dispatch_minimum( const NodesColumnFunctionSpace& fs, const Field& field, std::vector<DATATYPE>& min )
{
  const size_t nvar = field.stride(1);
  min.resize(nvar);
  std::vector<DATATYPE> local_minimum(nvar,std::numeric_limits<DATATYPE>::max());
  const ArrayView<DATATYPE,3> arr( field.data<DATATYPE>(), make_shape(field.shape(0),field.shape(1),nvar) );
  atlas_omp_parallel
  {
    std::vector<DATATYPE> local_minimum_private(nvar,std::numeric_limits<DATATYPE>::max());
    atlas_omp_for( size_t n=0; n<arr.shape(0); ++n ) {
      for( size_t l=0; l<arr.shape(1); ++l ) {
        for( size_t j=0; j<nvar; ++j ) {
          local_minimum_private[j] = std::min(arr(n,l,j),local_minimum_private[j]);
        }
      }
    }
    atlas_omp_critical
    {
      for( size_t j=0; j<nvar; ++j ) {
        local_minimum[j] = std::min(local_minimum_private[j],local_minimum[j]);
      }
    }
  }
  eckit::mpi::all_reduce(local_minimum,min,eckit::mpi::min());
}

template< typename DATATYPE >
void minimum( const NodesColumnFunctionSpace& fs, const Field& field, std::vector<DATATYPE>& min )
{
  if( field.kind() == DataType::kind<DATATYPE>() ) {
    return dispatch_minimum(fs,field,min);
  }
  else
  {
    switch( field.kind() )
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

template< typename DATATYPE >
void dispatch_maximum( const NodesColumnFunctionSpace& fs, const Field& field, std::vector<DATATYPE>& max )
{
  const size_t nvar = field.stride(1);
  max.resize(nvar);
  std::vector<DATATYPE> local_maximum(nvar,-std::numeric_limits<DATATYPE>::max());
  const ArrayView<DATATYPE,3> arr( field.data<DATATYPE>(), make_shape(field.shape(0),field.shape(1),nvar) );
  atlas_omp_parallel
  {
    std::vector<DATATYPE> local_maximum_private(nvar,-std::numeric_limits<DATATYPE>::max());
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

template< typename DATATYPE >
void maximum( const NodesColumnFunctionSpace& fs, const Field& field, std::vector<DATATYPE>& max )
{
  if( field.kind() == DataType::kind<DATATYPE>() ) {
    return dispatch_maximum(fs,field,max);
  }
  else
  {
    switch( field.kind() )
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



template< typename DATATYPE >
void dispatch_minimum_per_level( const NodesColumnFunctionSpace& fs, const Field& field, Field& min_field )
{
  ArrayShape shape;
  shape.reserve(field.rank()-1);
  for( size_t j=1; j<field.rank(); ++j )
    shape.push_back(field.shape(j));
  min_field.resize(shape);
  const size_t nvar = field.stride(1);
  ArrayView<DATATYPE,2> min( min_field.data<DATATYPE>(), make_shape(min_field.shape(0),min_field.stride(0)) );
  min = std::numeric_limits<DATATYPE>::max();
  const ArrayView<DATATYPE,3> arr( field.data<DATATYPE>(), make_shape(field.shape(0),field.shape(1),nvar) );
  atlas_omp_parallel
  {
    Array<DATATYPE> min_private(min.shape(0),min.shape(1));
    ArrayView<DATATYPE> min_private_view(min_private); min_private_view = std::numeric_limits<DATATYPE>::max();
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
  eckit::mpi::all_reduce(min_field.data<DATATYPE>(),min_field.size(),eckit::mpi::min());
}

void minimum_per_level( const NodesColumnFunctionSpace& fs, const Field& field, Field& min )
{
  if( field.kind() != min.kind() ) {
    throw eckit::Exception("Field and min are not of same datatype.",Here());
  }
  switch( field.kind() )
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

template< typename DATATYPE >
void dispatch_maximum_per_level( const NodesColumnFunctionSpace& fs, const Field& field, Field& max_field )
{
  ArrayShape shape;
  shape.reserve(field.rank()-1);
  for( size_t j=1; j<field.rank(); ++j )
    shape.push_back(field.shape(j));
  max_field.resize(shape);
  const size_t nvar = field.stride(1);
  ArrayView<DATATYPE,2> max( max_field.data<DATATYPE>(), make_shape(max_field.shape(0),max_field.stride(0)) );
  max = -std::numeric_limits<DATATYPE>::max();
  const ArrayView<DATATYPE,3> arr( field.data<DATATYPE>(), make_shape(field.shape(0),field.shape(1),nvar) );
  atlas_omp_parallel
  {
    Array<DATATYPE> max_private(max.shape(0),max.shape(1));
    ArrayView<DATATYPE> max_private_view(max_private); max_private_view = -std::numeric_limits<DATATYPE>::max();
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
  eckit::mpi::all_reduce(max_field.data<DATATYPE>(),max_field.size(),eckit::mpi::max());
}

void maximum_per_level( const NodesColumnFunctionSpace& fs, const Field& field, Field& max )
{
  if( field.kind() != max.kind() ) {
    throw eckit::Exception("Field and max are not of same datatype.",Here());
  }
  switch( field.kind() )
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

template< typename DATATYPE >
void dispatch_minimum_and_location( const NodesColumnFunctionSpace& fs, const Field& field, DATATYPE& min, gidx_t& glb_idx, size_t& level )
{
  DATATYPE local_minimum(std::numeric_limits<DATATYPE>::max());
  size_t loc_node;
  size_t loc_level;
  ArrayView<DATATYPE,2> arr( field.data<DATATYPE>(), make_shape(field.shape(0),field.shape(1)) );
  atlas_omp_parallel
  {
    DATATYPE local_minimum_private(std::numeric_limits<DATATYPE>::max());
    size_t loc_node_private;
    size_t loc_level_private;
    atlas_omp_for( size_t n=0; n<arr.shape(0); ++n ) {
      for( size_t l=0; l<arr.shape(1); ++l ) {
        if( arr(n,l) < local_minimum_private ) {
          local_minimum_private = arr(n,l);
          loc_node_private = n;
          loc_level_private = l;
        }
      }
    }
    atlas_omp_critical_ordered
    {
      for( size_t l=0; l<arr.shape(1); ++l ) {
        if( local_minimum_private < local_minimum ) {
          local_minimum = local_minimum_private;
          loc_node = loc_node_private;
          loc_level = loc_level_private;
        }
      }
    }
  }
  std::pair<DATATYPE,int> min_and_gidx_loc;
  std::pair<DATATYPE,int> min_and_level_loc;
  std::pair<DATATYPE,int> min_and_gidx_glb;
  std::pair<DATATYPE,int> min_and_level_glb;
  const gidx_t* global_index = fs.nodes().global_index().data<gidx_t>();
  glb_idx = global_index[loc_node];
  ASSERT( glb_idx < std::numeric_limits<int>::max() ); // pairs with 64bit integer for second not implemented
  min_and_gidx_loc = std::make_pair(local_minimum,glb_idx);
  min_and_level_loc = std::make_pair(local_minimum,loc_level);
  eckit::mpi::all_reduce(min_and_gidx_loc, min_and_gidx_glb, eckit::mpi::minloc());
  eckit::mpi::all_reduce(min_and_level_loc,min_and_level_glb,eckit::mpi::minloc());
  min     = min_and_gidx_glb.first;
  glb_idx = min_and_gidx_glb.second;
  level   = min_and_level_glb.second;
}

template< typename DATATYPE >
void minimum_and_location( const NodesColumnFunctionSpace& fs, const Field& field, DATATYPE& min, gidx_t& glb_idx, size_t& level )
{
  if( field.kind() == DataType::kind<DATATYPE>() ) {
    return dispatch_minimum_and_location(fs,field,min,glb_idx,level);
  }
  else
  {
    switch( field.kind() )
    {
      case DataType::KIND_INT32 : {
        int tmp;
        dispatch_minimum_and_location(fs,field,tmp,glb_idx,level);
        min = tmp;
        return;
      }
      case DataType::KIND_INT64 : {
        long tmp;
        dispatch_minimum_and_location(fs,field,tmp,glb_idx,level);
        min = tmp;
        return;
      }
      case DataType::KIND_REAL32 : {
        float tmp;
        dispatch_minimum_and_location(fs,field,tmp,glb_idx,level);
        min = tmp;
        return;
      }
      case DataType::KIND_REAL64 : {
        double tmp;
        dispatch_minimum_and_location(fs,field,tmp,glb_idx,level);
        min = tmp;
        return;
      }
      default: throw eckit::Exception("datatype not supported",Here());
    }
  }
}

template< typename DATATYPE >
void dispatch_maximum_and_location( const NodesColumnFunctionSpace& fs, const Field& field, DATATYPE& max, gidx_t& glb_idx, size_t& level )
{
  DATATYPE local_maximum(-std::numeric_limits<DATATYPE>::max());
  size_t loc_node;
  size_t loc_level;
  ArrayView<DATATYPE,2> arr( field.data<DATATYPE>(), make_shape(field.shape(0),field.shape(1)) );
  atlas_omp_parallel
  {
    DATATYPE local_maximum_private(-std::numeric_limits<DATATYPE>::max());
    size_t loc_node_private;
    size_t loc_level_private;
    atlas_omp_for( size_t n=0; n<arr.shape(0); ++n ) {
      for( size_t l=0; l<arr.shape(1); ++l ) {
        if( arr(n,l) > local_maximum ) {
          local_maximum = arr(n,l);
          loc_node = n;
          loc_level = l;
        }
      }
    }
    atlas_omp_critical_ordered
    {
      for( size_t l=0; l<arr.shape(1); ++l ) {
        if( local_maximum_private > local_maximum ) {
          local_maximum = local_maximum_private;
          loc_node = loc_node_private;
          loc_level = loc_level_private;
        }
      }

    }
  }
  std::pair<DATATYPE,int> max_and_gidx_loc;
  std::pair<DATATYPE,int> max_and_level_loc;
  std::pair<DATATYPE,int> max_and_gidx_glb;
  std::pair<DATATYPE,int> max_and_level_glb;
  const gidx_t* global_index = fs.nodes().global_index().data<gidx_t>();
  glb_idx = global_index[loc_node];
  ASSERT( glb_idx < std::numeric_limits<int>::max() ); // pairs with 64bit integer for second not implemented
  max_and_gidx_loc = std::make_pair(local_maximum,glb_idx);
  max_and_level_loc = std::make_pair(local_maximum,loc_level);
  eckit::mpi::all_reduce(max_and_gidx_loc, max_and_gidx_glb, eckit::mpi::maxloc());
  eckit::mpi::all_reduce(max_and_level_loc,max_and_level_glb,eckit::mpi::maxloc());
  max     = max_and_gidx_glb.first;
  glb_idx = max_and_gidx_glb.second;
  level   = max_and_level_glb.second;
}

template< typename DATATYPE >
void maximum_and_location( const NodesColumnFunctionSpace& fs, const Field& field, DATATYPE& max, gidx_t& glb_idx, size_t& level )
{
  if( field.kind() == DataType::kind<DATATYPE>() ) {
    return dispatch_maximum_and_location(fs,field,max,glb_idx,level);
  }
  else
  {
    switch( field.kind() )
    {
      case DataType::KIND_INT32 : {
        int tmp;
        dispatch_maximum_and_location(fs,field,tmp,glb_idx,level);
        max = tmp;
        return;
      }
      case DataType::KIND_INT64 : {
        long tmp;
        dispatch_maximum_and_location(fs,field,tmp,glb_idx,level);
        max = tmp;
        return;
      }
      case DataType::KIND_REAL32 : {
        float tmp;
        dispatch_maximum_and_location(fs,field,tmp,glb_idx,level);
        max = tmp;
        return;
      }
      case DataType::KIND_REAL64 : {
        double tmp;
        dispatch_maximum_and_location(fs,field,tmp,glb_idx,level);
        max = tmp;
        return;
      }
      default: throw eckit::Exception("datatype not supported",Here());
    }
  }
}


template< typename DATATYPE >
void dispatch_minimum_and_location( const NodesColumnFunctionSpace& fs, const Field& field, std::vector<DATATYPE>& min, std::vector<gidx_t>& glb_idx, std::vector<size_t>& level )
{
  size_t nvar = field.stride(1);
  min.resize(nvar);
  glb_idx.resize(nvar);
  level.resize(nvar);
  std::vector<DATATYPE> local_minimum(nvar,std::numeric_limits<DATATYPE>::max());
  std::vector<size_t> loc_node(nvar);
  std::vector<size_t> loc_level(nvar);
  ArrayView<DATATYPE,3> arr( field.data<DATATYPE>(), make_shape(field.shape(0),field.shape(1),nvar) );
  atlas_omp_parallel
  {
    std::vector<DATATYPE> local_minimum_private(nvar,std::numeric_limits<DATATYPE>::max());
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
  std::vector< std::pair<DATATYPE,int> > min_and_gidx_loc(nvar);
  std::vector< std::pair<DATATYPE,int> > min_and_level_loc(nvar);
  std::vector< std::pair<DATATYPE,int> > min_and_gidx_glb(nvar);
  std::vector< std::pair<DATATYPE,int> > min_and_level_glb(nvar);
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

template< typename DATATYPE >
void minimum_and_location( const NodesColumnFunctionSpace& fs, const Field& field, std::vector<DATATYPE>& min, std::vector<gidx_t>& glb_idx, std::vector<size_t>& level )
{
  if( field.kind() == DataType::kind<DATATYPE>() ) {
    return dispatch_minimum_and_location(fs,field,min,glb_idx,level);
  }
  else
  {
    switch( field.kind() )
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



template< typename DATATYPE >
void dispatch_maximum_and_location( const NodesColumnFunctionSpace& fs, const Field& field, std::vector<DATATYPE>& max, std::vector<gidx_t>& glb_idx, std::vector<size_t>& level )
{
  size_t nvar = field.stride(1);
  max.resize(nvar);
  glb_idx.resize(nvar);
  level.resize(nvar);
  std::vector<DATATYPE> local_maximum(nvar,-std::numeric_limits<DATATYPE>::max());
  std::vector<size_t> loc_node(nvar);
  std::vector<size_t> loc_level(nvar);
  ArrayView<DATATYPE,3> arr( field.data<DATATYPE>(), make_shape(field.shape(0),field.shape(1),nvar) );
  atlas_omp_parallel
  {
    std::vector<DATATYPE> local_maximum_private(nvar,-std::numeric_limits<DATATYPE>::max());
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
  std::vector< std::pair<DATATYPE,int> > max_and_gidx_loc(nvar);
  std::vector< std::pair<DATATYPE,int> > max_and_level_loc(nvar);
  std::vector< std::pair<DATATYPE,int> > max_and_gidx_glb(nvar);
  std::vector< std::pair<DATATYPE,int> > max_and_level_glb(nvar);
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

template< typename DATATYPE >
void maximum_and_location( const NodesColumnFunctionSpace& fs, const Field& field, std::vector<DATATYPE>& max, std::vector<gidx_t>& glb_idx, std::vector<size_t>& level )
{
  if( field.kind() == DataType::kind<DATATYPE>() ) {
    return dispatch_maximum_and_location(fs,field,max,glb_idx,level);
  }
  else
  {
    switch( field.kind() )
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

template< typename DATATYPE >
void dispatch_minimum_and_location_per_level( const NodesColumnFunctionSpace& fs, const Field& field, Field& min_field, Field& glb_idx_field )
{
  ArrayShape shape;
  shape.reserve(field.rank()-1);
  for( size_t j=1; j<field.rank(); ++j )
    shape.push_back(field.shape(j));
  min_field.resize(shape);
  glb_idx_field.resize(shape);
  const size_t nvar = field.stride(1);
  ArrayView<DATATYPE,2> min( min_field.data<DATATYPE>(), make_shape(min_field.shape(0),min_field.stride(0)) );
  min = std::numeric_limits<DATATYPE>::max();
  ArrayView<gidx_t,2> glb_idx( glb_idx_field.data<gidx_t>(), make_shape(glb_idx_field.shape(0),glb_idx_field.stride(0)) );

  const ArrayView<DATATYPE,3> arr( field.data<DATATYPE>(), make_shape(field.shape(0),field.shape(1),nvar) );
  atlas_omp_parallel
  {
    Array<DATATYPE> min_private(min.shape(0),min.shape(1));
    ArrayView<DATATYPE> min_private_view(min_private); min_private_view = std::numeric_limits<DATATYPE>::max();
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
  const size_t nlev = field.shape(1);
  std::vector< std::pair<DATATYPE,int> > min_and_gidx_loc(nlev*nvar);
  std::vector< std::pair<DATATYPE,int> > min_and_gidx_glb(nlev*nvar);
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

void minimum_and_location_per_level( const NodesColumnFunctionSpace& fs, const Field& field, Field& min, Field& glb_idx )
{
  if( field.kind() != min.kind() ) {
    throw eckit::Exception("Field and min are not of same datatype.",Here());
  }
  if( glb_idx.kind() != DataType::kind<gidx_t>() ) {
    throw eckit::Exception("glb_idx Field is not of correct datatype",Here());
  }
  switch( field.kind() )
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


template< typename DATATYPE >
void dispatch_maximum_and_location_per_level( const NodesColumnFunctionSpace& fs, const Field& field, Field& max_field, Field& glb_idx_field )
{
  ArrayShape shape;
  shape.reserve(field.rank()-1);
  for( size_t j=1; j<field.rank(); ++j )
    shape.push_back(field.shape(j));
  max_field.resize(shape);
  glb_idx_field.resize(shape);
  const size_t nvar = field.stride(1);
  ArrayView<DATATYPE,2> max( max_field.data<DATATYPE>(), make_shape(max_field.shape(0),max_field.stride(0)) );
  max = -std::numeric_limits<DATATYPE>::max();
  ArrayView<gidx_t,2> glb_idx( glb_idx_field.data<gidx_t>(), make_shape(glb_idx_field.shape(0),glb_idx_field.stride(0)) );

  const ArrayView<DATATYPE,3> arr( field.data<DATATYPE>(), make_shape(field.shape(0),field.shape(1),nvar) );
  atlas_omp_parallel
  {
    Array<DATATYPE> max_private(max.shape(0),max.shape(1));
    ArrayView<DATATYPE> max_private_view(max_private); max_private_view = -std::numeric_limits<DATATYPE>::max();
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

  const size_t nlev = field.shape(1);
  std::vector< std::pair<DATATYPE,int> > max_and_gidx_loc(nlev*nvar);
  std::vector< std::pair<DATATYPE,int> > max_and_gidx_glb(nlev*nvar);
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

void maximum_and_location_per_level( const NodesColumnFunctionSpace& fs, const Field& field, Field& max, Field& glb_idx )
{
  if( field.kind() != max.kind() ) {
    throw eckit::Exception("Field and max are not of same datatype.",Here());
  }
  if( glb_idx.kind() != DataType::kind<gidx_t>() ) {
    throw eckit::Exception("glb_idx Field is not of correct datatype",Here());
  }
  switch( field.kind() )
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



template< typename DATATYPE >
void mean( const NodesColumnFunctionSpace& fs, const Field& field, DATATYPE& result, size_t& N )
{
  sum(fs,field,result,N);
  result /= static_cast<double>(N);
}

template< typename DATATYPE >
void mean( const NodesColumnFunctionSpace& fs, const Field& field, std::vector<DATATYPE>& result, size_t& N )
{
  sum(fs,field,result,N);
  for( size_t j=0; j<result.size(); ++j ) {
    result[j] /= static_cast<double>(N);
  }
}

template< typename DATATYPE >
void dispatch_mean_per_level( const NodesColumnFunctionSpace& fs, const Field& field, Field& mean, size_t& N )
{
  dispatch_sum_per_level<DATATYPE>(fs,field,mean,N);
  DATATYPE* rawdata = mean.data<DATATYPE>();
  for( size_t j=0; j<mean.size(); ++j ) {
    rawdata[j] /= static_cast<double>(N);
  }
}


void mean_per_level( const NodesColumnFunctionSpace& fs, const Field& field, Field& mean, size_t& N )
{
  if( field.kind() != mean.kind() ) {
    throw eckit::Exception("Field and sum are not of same datatype.",Here());
  }
  switch( field.kind() )
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

template< typename DATATYPE >
void mean_and_standard_deviation( const NodesColumnFunctionSpace& fs, const Field& field, DATATYPE& mu, DATATYPE& sigma, size_t& N )
{
  mean(fs,field,mu,N);
  Field::Ptr squared_diff_field( fs.createField(field) );
  ArrayView<DATATYPE,2> squared_diff( *squared_diff_field );
  ArrayView<DATATYPE,2> values( field );

  atlas_omp_parallel_for( size_t n=0; n<field.shape(0); ++n ) {
    for( size_t l=0; l<field.shape(1); ++l ) {
      squared_diff(n,l) = sqr( values(n,l) - mu );
    }
  }
  mean(fs,*squared_diff_field,sigma,N);
  sigma = std::sqrt(sigma);
}

template< typename DATATYPE >
void mean_and_standard_deviation( const NodesColumnFunctionSpace& fs, const Field& field, std::vector<DATATYPE>& mu, std::vector<DATATYPE>& sigma, size_t& N )
{
  mean(fs,field,mu,N);
  Field::Ptr squared_diff_field( fs.createField(field) );
  ArrayView<DATATYPE,3> squared_diff( squared_diff_field->data<DATATYPE>(), make_shape(squared_diff_field->shape(0),squared_diff_field->shape(1),squared_diff_field->stride(1)) );
  ArrayView<DATATYPE,3> values( field.data<DATATYPE>(), make_shape(field.shape(0),field.shape(1),field.stride(1)) );

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

template< typename DATATYPE >
void dispatch_mean_and_standard_deviation_per_level( const NodesColumnFunctionSpace& fs, const Field& field, Field& mean, Field& stddev, size_t& N )
{
  dispatch_mean_per_level<DATATYPE>(fs,field,mean,N);
  const size_t nvar = mean.stride(0);
  ArrayView<DATATYPE,2> mu( mean.data<DATATYPE>(), make_shape(mean.shape(0),nvar) );
  Field::Ptr squared_diff_field( fs.createField(field) );
  ArrayView<DATATYPE,3> squared_diff( squared_diff_field->data<DATATYPE>(), make_shape(squared_diff_field->shape(0),squared_diff_field->shape(1),squared_diff_field->stride(1)) );
  ArrayView<DATATYPE,3> values( field.data<DATATYPE>(), make_shape(field.shape(0),field.shape(1),field.stride(1)) );

  atlas_omp_parallel_for( size_t n=0; n<values.shape(0); ++n ) {
    for( size_t l=0; l<values.shape(1); ++l ) {
      for( size_t j=0; j<nvar; ++j ) {
        squared_diff(n,l,j) = sqr( values(n,l,j) - mu(l,j) );
      }
    }
  }
  dispatch_mean_per_level<DATATYPE>(fs,*squared_diff_field,stddev,N);
  DATATYPE* sigma = stddev.data<DATATYPE>();
  atlas_omp_parallel_for( size_t j=0; j<stddev.size(); ++j ) {
    sigma[j] = std::sqrt(sigma[j]);
  }
}


void mean_and_standard_deviation_per_level( const NodesColumnFunctionSpace& fs, const Field& field, Field& mean, Field& stddev, size_t& N )
{
  if( field.kind() != mean.kind() ) {
    throw eckit::Exception("Field and mean are not of same datatype.",Here());
  }
  if( field.kind() != stddev.kind() ) {
    throw eckit::Exception("Field and stddev are not of same datatype.",Here());
  }
  switch( field.kind() )
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

} // Collectives implementation

template<> void NodesColumnFunctionSpace::sum( const Field& field, int&    result, size_t& N ) const { return detail::sum(*this,field,result,N); }
template<> void NodesColumnFunctionSpace::sum( const Field& field, long&   result, size_t& N ) const { return detail::sum(*this,field,result,N); }
template<> void NodesColumnFunctionSpace::sum( const Field& field, float&  result, size_t& N ) const { return detail::sum(*this,field,result,N); }
template<> void NodesColumnFunctionSpace::sum( const Field& field, double& result, size_t& N ) const { return detail::sum(*this,field,result,N); }

template<> void NodesColumnFunctionSpace::sum( const Field& field, std::vector<int>&    result, size_t& N ) const { return detail::sum(*this,field,result,N); }
template<> void NodesColumnFunctionSpace::sum( const Field& field, std::vector<long>&   result, size_t& N ) const { return detail::sum(*this,field,result,N); }
template<> void NodesColumnFunctionSpace::sum( const Field& field, std::vector<float>&  result, size_t& N ) const { return detail::sum(*this,field,result,N); }
template<> void NodesColumnFunctionSpace::sum( const Field& field, std::vector<double>& result, size_t& N ) const { return detail::sum(*this,field,result,N); }

void NodesColumnFunctionSpace::sumPerLevel(const Field &field, Field &sum, size_t &N) const { return detail::sum_per_level(*this,field,sum,N); }

template<> void NodesColumnFunctionSpace::orderIndependentSum( const Field& field, int&    result, size_t& N ) const { return detail::order_independent_sum(*this,field,result,N); }
template<> void NodesColumnFunctionSpace::orderIndependentSum( const Field& field, long&   result, size_t& N ) const { return detail::order_independent_sum(*this,field,result,N); }
template<> void NodesColumnFunctionSpace::orderIndependentSum( const Field& field, float&  result, size_t& N ) const { return detail::order_independent_sum(*this,field,result,N); }
template<> void NodesColumnFunctionSpace::orderIndependentSum( const Field& field, double& result, size_t& N ) const { return detail::order_independent_sum(*this,field,result,N); }

template<> void NodesColumnFunctionSpace::orderIndependentSum( const Field& field, std::vector<int>&    result, size_t& N ) const { return detail::order_independent_sum(*this,field,result,N); }
template<> void NodesColumnFunctionSpace::orderIndependentSum( const Field& field, std::vector<long>&   result, size_t& N ) const { return detail::order_independent_sum(*this,field,result,N); }
template<> void NodesColumnFunctionSpace::orderIndependentSum( const Field& field, std::vector<float>&  result, size_t& N ) const { return detail::order_independent_sum(*this,field,result,N); }
template<> void NodesColumnFunctionSpace::orderIndependentSum( const Field& field, std::vector<double>& result, size_t& N ) const { return detail::order_independent_sum(*this,field,result,N); }

void NodesColumnFunctionSpace::orderIndependentSumPerLevel(const Field &field, Field &sum, size_t &N) const { return detail::order_independent_sum_per_level(*this,field,sum,N); }

template<> void NodesColumnFunctionSpace::minimum( const Field& field, int&    min ) const { return detail::minimum(*this,field,min); }
template<> void NodesColumnFunctionSpace::minimum( const Field& field, long&   min ) const { return detail::minimum(*this,field,min); }
template<> void NodesColumnFunctionSpace::minimum( const Field& field, float&  min ) const { return detail::minimum(*this,field,min); }
template<> void NodesColumnFunctionSpace::minimum( const Field& field, double& min ) const { return detail::minimum(*this,field,min); }

template<> void NodesColumnFunctionSpace::maximum( const Field& field, int&    max ) const { return detail::maximum(*this,field,max); }
template<> void NodesColumnFunctionSpace::maximum( const Field& field, long&   max ) const { return detail::maximum(*this,field,max); }
template<> void NodesColumnFunctionSpace::maximum( const Field& field, float&  max ) const { return detail::maximum(*this,field,max); }
template<> void NodesColumnFunctionSpace::maximum( const Field& field, double& max ) const { return detail::maximum(*this,field,max); }

template<> void NodesColumnFunctionSpace::minimum( const Field& field, std::vector<int>&    min ) const { return detail::minimum(*this,field,min); }
template<> void NodesColumnFunctionSpace::minimum( const Field& field, std::vector<long>&   min ) const { return detail::minimum(*this,field,min); }
template<> void NodesColumnFunctionSpace::minimum( const Field& field, std::vector<float>&  min ) const { return detail::minimum(*this,field,min); }
template<> void NodesColumnFunctionSpace::minimum( const Field& field, std::vector<double>& min ) const { return detail::minimum(*this,field,min); }

template<> void NodesColumnFunctionSpace::maximum( const Field& field, std::vector<int>&    max ) const { return detail::maximum(*this,field,max); }
template<> void NodesColumnFunctionSpace::maximum( const Field& field, std::vector<long>&   max ) const { return detail::maximum(*this,field,max); }
template<> void NodesColumnFunctionSpace::maximum( const Field& field, std::vector<float>&  max ) const { return detail::maximum(*this,field,max); }
template<> void NodesColumnFunctionSpace::maximum( const Field& field, std::vector<double>& max ) const { return detail::maximum(*this,field,max); }

void NodesColumnFunctionSpace::minimumPerLevel(const Field &field, Field &min) const { return detail::minimum_per_level(*this,field,min); }
void NodesColumnFunctionSpace::maximumPerLevel(const Field &field, Field &max) const { return detail::maximum_per_level(*this,field,max);}

template<> void NodesColumnFunctionSpace::minimumAndLocation( const Field& field, int&    min, gidx_t& glb_idx, size_t& level ) const { return detail::minimum_and_location(*this,field,min,glb_idx,level); }
template<> void NodesColumnFunctionSpace::minimumAndLocation( const Field& field, long&   min, gidx_t& glb_idx, size_t& level ) const { return detail::minimum_and_location(*this,field,min,glb_idx,level); }
template<> void NodesColumnFunctionSpace::minimumAndLocation( const Field& field, float&  min, gidx_t& glb_idx, size_t& level ) const { return detail::minimum_and_location(*this,field,min,glb_idx,level); }
template<> void NodesColumnFunctionSpace::minimumAndLocation( const Field& field, double& min, gidx_t& glb_idx, size_t& level ) const { return detail::minimum_and_location(*this,field,min,glb_idx,level); }

template<> void NodesColumnFunctionSpace::maximumAndLocation( const Field& field, int&    max, gidx_t& glb_idx, size_t& level ) const { return detail::maximum_and_location(*this,field,max,glb_idx,level); }
template<> void NodesColumnFunctionSpace::maximumAndLocation( const Field& field, long&   max, gidx_t& glb_idx, size_t& level ) const { return detail::maximum_and_location(*this,field,max,glb_idx,level); }
template<> void NodesColumnFunctionSpace::maximumAndLocation( const Field& field, float&  max, gidx_t& glb_idx, size_t& level ) const { return detail::maximum_and_location(*this,field,max,glb_idx,level); }
template<> void NodesColumnFunctionSpace::maximumAndLocation( const Field& field, double& max, gidx_t& glb_idx, size_t& level ) const { return detail::maximum_and_location(*this,field,max,glb_idx,level); }

template<> void NodesColumnFunctionSpace::minimumAndLocation( const Field& field, std::vector<int>&    min, std::vector<gidx_t>& glb_idx, std::vector<size_t>& level ) const { return detail::minimum_and_location(*this,field,min,glb_idx,level); }
template<> void NodesColumnFunctionSpace::minimumAndLocation( const Field& field, std::vector<long>&   min, std::vector<gidx_t>& glb_idx, std::vector<size_t>& level ) const { return detail::minimum_and_location(*this,field,min,glb_idx,level); }
template<> void NodesColumnFunctionSpace::minimumAndLocation( const Field& field, std::vector<float>&  min, std::vector<gidx_t>& glb_idx, std::vector<size_t>& level ) const { return detail::minimum_and_location(*this,field,min,glb_idx,level); }
template<> void NodesColumnFunctionSpace::minimumAndLocation( const Field& field, std::vector<double>& min, std::vector<gidx_t>& glb_idx, std::vector<size_t>& level ) const { return detail::minimum_and_location(*this,field,min,glb_idx,level); }

template<> void NodesColumnFunctionSpace::maximumAndLocation( const Field& field, std::vector<int>&    max, std::vector<gidx_t>& glb_idx, std::vector<size_t>& level ) const { return detail::maximum_and_location(*this,field,max,glb_idx,level); }
template<> void NodesColumnFunctionSpace::maximumAndLocation( const Field& field, std::vector<long>&   max, std::vector<gidx_t>& glb_idx, std::vector<size_t>& level ) const { return detail::maximum_and_location(*this,field,max,glb_idx,level); }
template<> void NodesColumnFunctionSpace::maximumAndLocation( const Field& field, std::vector<float>&  max, std::vector<gidx_t>& glb_idx, std::vector<size_t>& level ) const { return detail::maximum_and_location(*this,field,max,glb_idx,level); }
template<> void NodesColumnFunctionSpace::maximumAndLocation( const Field& field, std::vector<double>& max, std::vector<gidx_t>& glb_idx, std::vector<size_t>& level ) const { return detail::maximum_and_location(*this,field,max,glb_idx,level); }

void NodesColumnFunctionSpace::minimumAndLocationPerLevel(const Field &field, Field &min, Field &glb_idx) const { detail::minimum_and_location_per_level(*this,field,min,glb_idx); }
void NodesColumnFunctionSpace::maximumAndLocationPerLevel(const Field &field, Field &max, Field &glb_idx) const { detail::maximum_and_location_per_level(*this,field,max,glb_idx); }

template<> void NodesColumnFunctionSpace::mean( const Field& field, float&  result, size_t& N ) const { return detail::mean(*this,field,result,N); }
template<> void NodesColumnFunctionSpace::mean( const Field& field, double& result, size_t& N ) const { return detail::mean(*this,field,result,N); }

template<> void NodesColumnFunctionSpace::mean( const Field& field, std::vector<float>&  result, size_t& N ) const { return detail::mean(*this,field,result,N); }
template<> void NodesColumnFunctionSpace::mean( const Field& field, std::vector<double>& result, size_t& N ) const { return detail::mean(*this,field,result,N); }

void NodesColumnFunctionSpace::meanPerLevel(const Field &field, Field &mean, size_t &N) const { return detail::mean_per_level(*this,field,mean,N); }

template<> void NodesColumnFunctionSpace::meanAndStandardDeviation( const Field& field, float&  mu, float&  sigma, size_t& N ) const { return detail::mean_and_standard_deviation(*this,field,mu,sigma,N); }
template<> void NodesColumnFunctionSpace::meanAndStandardDeviation( const Field& field, double& mu, double& sigma, size_t& N ) const { return detail::mean_and_standard_deviation(*this,field,mu,sigma,N); }

template<> void NodesColumnFunctionSpace::meanAndStandardDeviation( const Field& field, std::vector<float>&  mu, std::vector<float>&  sigma, size_t& N ) const { return detail::mean_and_standard_deviation(*this,field,mu,sigma,N); }
template<> void NodesColumnFunctionSpace::meanAndStandardDeviation( const Field& field, std::vector<double>& mu, std::vector<double>& sigma, size_t& N ) const { return detail::mean_and_standard_deviation(*this,field,mu,sigma,N); }

void NodesColumnFunctionSpace::meanAndStandardDeviationPerLevel(const Field &field, Field &mean, Field &stddev, size_t &N) const { return detail::mean_and_standard_deviation_per_level(*this,field,mean,stddev,N); }


// ----------------------------------------------------------------------

extern "C" {
NodesFunctionSpace* atlas__NodesFunctionSpace__new (const char* name, Mesh* mesh, int halo)
{
  ASSERT(mesh);
  return new NodesFunctionSpace(std::string(name),*mesh,Halo(halo));
}

void atlas__NodesFunctionSpace__delete (NodesFunctionSpace* This)
{
  ASSERT(This);
  delete(This);
}

int atlas__NodesFunctionSpace__nb_nodes(const NodesFunctionSpace* This)
{
  ASSERT(This);
  return This->nb_nodes();
}

Mesh* atlas__NodesFunctionSpace__mesh(NodesFunctionSpace* This)
{
  ASSERT(This);
  return &This->mesh();
}

Nodes* atlas__NodesFunctionSpace__nodes(NodesFunctionSpace* This)
{
  ASSERT(This);
  return &This->nodes();
}

Field* atlas__NodesFunctionSpace__create_field (const NodesFunctionSpace* This, const char* name, int kind )
{
  ASSERT(This);
  return This->createField(std::string(name),DataType::kind_to_datatype(kind));
}

Field* atlas__NodesFunctionSpace__create_field_vars (const NodesFunctionSpace* This, const char* name, int variables[], int variables_size, int fortran_ordering, int kind)
{
  ASSERT(This);
  ASSERT(variables_size);
  std::vector<size_t> variables_(variables_size);
  if( fortran_ordering )
    std::reverse_copy( variables, variables+variables_size,variables_.begin() );
  else
    variables_.assign(variables,variables+variables_size);
  return This->createField(std::string(name),variables_,DataType::kind_to_datatype(kind));
}

Field* atlas__NodesFunctionSpace__create_field_template (const NodesFunctionSpace* This, const char* name, const Field* field_template )
{
  ASSERT(This);
  return This->createField(std::string(name),*field_template);
}

Field* atlas__NodesFunctionSpace__create_global_field (const NodesFunctionSpace* This, const char* name, int kind )
{
  ASSERT(This);
  return This->createGlobalField(std::string(name),DataType::kind_to_datatype(kind));
}

Field* atlas__NodesFunctionSpace__create_global_field_vars (const NodesFunctionSpace* This, const char* name, int variables[], int variables_size, int fortran_ordering, int kind)
{
  ASSERT(This);
  ASSERT(variables_size);
  std::vector<size_t> variables_(variables_size);
  if( fortran_ordering )
    std::reverse_copy( variables, variables+variables_size, variables_.begin() );
  else
    variables_.assign(variables,variables+variables_size);
  return This->createGlobalField(std::string(name),variables_,DataType::kind_to_datatype(kind));
}

Field* atlas__NodesFunctionSpace__create_global_field_template (const NodesFunctionSpace* This, const char* name, const Field* field_template )
{
  ASSERT(This);
  return This->createGlobalField(std::string(name),*field_template);
}


void atlas__NodesFunctionSpace__halo_exchange_fieldset(const NodesFunctionSpace* This, FieldSet* fieldset)
{
  ASSERT(This);
  ASSERT(fieldset);
  ATLAS_ERROR_HANDLING( This->haloExchange(*fieldset); );
}

void atlas__NodesFunctionSpace__halo_exchange_field(const NodesFunctionSpace* This, Field* field)
{
  ASSERT(This);
  ASSERT(field);
  ATLAS_ERROR_HANDLING( This->haloExchange(*field); );
}

void atlas__NodesFunctionSpace__gather_fieldset(const NodesFunctionSpace* This, const FieldSet* local, FieldSet* global)
{
  ASSERT(This);
  ASSERT(local);
  ASSERT(global);
  ATLAS_ERROR_HANDLING( This->gather(*local,*global); );
}

void atlas__NodesFunctionSpace__gather_field(const NodesFunctionSpace* This, const Field* local, Field* global)
{
  ASSERT(This);
  ASSERT(local);
  ASSERT(global);
  ATLAS_ERROR_HANDLING( This->gather(*local,*global); );
}

void atlas__NodesFunctionSpace__scatter_fieldset(const NodesFunctionSpace* This, const FieldSet* global, FieldSet* local)
{
  ASSERT(This);
  ASSERT(local);
  ASSERT(global);
  ATLAS_ERROR_HANDLING( This->scatter(*global,*local); );
}

void atlas__NodesFunctionSpace__scatter_field(const NodesFunctionSpace* This, const Field* global, Field* local)
{
  ASSERT(This);
  ASSERT(global);
  ASSERT(local);
  ATLAS_ERROR_HANDLING( This->scatter(*global,*local); );
}

void atlas__NodesFunctionSpace__checksum_fieldset(const NodesFunctionSpace* This, const FieldSet* fieldset, char* &checksum, int &size, int &allocated)
{
  ASSERT(This);
  ASSERT(fieldset);
  ATLAS_ERROR_HANDLING(
    std::string checksum_str (This->checksum(*fieldset));
    size = checksum_str.size();
    checksum = new char[size+1]; allocated = true;
    strcpy(checksum,checksum_str.c_str());
  );
}

void atlas__NodesFunctionSpace__checksum_field(const NodesFunctionSpace* This, const Field* field, char* &checksum, int &size, int &allocated)
{
  ASSERT(This);
  ASSERT(field);
  ATLAS_ERROR_HANDLING(
    std::string checksum_str (This->checksum(*field));
    size = checksum_str.size();
    checksum = new char[size+1]; allocated = true;
    strcpy(checksum,checksum_str.c_str());
  );
}

void atlas__NodesFunctionSpace__sum_double(const NodesFunctionSpace* This, const Field* field, double &sum, int &N)
{
  ASSERT(This);
  ASSERT(field);
  size_t size_t_N;
  ATLAS_ERROR_HANDLING( This->sum(*field,sum,size_t_N) );
  N = size_t_N;
}

void atlas__NodesFunctionSpace__sum_arr_double(const NodesFunctionSpace* This, const Field* field, double* &sum, int &size, int &N)
{
  ASSERT(This);
  ASSERT(field);
  size_t size_t_N;
  ATLAS_ERROR_HANDLING(
        std::vector<double> sumvec;
        This->orderIndependentSum(*field,sumvec,size_t_N);
        size = sumvec.size();
        sum = new double[size];
        for( size_t j=0; j<size; ++j ) sum[j] = sumvec[j];
  );
  N = size_t_N;
}

void atlas__NodesFunctionSpace__oisum_double(const NodesFunctionSpace* This, const Field* field, double &sum, int &N)
{
  ASSERT(This);
  ASSERT(field);
  size_t size_t_N;
  ATLAS_ERROR_HANDLING( This->sum(*field,sum,size_t_N) );
  N = size_t_N;
}

void atlas__NodesFunctionSpace__oisum_arr_double(const NodesFunctionSpace* This, const Field* field, double* &sum, int &size, int &N)
{
  ASSERT(This);
  ASSERT(field);
  size_t size_t_N;
  ATLAS_ERROR_HANDLING(
        std::vector<double> sumvec;
        This->orderIndependentSum(*field,sumvec,size_t_N);
        size = sumvec.size();
        sum = new double[size];
        for( size_t j=0; j<size; ++j ) sum[j] = sumvec[j];
  );
  N = size_t_N;
}

void atlas__NodesFunctionSpace__min_double(const NodesFunctionSpace* This, const Field* field, double &minimum)
{
  ASSERT(This);
  ASSERT(field);
  ATLAS_ERROR_HANDLING( This->minimum(*field,minimum) );
}

void atlas__NodesFunctionSpace__max_double(const NodesFunctionSpace* This, const Field* field, double &maximum)
{
  ASSERT(This);
  ASSERT(field);
  ATLAS_ERROR_HANDLING( This->maximum(*field,maximum) );
}

void atlas__NodesFunctionSpace__min_arr_double(const NodesFunctionSpace* This, const Field* field, double* &minimum, int &size)
{
  ASSERT(This);
  ASSERT(field);
  ATLAS_ERROR_HANDLING(
        std::vector<double> minvec;
        This->minimum(*field,minvec);
        size = minvec.size();
        minimum = new double[size];
        for( size_t j=0; j<size; ++j ) minimum[j] = minvec[j];
  );
}

void atlas__NodesFunctionSpace__max_arr_double(const NodesFunctionSpace* This, const Field* field, double* &maximum, int &size)
{
  ASSERT(This);
  ASSERT(field);
  ATLAS_ERROR_HANDLING(
        std::vector<double> maxvec;
        This->maximum(*field,maxvec);
        size = maxvec.size();
        maximum = new double[size];
        for( size_t j=0; j<size; ++j ) maximum[j] = maxvec[j];
  );
}

void atlas__NodesFunctionSpace__minloc_double(const NodesFunctionSpace* This, const Field* field, double &minimum, long &glb_idx)
{
  ASSERT(This);
  ASSERT(field);
  gidx_t gidx;
  ATLAS_ERROR_HANDLING( This->minimumAndLocation(*field,minimum,gidx) );
  glb_idx = gidx;
}

void atlas__NodesFunctionSpace__maxloc_double(const NodesFunctionSpace* This, const Field* field, double &maximum, long &glb_idx)
{
  ASSERT(This);
  ASSERT(field);
  gidx_t gidx;
  ATLAS_ERROR_HANDLING( This->maximumAndLocation(*field,maximum,gidx) );
  glb_idx = gidx;
}

void atlas__NodesFunctionSpace__minloc_arr_double(const NodesFunctionSpace* This, const Field* field, double* &minimum, long* &glb_idx, int &size)
{
  ASSERT(This);
  ASSERT(field);
  ATLAS_ERROR_HANDLING(
        std::vector<double> minvec;
        std::vector<gidx_t> gidxvec;
        This->minimumAndLocation(*field,minvec,gidxvec);
        size = minvec.size();
        minimum = new double[size];
        glb_idx = new long[size];
        for( size_t j=0; j<size; ++j ) {
          minimum[j] = minvec[j];
          glb_idx[j] = gidxvec[j];
        }
  );
}

void atlas__NodesFunctionSpace__maxloc_arr_double(const NodesFunctionSpace* This, const Field* field, double* &maximum, long* &glb_idx, int &size)
{
  ASSERT(This);
  ASSERT(field);
  ATLAS_ERROR_HANDLING(
        std::vector<double> maxvec;
        std::vector<gidx_t> gidxvec;
        This->maximumAndLocation(*field,maxvec,gidxvec);
        size = maxvec.size();
        maximum = new double[size];
        glb_idx = new long[size];
        for( size_t j=0; j<size; ++j ) {
          maximum[j] = maxvec[j];
          glb_idx[j] = gidxvec[j];
        }
  );
}

void atlas__NodesFunctionSpace__mean_double(const NodesFunctionSpace* This, const Field* field, double &mean, int &N)
{
  ASSERT(This);
  ASSERT(field);
  size_t size_t_N;
  ATLAS_ERROR_HANDLING( This->mean(*field,mean,size_t_N) );
  N = size_t_N;
}

void atlas__NodesFunctionSpace__mean_arr_double(const NodesFunctionSpace* This, const Field* field, double* &mean, int &size, int &N)
{
  ASSERT(This);
  ASSERT(field);
  size_t size_t_N;
  ATLAS_ERROR_HANDLING(
        std::vector<double> meanvec;
        This->mean(*field,meanvec,size_t_N);
        size = meanvec.size();
        mean = new double[size];
        for( size_t j=0; j<size; ++j ) mean[j] = meanvec[j];
  );
  N = size_t_N;
}

void atlas__NodesFunctionSpace__mean_and_stddev_double(const NodesFunctionSpace* This, const Field* field, double &mean, double &stddev, int &N)
{
  ASSERT(This);
  ASSERT(field);
  size_t size_t_N;
  ATLAS_ERROR_HANDLING( This->meanAndStandardDeviation(*field,mean,stddev,size_t_N) );
  N = size_t_N;
}

void atlas__NodesFunctionSpace__mean_and_stddev_arr_double(const NodesFunctionSpace* This, const Field* field, double* &mean, double* &stddev, int &size, int &N)
{
  ASSERT(This);
  ASSERT(field);
  size_t size_t_N;
  ATLAS_ERROR_HANDLING(
        std::vector<double> meanvec;
        std::vector<double> stddevvec;
        This->meanAndStandardDeviation(*field,meanvec,stddevvec,size_t_N);
        size = meanvec.size();
        mean = new double[size];
        stddev = new double[size];
        for( size_t j=0; j<size; ++j ) {
          mean[j] = meanvec[j];
          stddev[j] = stddevvec[j];
        }
  );
  N = size_t_N;
}





NodesColumnFunctionSpace* atlas__NodesColumnFunctionSpace__new (const char* name, Mesh* mesh, int nb_levels, int halo)
{
  ASSERT(mesh);
  return new NodesColumnFunctionSpace(std::string(name),*mesh,nb_levels,Halo(halo));
}

void atlas__NodesColumnFunctionSpace__delete (NodesColumnFunctionSpace* This)
{
  ASSERT(This);
  delete(This);
}

int atlas__NodesColumnFunctionSpace__nb_levels(const NodesColumnFunctionSpace* This)
{
  ASSERT(This);
  return This->nb_levels();
}

Field* atlas__NodesColumnFunctionSpace__create_field (const NodesColumnFunctionSpace* This, const char* name, int kind )
{
  ASSERT(This);
  return This->createField(std::string(name),DataType::kind_to_datatype(kind));
}

Field* atlas__NodesColumnFunctionSpace__create_field_vars (const NodesColumnFunctionSpace* This, const char* name, int variables[], int variables_size, int fortran_ordering, int kind)
{
  ASSERT(This);
  ASSERT(variables_size);
  std::vector<size_t> variables_(variables_size);
  if( fortran_ordering )
    std::reverse_copy( variables, variables+variables_size, variables_.begin() );
  else
    variables_.assign(variables,variables+variables_size);
  return This->createField(std::string(name),variables_,DataType::kind_to_datatype(kind));
}

Field* atlas__NodesColumnFunctionSpace__create_field_template (const NodesColumnFunctionSpace* This, const char* name, const Field* field_template )
{
  ASSERT(This);
  return This->createField(std::string(name),*field_template);
}

Field* atlas__NodesColumnFunctionSpace__create_global_field (const NodesColumnFunctionSpace* This, const char* name, int kind )
{
  ASSERT(This);
  return This->createGlobalField(std::string(name),DataType::kind_to_datatype(kind));
}

Field* atlas__NodesColumnFunctionSpace__create_global_field_vars (const NodesColumnFunctionSpace* This, const char* name, int variables[], int variables_size, int fortran_ordering, int kind)
{
  ASSERT(This);
  ASSERT(variables_size);
  std::vector<size_t> variables_(variables_size);
  if( fortran_ordering )
    std::reverse_copy( variables, variables+variables_size, variables_.begin() );
  else
    variables_.assign(variables,variables+variables_size);
  return This->createGlobalField(std::string(name),variables_,DataType::kind_to_datatype(kind));
}

Field* atlas__NodesColumnFunctionSpace__create_global_field_template (const NodesColumnFunctionSpace* This, const char* name, const Field* field_template )
{
  ASSERT(This);
  return This->createGlobalField(std::string(name),*field_template);
}


void atlas__NodesColumnFunctionSpace__checksum_fieldset(const NodesColumnFunctionSpace* This, const FieldSet* fieldset, char* &checksum, int &size, int &allocated)
{
  ASSERT(This);
  ASSERT(fieldset);
  ATLAS_ERROR_HANDLING(
    std::string checksum_str (This->checksum(*fieldset));
    size = checksum_str.size();
    checksum = new char[size+1]; allocated = true;
    strcpy(checksum,checksum_str.c_str());
  );
}

void atlas__NodesColumnFunctionSpace__checksum_field(const NodesColumnFunctionSpace* This, const Field* field, char* &checksum, int &size, int &allocated)
{
  ASSERT(This);
  ASSERT(field);
  ATLAS_ERROR_HANDLING(
    std::string checksum_str (This->checksum(*field));
    size = checksum_str.size();
    checksum = new char[size+1]; allocated = true;
    strcpy(checksum,checksum_str.c_str());
  );
}





}


} // namespace functionspace
} // namespace atlas

