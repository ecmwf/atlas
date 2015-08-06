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
#include "atlas/util/Bitflags.h"
#include "atlas/functionspace/NodesFunctionSpace.h"
#include "atlas/field/FieldT.h"
#include "atlas/actions/BuildParallelFields.h"
#include "atlas/actions/BuildHalo.h"
#include "atlas/actions/BuildPeriodicBoundaries.h"

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

    ArrayView<int,1> flags ( mesh_.nodes().field("flags") );
    std::vector<int> mask(mesh_.nodes().size());
    for( size_t j=0; j<mask.size(); ++j ) {
      mask[j] = util::Topology::check(flags(j),util::Topology::GHOST) ? 1 : 0;

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

    ArrayView<int,1> flags ( mesh_.nodes().field("flags") );
    std::vector<int> mask(mesh_.nodes().size());
    for( size_t j=0; j<mask.size(); ++j ) {
      mask[j] = util::Topology::check(flags(j),util::Topology::GHOST) ? 1 : 0;
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
      if     ( field.datatype() == DataType::datatype<int>() ) {
        ArrayView<int,2> view(field.data<int>(),strides.data(),shape.data());
        halo_exchange.execute( view );
      }
      else if( field.datatype() == DataType::datatype<long>() ) {
        ArrayView<long,2> view(field.data<long>(),strides.data(),shape.data());
        halo_exchange.execute( view );
      }
      else if( field.datatype() == DataType::datatype<float>() ) {
        ArrayView<float,2> view(field.data<float>(),strides.data(),shape.data());
        halo_exchange.execute( view );
      }
      else if( field.datatype() == DataType::datatype<double>() ) {
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

    if     ( loc.datatype() == DataType::datatype<int>() ) {
      mpl::Field<int const> loc_field(loc.data<int>(),loc.stride(0));
      mpl::Field<int      > glb_field(glb.data<int>(),glb.stride(0));
      gather_scatter.gather( &loc_field, &glb_field, 1 );
    }
    else if( loc.datatype() == DataType::datatype<long>() ) {
      mpl::Field<long const> loc_field(loc.data<long>(),loc.stride(0));
      mpl::Field<long      > glb_field(glb.data<long>(),glb.stride(0));
      gather_scatter.gather( &loc_field, &glb_field, 1 );
    }
    else if( loc.datatype() == DataType::datatype<float>() ) {
      mpl::Field<float const> loc_field(loc.data<float>(),loc.stride(0));
      mpl::Field<float      > glb_field(glb.data<float>(),glb.stride(0));
      gather_scatter.gather( &loc_field, &glb_field, 1 );
    }
    else if( loc.datatype() == DataType::datatype<double>() ) {
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

    if     ( loc.datatype() == DataType::datatype<int>() ) {
      mpl::Field<int const> glb_field(glb.data<int>(),glb.stride(0));
      mpl::Field<int      > loc_field(loc.data<int>(),loc.stride(0));
      gather_scatter.scatter( &glb_field, &loc_field, 1 );
    }
    else if( loc.datatype() == DataType::datatype<long>() ) {
      mpl::Field<long const> glb_field(glb.data<long>(),glb.stride(0));
      mpl::Field<long      > loc_field(loc.data<long>(),loc.stride(0));
      gather_scatter.scatter( &glb_field, &loc_field, 1 );
    }
    else if( loc.datatype() == DataType::datatype<float>() ) {
      mpl::Field<float const> glb_field(glb.data<float>(),glb.stride(0));
      mpl::Field<float      > loc_field(loc.data<float>(),loc.stride(0));
      gather_scatter.scatter( &glb_field, &loc_field, 1 );
    }
    else if( loc.datatype() == DataType::datatype<double>() ) {
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
    if     ( field.datatype() == DataType::datatype<int>() )
      md5 << checksum.execute( field.data<int>(), field.stride(0) );
    else if( field.datatype() == DataType::datatype<long>() )
      md5 << checksum.execute( field.data<long>(), field.stride(0) );
    else if( field.datatype() == DataType::datatype<float>() )
      md5 << checksum.execute( field.data<float>(), field.stride(0) );
    else if( field.datatype() == DataType::datatype<double>() )
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

template<> void NodesFunctionSpace::sum( const Field& field, double& result, size_t& N ) const
{
  size_t root = 0;
  Field::Ptr global( createGlobalField(field) );
  gather(field,*global);
  result = std::accumulate(global->data<double>(),global->data<double>()+global->size(),0.);
  eckit::mpi::broadcast(result,root);
  N = nb_nodes_global_broadcasted_;
}

template<> void NodesFunctionSpace::sum( const Field& field, std::vector<double>& result, size_t& N ) const
{
  size_t nvar = field.stride(0);
  result.resize(nvar,0);
  size_t root = 0;
  Field::Ptr global( createGlobalField(field) );
  gather(field,*global);
  if( eckit::mpi::rank() == 0 ) {
    const ArrayView<double,2> glb( global->data<double>(), make_shape(global->shape(0),global->stride(0)) );
    for( size_t n=0; n<nb_nodes_global(); ++n ) {
      for( size_t j=0; j<nvar; ++j ) {
        result[j] += glb(n,j);
      }
    }
  }
  eckit::mpi::broadcast(result,root);
  N = nb_nodes_global_broadcasted_;
}

namespace { inline double sqr(const double& val) { return val*val; } }

template<> void NodesFunctionSpace::maximum( const Field& field, double& result, gidx_t& glb_idx ) const
{
  throw eckit::Exception("Untested and no glb_idx return",Here());
  double local_maximum = *std::max_element(field.data<double>(),field.data<double>()+field.size());
  eckit::mpi::all_reduce(local_maximum,result,eckit::mpi::max());
}

template<> void NodesFunctionSpace::maximum( const Field& field, std::vector<double>& result, gidx_t& glb_idx ) const
{
  throw eckit::Exception("Untested and no glb_idx return",Here());
  size_t nvar = field.stride(0);
  std::vector<double> local_maximum(nvar,std::numeric_limits<double>::min());
  ArrayView<double,2> values( field.data<double>(), make_shape(field.shape(0),nvar) );
  for( size_t n=0; n<values.shape(0); ++n ) {
    for( size_t j=0; j<nvar; ++j ) {
      local_maximum[j] = std::max(values(n,j),local_maximum[j]);
    }
  }
  eckit::mpi::all_reduce(local_maximum,result,eckit::mpi::max());
}

template<> void NodesFunctionSpace::minimum( const Field& field, double& result, gidx_t& glb_idx ) const
{
  throw eckit::Exception("Untested and no glb_idx return",Here());
  double local_minimum = *std::min_element(field.data<double>(),field.data<double>()+field.size());
  eckit::mpi::all_reduce(local_minimum,result,eckit::mpi::min());
}

template<> void NodesFunctionSpace::minimum( const Field& field, std::vector<double>& result, std::vector<gidx_t>& glb_idx ) const
{
  throw eckit::Exception("Untested and no glb_idx return",Here());
  size_t nvar = field.stride(0);
  std::vector<double> local_minimum(nvar,std::numeric_limits<double>::max());
  ArrayView<double,2> values( field.data<double>(), make_shape(field.shape(0),nvar) );
  for( size_t n=0; n<values.shape(0); ++n ) {
    for( size_t j=0; j<nvar; ++j ) {
      local_minimum[j] = std::min(values(n,j),local_minimum[j]);
    }
  }
  eckit::mpi::all_reduce(local_minimum,result,eckit::mpi::min());
}

template<> void NodesFunctionSpace::mean( const Field& field, double& result, size_t& N ) const
{
  throw eckit::Exception("Untested",Here());
  sum(field,result,N);
  result /= static_cast<double>(N);
}

template<> void NodesFunctionSpace::mean( const Field& field, std::vector<double>& result, size_t& N ) const
{
  throw eckit::Exception("Untested",Here());
  sum(field,result,N);
  for( size_t j=0; j<result.size(); ++j ) {
    result[j] /= static_cast<double>(N);
  }
}

template<> void NodesFunctionSpace::mean_and_standard_deviation( const Field& field, double& mu, double& sigma, size_t& N ) const
{
  throw eckit::Exception("Untested",Here());
  mean(field,mu,N);
  Field::Ptr squared_diff_field( createField(field) );
  ArrayView<double,1> squared_diff( *squared_diff_field );
  ArrayView<double,1> values( field );
  for( size_t n=0; n<field.size(); ++n )
    squared_diff(n) = sqr( values(n) - mu );

  mean(*squared_diff_field,sigma,N);
  sigma = std::sqrt(sigma);
}

template<> void NodesFunctionSpace::mean_and_standard_deviation( const Field& field, std::vector<double>& mu, std::vector<double>& sigma, size_t& N ) const
{
  throw eckit::Exception("Untested",Here());
  mean(field,mu,N);
  Field::Ptr squared_diff_field( createField(field) );
  ArrayView<double,2> squared_diff( squared_diff_field->data<double>(), make_shape(squared_diff_field->shape(0),squared_diff_field->stride(0)) );
  ArrayView<double,2> values( field.data<double>(), make_shape(field.shape(0),field.stride(0)) );
  for( size_t n=0; n<values.shape(0); ++n ) {
    for( size_t j=0; j<values.shape(1); ++j ) {
      squared_diff(n,j) = sqr( values(n,j) - mu[j] );
    }
  }
  mean(*squared_diff_field,sigma,N);
  for( size_t j=0; j<sigma.size(); ++j ) {
    sigma[j] = std::sqrt(sigma[j]);
  }
}


// ----------------------------------------------------------------------

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
}


} // namespace functionspace
} // namespace atlas

