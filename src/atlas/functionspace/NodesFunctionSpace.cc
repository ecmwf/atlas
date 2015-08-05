/*
 * (C) Copyright 1996-2014 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "atlas/atlas_config.h"
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
    halo_(halo.size()),
    nb_nodes_(0),
    nb_nodes_global_(0)
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
      if( mask[j] == 1 && util::Topology::check(flags(j),util::Topology::BC) ) {
        mask[j] = 0;
      }
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

template<> Field* NodesFunctionSpace::createField<int>(const std::string& name) const {
  return new field::FieldT<int>(name, make_shape(nb_nodes()));
}

template<> Field* NodesFunctionSpace::createField<long>(const std::string& name) const {
  return new field::FieldT<long>(name, make_shape(nb_nodes()) );
}

template<> Field* NodesFunctionSpace::createField<float>(const std::string& name) const {
  return new field::FieldT<float>(name, make_shape(nb_nodes()) );
}

template<> Field* NodesFunctionSpace::createField<double>(const std::string& name) const {
  return new field::FieldT<double>(name, make_shape(nb_nodes()) );
}

template<> Field* NodesFunctionSpace::createField<int>(const std::string& name, const std::vector<size_t>& variables) const {
  std::vector<size_t> shape(1,nb_nodes());
  for( size_t i=0; i<variables.size(); ++i ) shape.push_back(variables[i]);
  return new field::FieldT<int>(name, shape);
}

template<> Field* NodesFunctionSpace::createField<long>(const std::string& name, const std::vector<size_t>& variables) const {
  std::vector<size_t> shape(1,nb_nodes());
  for( size_t i=0; i<variables.size(); ++i ) shape.push_back(variables[i]);
  return new field::FieldT<long>(name, shape);
}

template<> Field* NodesFunctionSpace::createField<float>(const std::string& name, const std::vector<size_t>& variables) const {
  std::vector<size_t> shape(1,nb_nodes());
  for( size_t i=0; i<variables.size(); ++i ) shape.push_back(variables[i]);
  return new field::FieldT<float>(name, shape);
}

template<> Field* NodesFunctionSpace::createField<double>(const std::string& name, const std::vector<size_t>& variables) const {
  std::vector<size_t> shape(1,nb_nodes());
  for( size_t i=0; i<variables.size(); ++i ) shape.push_back(variables[i]);
  return new field::FieldT<double>(name, shape);
}

Field* NodesFunctionSpace::createField(const std::string& name, const Field& other) const {
  ArrayShape shape = other.shape();
  shape[0] = nb_nodes();
  if( other.datatype() == DataType::datatype<int>() )
    return new field::FieldT<int>(name,shape);
  else if( other.datatype() == DataType::datatype<long>() )
    return new field::FieldT<long>(name,shape);
  else if( other.datatype() == DataType::datatype<float>() )
    return new field::FieldT<float>(name,shape);
  else if( other.datatype() == DataType::datatype<double>() )
    return new field::FieldT<double>(name,shape);
  else
    eckit::Exception("Should not be here",Here());
  return 0;
}

template<> Field* NodesFunctionSpace::createGlobalField<int>(const std::string& name) const {
  return new field::FieldT<int>(name, make_shape(nb_nodes_global()));
}

template<> Field* NodesFunctionSpace::createGlobalField<long>(const std::string& name) const {
  return new field::FieldT<long>(name, make_shape(nb_nodes_global()) );
}

template<> Field* NodesFunctionSpace::createGlobalField<float>(const std::string& name) const {
  return new field::FieldT<float>(name, make_shape(nb_nodes_global()) );
}

template<> Field* NodesFunctionSpace::createGlobalField<double>(const std::string& name) const {
  return new field::FieldT<double>(name, make_shape(nb_nodes_global()) );
}

template<> Field* NodesFunctionSpace::createGlobalField<int>(const std::string& name, const std::vector<size_t>& variables) const {
  std::vector<size_t> shape(1,nb_nodes_global());
  for( size_t i=0; i<variables.size(); ++i ) shape.push_back(variables[i]);
  return new field::FieldT<int>(name, shape);
}

template<> Field* NodesFunctionSpace::createGlobalField<long>(const std::string& name, const std::vector<size_t>& variables) const {
  std::vector<size_t> shape(1,nb_nodes_global());
  for( size_t i=0; i<variables.size(); ++i ) shape.push_back(variables[i]);
  return new field::FieldT<long>(name, shape);
}

template<> Field* NodesFunctionSpace::createGlobalField<float>(const std::string& name, const std::vector<size_t>& variables) const {
  std::vector<size_t> shape(1,nb_nodes_global());
  for( size_t i=0; i<variables.size(); ++i ) shape.push_back(variables[i]);
  return new field::FieldT<float>(name, shape);
}

template<> Field* NodesFunctionSpace::createGlobalField<double>(const std::string& name, const std::vector<size_t>& variables) const {
  std::vector<size_t> shape(1,nb_nodes_global());
  for( size_t i=0; i<variables.size(); ++i ) shape.push_back(variables[i]);
  return new field::FieldT<double>(name, shape);
}

Field* NodesFunctionSpace::createGlobalField(const std::string& name, const Field& other) const
{
  ArrayShape shape = other.shape();
  shape[0] = nb_nodes_global();
  if( other.datatype() == DataType::datatype<int>() )
    return new field::FieldT<int>(name,shape);
  else if( other.datatype() == DataType::datatype<long>() )
    return new field::FieldT<long>(name,shape);
  else if( other.datatype() == DataType::datatype<float>() )
    return new field::FieldT<float>(name,shape);
  else if( other.datatype() == DataType::datatype<double>() )
    return new field::FieldT<double>(name,shape);
  else
    eckit::Exception("Should not be here",Here());
  return 0;
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

template<> Field* NodesColumnFunctionSpace::createField<int>(const std::string& name) const {
  return new field::FieldT<int>(name, make_shape(nb_nodes(),nb_levels_) );
}

template<> Field* NodesColumnFunctionSpace::createField<long>(const std::string& name) const {
  return new field::FieldT<long>(name, make_shape(nb_nodes(),nb_levels_) );
}

template<> Field* NodesColumnFunctionSpace::createField<float>(const std::string& name) const {
  return new field::FieldT<float>(name, make_shape(nb_nodes(),nb_levels_) );
}

template<> Field* NodesColumnFunctionSpace::createField<double>(const std::string& name) const {
  return new field::FieldT<double>(name, make_shape(nb_nodes(),nb_levels_) );
}

template<> Field* NodesColumnFunctionSpace::createField<int>(const std::string& name, size_t var1) const {
  return new field::FieldT<int>(name, make_shape(nb_nodes(),nb_levels_,var1) );
}

template<> Field* NodesColumnFunctionSpace::createField<long>(const std::string& name, size_t var1) const {
  return new field::FieldT<long>(name, make_shape(nb_nodes(),nb_levels_,var1) );
}

template<> Field* NodesColumnFunctionSpace::createField<float>(const std::string& name, size_t var1) const {
  return new field::FieldT<float>(name, make_shape(nb_nodes(),nb_levels_,var1) );
}

template<> Field* NodesColumnFunctionSpace::createField<double>(const std::string& name, size_t var1) const {
  return new field::FieldT<double>(name, make_shape(nb_nodes(),nb_levels_,var1) );
}

template<> Field* NodesColumnFunctionSpace::createField<int>(const std::string& name, size_t var1, size_t var2) const {
  return new field::FieldT<int>(name, make_shape(nb_nodes(),nb_levels_,var1,var2) );
}

template<> Field* NodesColumnFunctionSpace::createField<long>(const std::string& name, size_t var1, size_t var2) const {
  return new field::FieldT<long>(name, make_shape(nb_nodes(),nb_levels_,var1,var2) );
}

template<> Field* NodesColumnFunctionSpace::createField<float>(const std::string& name, size_t var1, size_t var2) const {
  return new field::FieldT<float>(name, make_shape(nb_nodes(),nb_levels_,var1,var2) );
}

template<> Field* NodesColumnFunctionSpace::createField<double>(const std::string& name, size_t var1, size_t var2) const {
  return new field::FieldT<double>(name, make_shape(nb_nodes(),nb_levels_,var1,var2) );
}

template<> Field* NodesColumnFunctionSpace::createField<int>(const std::string& name, const std::vector<size_t>& variables) const {
  std::vector<size_t> shape(1,nb_nodes());   shape.push_back(nb_levels_);
  for( size_t i=0; i<variables.size(); ++i ) shape.push_back(variables[i]);
  return new field::FieldT<int>(name, shape);
}

template<> Field* NodesColumnFunctionSpace::createField<long>(const std::string& name, const std::vector<size_t>& variables) const {
  std::vector<size_t> shape(1,nb_nodes());   shape.push_back(nb_levels_);
  for( size_t i=0; i<variables.size(); ++i ) shape.push_back(variables[i]);
  return new field::FieldT<long>(name, shape);
}

template<> Field* NodesColumnFunctionSpace::createField<float>(const std::string& name, const std::vector<size_t>& variables) const {
  std::vector<size_t> shape(1,nb_nodes());   shape.push_back(nb_levels_);
  for( size_t i=0; i<variables.size(); ++i ) shape.push_back(variables[i]);
  return new field::FieldT<float>(name, shape);
}

template<> Field* NodesColumnFunctionSpace::createField<double>(const std::string& name, const std::vector<size_t>& variables) const {
  std::vector<size_t> shape(1,nb_nodes());   shape.push_back(nb_levels_);
  for( size_t i=0; i<variables.size(); ++i ) shape.push_back(variables[i]);
  return new field::FieldT<double>(name, shape);
}

template<> Field* NodesColumnFunctionSpace::createGlobalField<int>(const std::string& name) const {
  return new field::FieldT<int>(name, make_shape(nb_nodes_global(),nb_levels_));
}

template<> Field* NodesColumnFunctionSpace::createGlobalField<long>(const std::string& name) const {
  return new field::FieldT<long>(name, make_shape(nb_nodes_global(),nb_levels_) );
}

template<> Field* NodesColumnFunctionSpace::createGlobalField<float>(const std::string& name) const {
  return new field::FieldT<float>(name, make_shape(nb_nodes_global(),nb_levels_) );
}

template<> Field* NodesColumnFunctionSpace::createGlobalField<double>(const std::string& name) const {
  return new field::FieldT<double>(name, make_shape(nb_nodes_global(),nb_levels_) );
}

template<> Field* NodesColumnFunctionSpace::createGlobalField<int>(const std::string& name, size_t var1) const {
  return new field::FieldT<int>(name, make_shape(nb_nodes_global(),nb_levels_,var1) );
}

template<> Field* NodesColumnFunctionSpace::createGlobalField<long>(const std::string& name, size_t var1) const {
  return new field::FieldT<long>(name, make_shape(nb_nodes_global(),nb_levels_,var1) );
}

template<> Field* NodesColumnFunctionSpace::createGlobalField<float>(const std::string& name, size_t var1) const {
  return new field::FieldT<float>(name, make_shape(nb_nodes_global(),nb_levels_,var1) );
}

template<> Field* NodesColumnFunctionSpace::createGlobalField<double>(const std::string& name, size_t var1) const {
  return new field::FieldT<double>(name, make_shape(nb_nodes_global(),nb_levels_,var1) );
}

template<> Field* NodesColumnFunctionSpace::createGlobalField<int>(const std::string& name, size_t var1, size_t var2) const {
  return new field::FieldT<int>(name, make_shape(nb_nodes_global(),nb_levels_,var1,var2) );
}

template<> Field* NodesColumnFunctionSpace::createGlobalField<long>(const std::string& name, size_t var1, size_t var2) const {
  return new field::FieldT<long>(name, make_shape(nb_nodes_global(),nb_levels_,var1,var2) );
}

template<> Field* NodesColumnFunctionSpace::createGlobalField<float>(const std::string& name, size_t var1, size_t var2) const {
  return new field::FieldT<float>(name, make_shape(nb_nodes_global(),nb_levels_,var1,var2) );
}

template<> Field* NodesColumnFunctionSpace::createGlobalField<double>(const std::string& name, size_t var1, size_t var2) const {
  return new field::FieldT<double>(name, make_shape(nb_nodes_global(),nb_levels_,var1,var2) );
}

template<> Field* NodesColumnFunctionSpace::createGlobalField<int>(const std::string& name, const std::vector<size_t>& variables) const {
  std::vector<size_t> shape(1,nb_nodes());   shape.push_back(nb_levels_);
  for( size_t i=0; i<variables.size(); ++i ) shape.push_back(variables[i]);
  return new field::FieldT<int>(name, shape);
}

template<> Field* NodesColumnFunctionSpace::createGlobalField<long>(const std::string& name, const std::vector<size_t>& variables) const {
  std::vector<size_t> shape(1,nb_nodes_global());   shape.push_back(nb_levels_);
  for( size_t i=0; i<variables.size(); ++i ) shape.push_back(variables[i]);
  return new field::FieldT<long>(name, shape);
}

template<> Field* NodesColumnFunctionSpace::createGlobalField<float>(const std::string& name, const std::vector<size_t>& variables) const {
  std::vector<size_t> shape(1,nb_nodes());   shape.push_back(nb_levels_);
  for( size_t i=0; i<variables.size(); ++i ) shape.push_back(variables[i]);
  return new field::FieldT<float>(name, shape);
}

template<> Field* NodesColumnFunctionSpace::createGlobalField<double>(const std::string& name, const std::vector<size_t>& variables) const {
  std::vector<size_t> shape(1,nb_nodes());   shape.push_back(nb_levels_);
  for( size_t i=0; i<variables.size(); ++i ) shape.push_back(variables[i]);
  return new field::FieldT<double>(name, shape);
}

Field* NodesColumnFunctionSpace::createGlobalField(const std::string& name, const Field& other) const
{
  ArrayShape shape = other.shape();
  shape[0] = nb_nodes_global();
  if( other.datatype() == DataType::datatype<int>() )
    return new field::FieldT<int>(name,shape);
  else if( other.datatype() == DataType::datatype<long>() )
    return new field::FieldT<long>(name,shape);
  else if( other.datatype() == DataType::datatype<float>() )
    return new field::FieldT<float>(name,shape);
  else if( other.datatype() == DataType::datatype<double>() )
    return new field::FieldT<double>(name,shape);
  else
    eckit::Exception("Should not be here",Here());
  return 0;
}


// ----------------------------------------------------------------------

} // namespace functionspace
} // namespace atlas

