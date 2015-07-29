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
#include "atlas/mpl/HaloExchange.h"
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

NodesFunctionSpace::NodesFunctionSpace(const std::string& name, Mesh& mesh, size_t halo)
  : next::FunctionSpace(name),
    mesh_(mesh),
    halo_(halo),
    nb_nodes_(0)
{
  if( ! mesh_.halo_exchange().has(halo_name()) && halo_ > 0)
  {
    // Create new halo-exchange
    mpl::HaloExchange* halo_exchange = new mpl::HaloExchange( halo_name() );

    // Set it up.
    actions::build_nodes_parallel_fields( mesh_.function_space("nodes") );

    actions::build_periodic_boundaries(mesh_);

    actions::build_halo(mesh_,halo_);

    actions::renumber_nodes_glb_idx(mesh_.function_space("nodes"));

    Field& ridx = mesh_.function_space("nodes").field("remote_idx");
    Field& part = mesh_.function_space("nodes").field("partition");
    Field& gidx = mesh_.function_space("nodes").field("glb_idx");

    ArrayView<int,1> flags ( mesh_.function_space("nodes").field("flags") );
    std::vector<int> mask(mesh_.function_space("nodes").shape(0));
    for( size_t j=0; j<mask.size(); ++j )
    {
      mask[j] = util::Topology::check(flags(j),util::Topology::GHOST) ? 1 : 0;
    }

    std::stringstream ss;
    ss << "nb_nodes_including_halo["<<halo_<<"]";
    mesh.metadata().get(ss.str(),nb_nodes_);

    halo_exchange->setup(part.data<int>(),ridx.data<int>(),REMOTE_IDX_BASE,nb_nodes_);

    // Store it in the mesh
    mesh_.halo_exchange().add(halo_exchange);
  }
  if( !nb_nodes_ )
    nb_nodes_ = mesh_.function_space("nodes").metadata().get<size_t>("nb_owned");
}

NodesFunctionSpace::~NodesFunctionSpace() {}

size_t NodesFunctionSpace::nb_nodes() const
{
  return nb_nodes_;
}

std::string NodesFunctionSpace::halo_name() const
{
  std::stringstream ss; ss << "nodes_" << halo_;
  return ss.str();
}

template<> Field* NodesFunctionSpace::createField<int>(const std::string& name) {
  return new field::FieldT<int>(name, make_shape(nb_nodes()));
}

template<> Field* NodesFunctionSpace::createField<long>(const std::string& name) {
  return new field::FieldT<long>(name, make_shape(nb_nodes()) );
}

template<> Field* NodesFunctionSpace::createField<float>(const std::string& name) {
  return new field::FieldT<float>(name, make_shape(nb_nodes()) );
}

template<> Field* NodesFunctionSpace::createField<double>(const std::string& name) {
  return new field::FieldT<double>(name, make_shape(nb_nodes()) );
}

template<> Field* NodesFunctionSpace::createField<int>(const std::string& name, size_t var1) {
  return new field::FieldT<int>(name, make_shape(nb_nodes(),var1) );
}

template<> Field* NodesFunctionSpace::createField<long>(const std::string& name, size_t var1) {
  return new field::FieldT<long>(name, make_shape(nb_nodes(),var1) );
}

template<> Field* NodesFunctionSpace::createField<float>(const std::string& name, size_t var1) {
  return new field::FieldT<float>(name, make_shape(nb_nodes(),var1) );
}

template<> Field* NodesFunctionSpace::createField<double>(const std::string& name, size_t var1) {
  return new field::FieldT<double>(name, make_shape(nb_nodes(),var1) );
}

template<> Field* NodesFunctionSpace::createField<int>(const std::string& name, size_t var1, size_t var2) {
  return new field::FieldT<int>(name, make_shape(nb_nodes(),var1,var2) );
}

template<> Field* NodesFunctionSpace::createField<long>(const std::string& name, size_t var1, size_t var2) {
  return new field::FieldT<long>(name, make_shape(nb_nodes(),var1,var2) );
}

template<> Field* NodesFunctionSpace::createField<float>(const std::string& name, size_t var1, size_t var2) {
  return new field::FieldT<float>(name, make_shape(nb_nodes(),var1,var2) );
}

template<> Field* NodesFunctionSpace::createField<double>(const std::string& name, size_t var1, size_t var2) {
  return new field::FieldT<double>(name, make_shape(nb_nodes(),var1,var2) );
}


void NodesFunctionSpace::haloExchange( FieldSet& fieldset )
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
void NodesFunctionSpace::haloExchange( Field& field )
{
  if( halo_ ) {
    FieldSet fieldset;
    fieldset.add(field);
    haloExchange(fieldset);
  }
}



// ----------------------------------------------------------------------

NodesColumnFunctionSpace::NodesColumnFunctionSpace(const std::string& name, Mesh& mesh, size_t nb_levels, size_t halo)
  : NodesFunctionSpace(name,mesh,halo),
    nb_levels_(nb_levels)
{
}

NodesColumnFunctionSpace::~NodesColumnFunctionSpace() {}

size_t NodesColumnFunctionSpace::nb_levels() const
{
  return nb_levels_;
}

template<> Field* NodesColumnFunctionSpace::createField<int>(const std::string& name) {
  return new field::FieldT<int>(name, make_shape(nb_nodes(),nb_levels_) );
}

template<> Field* NodesColumnFunctionSpace::createField<long>(const std::string& name) {
  return new field::FieldT<long>(name, make_shape(nb_nodes(),nb_levels_) );
}

template<> Field* NodesColumnFunctionSpace::createField<float>(const std::string& name) {
  return new field::FieldT<float>(name, make_shape(nb_nodes(),nb_levels_) );
}

template<> Field* NodesColumnFunctionSpace::createField<double>(const std::string& name) {
  return new field::FieldT<double>(name, make_shape(nb_nodes(),nb_levels_) );
}

template<> Field* NodesColumnFunctionSpace::createField<int>(const std::string& name, size_t var1) {
  return new field::FieldT<int>(name, make_shape(nb_nodes(),nb_levels_,var1) );
}

template<> Field* NodesColumnFunctionSpace::createField<long>(const std::string& name, size_t var1) {
  return new field::FieldT<long>(name, make_shape(nb_nodes(),nb_levels_,var1) );
}

template<> Field* NodesColumnFunctionSpace::createField<float>(const std::string& name, size_t var1) {
  return new field::FieldT<float>(name, make_shape(nb_nodes(),nb_levels_,var1) );
}

template<> Field* NodesColumnFunctionSpace::createField<double>(const std::string& name, size_t var1) {
  return new field::FieldT<double>(name, make_shape(nb_nodes(),nb_levels_,var1) );
}

template<> Field* NodesColumnFunctionSpace::createField<int>(const std::string& name, size_t var1, size_t var2) {
  return new field::FieldT<int>(name, make_shape(nb_nodes(),nb_levels_,var1,var2) );
}

template<> Field* NodesColumnFunctionSpace::createField<long>(const std::string& name, size_t var1, size_t var2) {
  return new field::FieldT<long>(name, make_shape(nb_nodes(),nb_levels_,var1,var2) );
}

template<> Field* NodesColumnFunctionSpace::createField<float>(const std::string& name, size_t var1, size_t var2) {
  return new field::FieldT<float>(name, make_shape(nb_nodes(),nb_levels_,var1,var2) );
}

template<> Field* NodesColumnFunctionSpace::createField<double>(const std::string& name, size_t var1, size_t var2) {
  return new field::FieldT<double>(name, make_shape(nb_nodes(),nb_levels_,var1,var2) );
}

// ----------------------------------------------------------------------

} // namespace functionspace
} // namespace atlas

