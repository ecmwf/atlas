/*
 * (C) Copyright 1996-2014 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "atlas/Mesh.h"
#include "atlas/functionspace/NodesFunctionSpace.h"
#include "atlas/field/FieldT.h"

namespace atlas {
namespace functionspace {

NodesFunctionSpace::NodesFunctionSpace(const std::string& name, Mesh& mesh)
  : next::FunctionSpace(name),
    mesh_(mesh)
{
}

NodesFunctionSpace::~NodesFunctionSpace() {}

size_t NodesFunctionSpace::nb_nodes() const
{
  return mesh_.function_space("nodes").shape(0);
}


template<> Field* NodesFunctionSpace::create_field<int>(const std::string& name) {
  return new field::FieldT<int>(name, make_shape(nb_nodes()));
}

template<> Field* NodesFunctionSpace::create_field<long>(const std::string& name) {
  return new field::FieldT<long>(name, make_shape(nb_nodes()) );
}

template<> Field* NodesFunctionSpace::create_field<float>(const std::string& name) {
  return new field::FieldT<float>(name, make_shape(nb_nodes()) );
}

template<> Field* NodesFunctionSpace::create_field<double>(const std::string& name) {
  return new field::FieldT<double>(name, make_shape(nb_nodes()) );
}

template<> Field* NodesFunctionSpace::create_field<int>(const std::string& name, const size_t var1) {
  return new field::FieldT<int>(name, make_shape(nb_nodes(),var1) );
}

template<> Field* NodesFunctionSpace::create_field<long>(const std::string& name, const size_t var1) {
  return new field::FieldT<long>(name, make_shape(nb_nodes(),var1) );
}

template<> Field* NodesFunctionSpace::create_field<float>(const std::string& name, const size_t var1) {
  return new field::FieldT<float>(name, make_shape(nb_nodes(),var1) );
}

template<> Field* NodesFunctionSpace::create_field<double>(const std::string& name, const size_t var1) {
  return new field::FieldT<double>(name, make_shape(nb_nodes(),var1) );
}

template<> Field* NodesFunctionSpace::create_field<int>(const std::string& name, const size_t var1, const size_t var2) {
  return new field::FieldT<int>(name, make_shape(nb_nodes(),var1,var2) );
}

template<> Field* NodesFunctionSpace::create_field<long>(const std::string& name, const size_t var1, const size_t var2) {
  return new field::FieldT<long>(name, make_shape(nb_nodes(),var1,var2) );
}

template<> Field* NodesFunctionSpace::create_field<float>(const std::string& name, const size_t var1, const size_t var2) {
  return new field::FieldT<float>(name, make_shape(nb_nodes(),var1,var2) );
}

template<> Field* NodesFunctionSpace::create_field<double>(const std::string& name, const size_t var1, const size_t var2) {
  return new field::FieldT<double>(name, make_shape(nb_nodes(),var1,var2) );
}



// ----------------------------------------------------------------------

NodesColumnFunctionSpace::NodesColumnFunctionSpace(const std::string& name, Mesh& mesh, const size_t levels)
  : next::FunctionSpace(name),
    mesh_(mesh),
    levels_(levels)
{
}

NodesColumnFunctionSpace::~NodesColumnFunctionSpace() {}

size_t NodesColumnFunctionSpace::nb_nodes() const
{
  return mesh_.function_space("nodes").shape(0);
}

template<> Field* NodesColumnFunctionSpace::create_field<int>(const std::string& name) {
  return new field::FieldT<int>(name, make_shape(nb_nodes(),levels_) );
}

template<> Field* NodesColumnFunctionSpace::create_field<long>(const std::string& name) {
  return new field::FieldT<long>(name, make_shape(nb_nodes(),levels_) );
}

template<> Field* NodesColumnFunctionSpace::create_field<float>(const std::string& name) {
  return new field::FieldT<float>(name, make_shape(nb_nodes(),levels_) );
}

template<> Field* NodesColumnFunctionSpace::create_field<double>(const std::string& name) {
  return new field::FieldT<double>(name, make_shape(nb_nodes(),levels_) );
}

template<> Field* NodesColumnFunctionSpace::create_field<int>(const std::string& name, const size_t var1) {
  return new field::FieldT<int>(name, make_shape(nb_nodes(),levels_,var1) );
}

template<> Field* NodesColumnFunctionSpace::create_field<long>(const std::string& name, const size_t var1) {
  return new field::FieldT<long>(name, make_shape(nb_nodes(),levels_,var1) );
}

template<> Field* NodesColumnFunctionSpace::create_field<float>(const std::string& name, const size_t var1) {
  return new field::FieldT<float>(name, make_shape(nb_nodes(),levels_,var1) );
}

template<> Field* NodesColumnFunctionSpace::create_field<double>(const std::string& name, const size_t var1) {
  return new field::FieldT<double>(name, make_shape(nb_nodes(),levels_,var1) );
}

template<> Field* NodesColumnFunctionSpace::create_field<int>(const std::string& name, const size_t var1, const size_t var2) {
  return new field::FieldT<int>(name, make_shape(nb_nodes(),levels_,var1,var2) );
}

template<> Field* NodesColumnFunctionSpace::create_field<long>(const std::string& name, const size_t var1, const size_t var2) {
  return new field::FieldT<long>(name, make_shape(nb_nodes(),levels_,var1,var2) );
}

template<> Field* NodesColumnFunctionSpace::create_field<float>(const std::string& name, const size_t var1, const size_t var2) {
  return new field::FieldT<float>(name, make_shape(nb_nodes(),levels_,var1,var2) );
}

template<> Field* NodesColumnFunctionSpace::create_field<double>(const std::string& name, const size_t var1, const size_t var2) {
  return new field::FieldT<double>(name, make_shape(nb_nodes(),levels_,var1,var2) );
}

// ----------------------------------------------------------------------

} // namespace functionspace
} // namespace atlas

