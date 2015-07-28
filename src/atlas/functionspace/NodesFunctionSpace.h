/*
 * (C) Copyright 1996-2015 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef atlas_functionspace_NodesFunctionSpace_h
#define atlas_functionspace_NodesFunctionSpace_h


#include "eckit/memory/SharedPtr.h"
#include "atlas/FunctionSpace.h"

namespace atlas { class Mesh; }

namespace atlas {
namespace functionspace {

// -------------------------------------------------------------------


class NodesFunctionSpace : public next::FunctionSpace
{
public:

  NodesFunctionSpace(const std::string& name, Mesh& mesh);

  virtual ~NodesFunctionSpace();

  /// @brief Create a scalar field
  template< typename DATATYPE >
  Field* create_field(const std::string& name);

  /// @brief Create a vector field
  template< typename DATATYPE >
  Field* create_field(const std::string& name, const size_t var1);

  /// @brief Create a tensor field
  template< typename DATATYPE >
  Field* create_field(const std::string& name, const size_t var1, const size_t var2);

private: // methods

  size_t nb_nodes() const;

private: // data

  Mesh& mesh_; // non-const because functionspace may modify mesh
};

// -------------------------------------------------------------------

class NodesColumnFunctionSpace : public next::FunctionSpace
{
public:

  NodesColumnFunctionSpace(const std::string& name, Mesh& mesh, const size_t levels);

  virtual ~NodesColumnFunctionSpace();

  /// @brief Create a scalar field
  template< typename DATATYPE >
  Field* create_field(const std::string& name);

  /// @brief Create a vector field
  template< typename DATATYPE >
  Field* create_field(const std::string& name, const size_t var1);

  /// @brief Create a tensor field
  template< typename DATATYPE >
  Field* create_field(const std::string& name, const size_t var1, const size_t var2);

private: // methods

  size_t nb_nodes() const;

private: // data

  Mesh& mesh_; // non-const because functionspace may modify mesh
  size_t levels_;
};

} // namespace functionspace
} // namespace atlas

#endif // atlas_functionspace_NodesFunctionSpace_h
