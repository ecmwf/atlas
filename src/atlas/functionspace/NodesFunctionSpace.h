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
#include "atlas/FieldSet.h"

namespace atlas { class Mesh; }
namespace atlas { class FieldSet; }

namespace atlas {
namespace functionspace {

// -------------------------------------------------------------------

class Halo
{
public:
  Halo(const size_t size) : size_(size) {}
  size_t size() const { return size_; }
private:
  size_t size_;
};

// -------------------------------------------------------------------

class NodesFunctionSpace : public next::FunctionSpace
{
public:

  NodesFunctionSpace(const std::string& name, Mesh& mesh, const Halo& = Halo(0) );

  virtual ~NodesFunctionSpace();

  /// @brief Create a scalar field
  template< typename DATATYPE >
  Field* createField(const std::string& name) const;

  /// @brief Create a field with specified dimensions for the variables
  template< typename DATATYPE >
  Field* createField(const std::string& name, const std::vector<size_t>& variables) const;

  /// @brief Create a field based on other field (datatype and dimensioning)
  Field* createField(const std::string& name, const Field&) const;

  /// @brief Create a global scalar field
  template< typename DATATYPE >
  Field* createGlobalField(const std::string& name) const;

  /// @brief Create a global field with specified dimensions for the variables
  template< typename DATATYPE >
  Field* createGlobalField(const std::string& name, const std::vector<size_t>& variables) const;

  /// @brief Create a global field based on other field (datatype and dimensioning)
  Field* createGlobalField(const std::string& name, const Field&) const;

  void haloExchange( FieldSet& ) const;
  void haloExchange( Field& ) const;

  void gather( const FieldSet&, FieldSet& ) const;
  void gather( const Field&, Field& ) const;

  void scatter( const FieldSet&, FieldSet& ) const;
  void scatter( const Field&, Field& ) const;

  std::string checksum( const FieldSet& ) const;
  std::string checksum( const Field& ) const;

  size_t nb_nodes() const;
  size_t nb_nodes_global() const;

  const Mesh& mesh() const { return mesh_; }
        Mesh& mesh()       { return mesh_; }

private: // methods

  std::string halo_name() const;
  std::string gather_scatter_name() const;
  std::string checksum_name() const;

private: // data

  Mesh& mesh_; // non-const because functionspace may modify mesh
  size_t halo_;
  size_t nb_nodes_;
  size_t nb_nodes_global_;
};

// -------------------------------------------------------------------

class NodesColumnFunctionSpace : public NodesFunctionSpace
{
public:

  NodesColumnFunctionSpace(const std::string& name, Mesh& mesh, size_t nb_levels, const Halo& = Halo(0) );

  virtual ~NodesColumnFunctionSpace();

  /// @brief Create a scalar field
  template< typename DATATYPE >
  Field* createField(const std::string& name) const;

  /// @brief Create a vector field
  template< typename DATATYPE >
  Field* createField(const std::string& name, size_t var1) const;

  /// @brief Create a tensor field
  template< typename DATATYPE >
  Field* createField(const std::string& name, size_t var1, size_t var2) const;

  /// @brief Create a field with specified dimensions for the variables
  template< typename DATATYPE >
  Field* createField(const std::string& name, const std::vector<size_t>& variables) const;

  /// @brief Create a global scalar field
  template< typename DATATYPE >
  Field* createGlobalField(const std::string& name) const;

  /// @brief Create a global vector field
  template< typename DATATYPE >
  Field* createGlobalField(const std::string& name, size_t var1) const;

  /// @brief Create a global tensor field
  template< typename DATATYPE >
  Field* createGlobalField(const std::string& name, size_t var1, size_t var2) const;

  /// @brief Create a global field with specified dimensions for the variables
  template< typename DATATYPE >
  Field* createGlobalField(const std::string& name, const std::vector<size_t>& variables) const;

  /// @brief Create a global field based on other field (datatype and dimensioning)
  Field* createGlobalField(const std::string& name, const Field&) const;

  size_t nb_levels() const;

private: // data

  size_t nb_levels_;
};

} // namespace functionspace
} // namespace atlas

#endif // atlas_functionspace_NodesFunctionSpace_h
