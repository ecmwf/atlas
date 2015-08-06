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

  size_t nb_nodes() const;
  size_t nb_nodes_global() const; // Only on MPI rank 0, will this be different 0

  const Mesh& mesh() const { return mesh_; }
        Mesh& mesh()       { return mesh_; }

  const Nodes& nodes() const { return nodes_; }
        Nodes& nodes()       { return nodes_; }



// -- Local Field creation methods

  /// @brief Create an unnamed scalar field
  template< typename DATATYPE >
  Field* createField() const;

  /// @brief Create an unnamed scalar field
  Field* createField(const std::string& datatype) const;

  /// @brief Create a named scalar field
  template< typename DATATYPE >
  Field* createField(const std::string& name) const;

  /// @brief Create a named scalar field
  Field* createField(const std::string& name, const std::string& datatype) const;

  /// @brief Create an unnamed field with specified dimensions for the variables
  template< typename DATATYPE >
  Field* createField(const std::vector<size_t>& variables) const;

  /// @brief Create an unnamed field with specified dimensions for the variables
  Field* createField(const std::vector<size_t>& variables, const std::string& datatype) const;

  /// @brief Create a named field with specified dimensions for the variables
  template< typename DATATYPE >
  Field* createField(const std::string& name, const std::vector<size_t>& variables) const;

  /// @brief Create a named field with specified dimensions for the variables
  Field* createField(const std::string& name, const std::vector<size_t>& variables, const std::string& datatype) const;

  /// @brief Create an unnamed field based on other field (datatype and dimensioning)
  Field* createField(const Field&) const;

  /// @brief Create a named field based on other field (datatype and dimensioning)
  Field* createField(const std::string& name, const Field&) const;



// -- Global Field creation methods

  /// @brief Create an unnamed global scalar field
  template< typename DATATYPE >
  Field* createGlobalField() const;

  /// @brief Create an unnamed global scalar field
  Field* createGlobalField(const std::string& datatype) const;

  /// @brief Create a named global scalar field
  template< typename DATATYPE >
  Field* createGlobalField(const std::string& name) const;

  /// @brief Create a named global scalar field
  Field* createGlobalField(const std::string& name, const std::string& datatype) const;

  /// @brief Create an unnamed global field with specified dimensions for the variables
  template< typename DATATYPE >
  Field* createGlobalField(const std::vector<size_t>& variables) const;

  /// @brief Create an unnamed global field with specified dimensions for the variables
  Field* createGlobalField(const std::vector<size_t>& variables, const std::string& datatype) const;

  /// @brief Create a named global field with specified dimensions for the variables
  template< typename DATATYPE >
  Field* createGlobalField(const std::string& name, const std::vector<size_t>& variables) const;

  /// @brief Create a named field with specified dimensions for the variables
  Field* createGlobalField(const std::string& name, const std::vector<size_t>& variables, const std::string& datatype) const;

  /// @brief Create an unnamed global field based on other field (datatype and dimensioning)
  Field* createGlobalField(const Field&) const;

  /// @brief Create a named global field based on other field (datatype and dimensioning)
  Field* createGlobalField(const std::string& name, const Field&) const;



// -- Parallelisation aware methods

  void haloExchange( FieldSet& ) const;
  void haloExchange( Field& ) const;

  void gather( const FieldSet&, FieldSet& ) const;
  void gather( const Field&, Field& ) const;

  void scatter( const FieldSet&, FieldSet& ) const;
  void scatter( const Field&, Field& ) const;

  std::string checksum( const FieldSet& ) const;
  std::string checksum( const Field& ) const;

  template< typename DATATYPE >
  void sum( const Field&, DATATYPE& ) const;

  template< typename DATATYPE >
  void sum( const Field&, std::vector<DATATYPE>& ) const;

  template< typename DATATYPE >
  void maximum( const Field&, DATATYPE& ) const;

  template< typename DATATYPE >
  void maximum( const Field&, std::vector<DATATYPE>& ) const;

  template< typename DATATYPE >
  void minimum( const Field&, DATATYPE& ) const;

  template< typename DATATYPE >
  void minimum( const Field&, std::vector<DATATYPE>& ) const;

  template< typename DATATYPE >
  void mean( const Field&, DATATYPE& ) const;

  template< typename DATATYPE >
  void mean( const Field&, std::vector<DATATYPE>& ) const;

  template< typename DATATYPE >
  void mean( const Field&, const DATATYPE& sum, DATATYPE& ) const;

  template< typename DATATYPE >
  void standard_deviation( const Field&, DATATYPE& ) const;

  template< typename DATATYPE >
  void standard_deviation( const Field&, std::vector<DATATYPE>& ) const;

  template< typename DATATYPE >
  void standard_deviation( const Field&, const DATATYPE& mean, DATATYPE& ) const;

private: // methods

  std::string halo_name() const;
  std::string gather_scatter_name() const;
  std::string checksum_name() const;

private: // data

  Mesh& mesh_; // non-const because functionspace may modify mesh
  Nodes& nodes_; // non-const because functionspace may modify mesh
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

  /// @brief Create a field with specified dimensions for the variables
  template< typename DATATYPE >
  Field* createField(const std::string& name, const std::vector<size_t>& variables) const;

  /// @brief Create a global scalar field
  template< typename DATATYPE >
  Field* createGlobalField(const std::string& name) const;

  /// @brief Create a global field with specified dimensions for the variables
  template< typename DATATYPE >
  Field* createGlobalField(const std::string& name, const std::vector<size_t>& variables) const;

  /// @brief Create a global field based on other field (datatype and dimensioning)
  Field* createGlobalField(const std::string& name, const Field&) const;

  /// @brief Create a global field based on other field (datatype and dimensioning)
  Field* createGlobalField(const Field&) const;

  size_t nb_levels() const;

private: // data

  size_t nb_levels_;
};

// -------------------------------------------------------------------

template< typename DATATYPE >
Field* NodesFunctionSpace::createField() const
{
  return createField(DataType::datatype<DATATYPE>());
}

template< typename DATATYPE >
Field* NodesFunctionSpace::createField(const std::string& name) const
{
  return createField(name,DataType::datatype<DATATYPE>());
}

template< typename DATATYPE >
Field* NodesFunctionSpace::createField(const std::vector<size_t>& variables) const
{
  return createField(variables,DataType::datatype<DATATYPE>());
}

template< typename DATATYPE >
Field* NodesFunctionSpace::createField(const std::string& name, const std::vector<size_t>& variables) const
{
  return createField(name,variables,DataType::datatype<DATATYPE>());
}


template< typename DATATYPE >
Field* NodesFunctionSpace::createGlobalField() const
{
  return createGlobalField(DataType::datatype<DATATYPE>());
}

template< typename DATATYPE >
Field* NodesFunctionSpace::createGlobalField(const std::string& name) const
{
  return createGlobalField(name,DataType::datatype<DATATYPE>());
}

template< typename DATATYPE >
Field* NodesFunctionSpace::createGlobalField(const std::vector<size_t>& variables) const
{
  return createGlobalField(variables,DataType::datatype<DATATYPE>());
}

template< typename DATATYPE >
Field* NodesFunctionSpace::createGlobalField(const std::string& name, const std::vector<size_t>& variables) const
{
  return createGlobalField(name,variables,DataType::datatype<DATATYPE>());
}

extern "C" {
NodesFunctionSpace* atlas__NodesFunctionSpace__new (const char* name, Mesh* mesh, int halo);
void atlas__NodesFunctionSpace__delete (NodesFunctionSpace* This);
NodesColumnFunctionSpace* atlas__NodesColumnFunctionSpace__new (const char* name, Mesh* mesh, int nb_levels, int halo);
void atlas__NodesColumnFunctionSpace__delete (NodesColumnFunctionSpace* This);
}


} // namespace functionspace
} // namespace atlas

#endif // atlas_functionspace_NodesFunctionSpace_h
