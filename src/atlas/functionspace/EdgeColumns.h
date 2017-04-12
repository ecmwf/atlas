/*
 * (C) Copyright 1996-2017 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#pragma once

#include "eckit/memory/SharedPtr.h"
#include "atlas/mesh/Halo.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/field/FieldSet.h"
#include "atlas/functionspace/FunctionSpace.h"
#include "atlas/util/Config.h"

// ----------------------------------------------------------------------------
// Forward declarations

namespace atlas {
  class FieldSet;
}

namespace atlas {
namespace parallel {
    class HaloExchange;
    class GatherScatter;
    class Checksum;
}
}

namespace atlas {
namespace functionspace {
namespace detail {

// ----------------------------------------------------------------------------

class EdgeColumns : public FunctionSpaceImpl
{
public:

    EdgeColumns( const mesh::Mesh&, const mesh::Halo &, const eckit::Parametrisation & );
    EdgeColumns( const mesh::Mesh&, const mesh::Halo & );
    EdgeColumns( const mesh::Mesh& );

    virtual ~EdgeColumns();

    virtual std::string name() const { return "Edges"; }

    size_t nb_edges() const;
    size_t nb_edges_global() const; // Only on MPI rank 0, will this be different from 0
    std::vector<size_t> nb_edges_global_foreach_rank() const;

    const mesh::Mesh& mesh() const { return mesh_; }
          mesh::Mesh& mesh()       { return mesh_; }

    const mesh::HybridElements& edges() const { return edges_; }
          mesh::HybridElements& edges()       { return edges_; }



// -- Field creation methods

    /// @brief Create a named scalar field
    template< typename DATATYPE > Field createField(
              const std::string& name,
              const eckit::Parametrisation& = util::NoConfig() ) const;

    template< typename DATATYPE > Field createField(
              const std::string& name,
              size_t levels,
              const eckit::Parametrisation& = util::NoConfig() ) const;

    /// @brief Create a named scalar field
    Field createField(
        const std::string& name,
        array::DataType,
        const eckit::Parametrisation& = util::NoConfig() ) const;

    Field createField(
        const std::string& name,
        array::DataType,
        size_t levels,
        const eckit::Parametrisation& = util::NoConfig() ) const;

    /// @brief Create a named field with specified dimensions for the variables
    template< typename DATATYPE >  Field createField(
        const std::string& name,
        const std::vector<size_t>& variables,
        const eckit::Parametrisation& = util::NoConfig() ) const;

    template< typename DATATYPE >  Field createField(
        const std::string& name,
        size_t levels,
        const std::vector<size_t>& variables,
        const eckit::Parametrisation& = util::NoConfig() ) const;

    /// @brief Create a named field with specified dimensions for the variables
    Field createField(const std::string& name,
                              array::DataType,
                              const std::vector<size_t>& variables,
                              const eckit::Parametrisation& = util::NoConfig() ) const;

    Field createField(const std::string& name,
                              array::DataType,
                              size_t levels,
                              const std::vector<size_t>& variables,
                              const eckit::Parametrisation& = util::NoConfig() ) const;

    /// @brief Create a named field based on other field (datatype and dimensioning)
    Field createField(const std::string& name,
                              const Field&,
                              const eckit::Parametrisation& = util::NoConfig() ) const;



// -- Parallelisation aware methods

    void haloExchange( FieldSet& ) const;
    void haloExchange( Field& ) const;
    const parallel::HaloExchange& halo_exchange() const;

    void gather( const FieldSet&, FieldSet& ) const;
    void gather( const Field&, Field& ) const;
    const parallel::GatherScatter& gather() const;

    void scatter( const FieldSet&, FieldSet& ) const;
    void scatter( const Field&, Field& ) const;
    const parallel::GatherScatter& scatter() const;

    std::string checksum( const FieldSet& ) const;
    std::string checksum( const Field& ) const;
    const parallel::Checksum& checksum() const;

private: // methods

    std::string gather_scatter_name() const;
    std::string checksum_name() const;
    void constructor();
    size_t config_size(const eckit::Parametrisation& config) const;
    size_t footprint() const;

private: // data

    mesh::Mesh mesh_; // non-const because functionspace may modify mesh
    mesh::HybridElements& edges_; // non-const because functionspace may modify mesh
    size_t nb_edges_;
    size_t nb_edges_global_;

    eckit::SharedPtr<parallel::GatherScatter> gather_scatter_; // without ghost
    eckit::SharedPtr<parallel::HaloExchange>  halo_exchange_;
    eckit::SharedPtr<parallel::Checksum>      checksum_;
};

// -------------------------------------------------------------------

template< typename DATATYPE >
Field EdgeColumns::createField(
    const std::string& name,
    const eckit::Parametrisation& options ) const
{
    return createField(name,array::DataType::create<DATATYPE>(),options);
}

template< typename DATATYPE >
Field EdgeColumns::createField(
    const std::string& name,
    size_t levels,
    const eckit::Parametrisation& options ) const
{
    return createField(name,array::DataType::create<DATATYPE>(),levels,options);
}

template< typename DATATYPE >
Field EdgeColumns::createField(
    const std::string& name,
    const std::vector<size_t>& variables,
    const eckit::Parametrisation& options ) const
{
    return createField(name,array::DataType::create<DATATYPE>(),variables,options);
}

template< typename DATATYPE >
Field EdgeColumns::createField(
    const std::string& name,
    size_t levels,
    const std::vector<size_t>& variables,
    const eckit::Parametrisation& options ) const
{
    return createField(name,array::DataType::create<DATATYPE>(),levels,variables,options);
}

// -------------------------------------------------------------------------------

extern "C" {

EdgeColumns* atlas__functionspace__Edges__new (mesh::Mesh::Implementation* mesh, int halo);
EdgeColumns* atlas__functionspace__Edges__new_mesh (mesh::Mesh::Implementation* mesh);
void atlas__functionspace__Edges__delete (EdgeColumns* This);
int atlas__functionspace__Edges__nb_edges(const EdgeColumns* This);
mesh::Mesh::Implementation* atlas__functionspace__Edges__mesh(EdgeColumns* This);
mesh::Edges* atlas__functionspace__Edges__edges(EdgeColumns* This);
field::FieldImpl* atlas__functionspace__Edges__create_field (const EdgeColumns* This, const char* name, int kind, const eckit::Parametrisation* options);
field::FieldImpl* atlas__functionspace__Edges__create_field_vars (const EdgeColumns* This, const char* name, int variables[], int variables_size, int fortran_ordering, int kind, const eckit::Parametrisation* options);

field::FieldImpl* atlas__functionspace__Edges__create_field_lev (const EdgeColumns* This, const char* name, int levels, int kind, const eckit::Parametrisation* options);
field::FieldImpl* atlas__functionspace__Edges__create_field_lev_vars (const EdgeColumns* This, const char* name, int levels, int variables[], int variables_size, int fortran_ordering, int kind, const eckit::Parametrisation* options);


field::FieldImpl* atlas__functionspace__Edges__create_field_template (const EdgeColumns* This, const char* name, const field::FieldImpl* field_template, const eckit::Parametrisation* options);

void atlas__functionspace__Edges__halo_exchange_fieldset(const EdgeColumns* This, field::FieldSetImpl* fieldset);
void atlas__functionspace__Edges__halo_exchange_field(const EdgeColumns* This, field::FieldImpl* field);
const parallel::HaloExchange* atlas__functionspace__Edges__get_halo_exchange(const EdgeColumns* This);

void atlas__functionspace__Edges__gather_fieldset(const EdgeColumns* This, const field::FieldSetImpl* local, field::FieldSetImpl* global);
void atlas__functionspace__Edges__gather_field(const EdgeColumns* This, const field::FieldImpl* local, field::FieldImpl* global);
const parallel::GatherScatter* atlas__functionspace__Edges__get_gather(const EdgeColumns* This);

void atlas__functionspace__Edges__scatter_fieldset(const EdgeColumns* This, const field::FieldSetImpl* global, field::FieldSetImpl* local);
void atlas__functionspace__Edges__scatter_field(const EdgeColumns* This, const field::FieldImpl* global, field::FieldImpl* local);
const parallel::GatherScatter* atlas__functionspace__Edges__get_scatter(const EdgeColumns* This);

void atlas__functionspace__Edges__checksum_fieldset(const EdgeColumns* This, const field::FieldSetImpl* fieldset, char* &checksum, int &size, int &allocated);
void atlas__functionspace__Edges__checksum_field(const EdgeColumns* This, const field::FieldImpl* field, char* &checksum, int &size, int &allocated);
const parallel::Checksum* atlas__functionspace__Edges__get_checksum(const EdgeColumns* This);
}

} // namespace detail

// -------------------------------------------------------------------

class EdgeColumns : public FunctionSpace {

public:

    EdgeColumns();
    EdgeColumns( const FunctionSpace& );
    EdgeColumns( const mesh::Mesh&, const mesh::Halo&, const eckit::Parametrisation& );
    EdgeColumns( const mesh::Mesh& mesh, const mesh::Halo& );
    EdgeColumns( const mesh::Mesh& mesh );

    operator bool() const { return valid(); }
    bool valid() const { return functionspace_; }

    size_t nb_edges() const;
    size_t nb_edges_global() const; // Only on MPI rank 0, will this be different from 0

    const mesh::Mesh& mesh() const;

    const mesh::HybridElements& edges() const;

// -- Field creation methods

    /// @brief Create a named scalar field
    template< typename DATATYPE > Field createField(
              const std::string& name,
              const eckit::Parametrisation& = util::NoConfig() ) const;

    template< typename DATATYPE > Field createField(
              const std::string& name,
              size_t levels,
              const eckit::Parametrisation& = util::NoConfig() ) const;

    /// @brief Create a named scalar field
    Field createField(
        const std::string& name,
        array::DataType,
        const eckit::Parametrisation& = util::NoConfig() ) const;

    Field createField(
        const std::string& name,
        array::DataType,
        size_t levels,
        const eckit::Parametrisation& = util::NoConfig() ) const;

    /// @brief Create a named field with specified dimensions for the variables
    template< typename DATATYPE >  Field createField(
        const std::string& name,
        const std::vector<size_t>& variables,
        const eckit::Parametrisation& = util::NoConfig() ) const;

    template< typename DATATYPE >  Field createField(
        const std::string& name,
        size_t levels,
        const std::vector<size_t>& variables,
        const eckit::Parametrisation& = util::NoConfig() ) const;

    /// @brief Create a named field with specified dimensions for the variables
    Field createField(const std::string& name,
                              array::DataType,
                              const std::vector<size_t>& variables,
                              const eckit::Parametrisation& = util::NoConfig() ) const;

    Field createField(const std::string& name,
                              array::DataType,
                              size_t levels,
                              const std::vector<size_t>& variables,
                              const eckit::Parametrisation& = util::NoConfig() ) const;

    /// @brief Create a named field based on other field (datatype and dimensioning)
    Field createField(const std::string& name,
                              const Field&,
                              const eckit::Parametrisation& = util::NoConfig() ) const;



// -- Parallelisation aware methods

    void haloExchange( FieldSet& ) const;
    void haloExchange( Field& ) const;
    const parallel::HaloExchange& halo_exchange() const;

    void gather( const FieldSet&, FieldSet& ) const;
    void gather( const Field&, Field& ) const;
    const parallel::GatherScatter& gather() const;

    void scatter( const FieldSet&, FieldSet& ) const;
    void scatter( const Field&, Field& ) const;
    const parallel::GatherScatter& scatter() const;

    std::string checksum( const FieldSet& ) const;
    std::string checksum( const Field& ) const;
    const parallel::Checksum& checksum() const;

private:

  const detail::EdgeColumns* functionspace_;
};

template< typename DATATYPE >
Field EdgeColumns::createField(
          const std::string& name,
          const eckit::Parametrisation& config ) const {
  return functionspace_->createField<DATATYPE>(name,config);
}

template< typename DATATYPE >
Field EdgeColumns::createField(
          const std::string& name,
          size_t levels,
          const eckit::Parametrisation& config ) const {
  return functionspace_->createField<DATATYPE>(name,levels,config);
}


template< typename DATATYPE >
Field EdgeColumns::createField(
    const std::string& name,
    const std::vector<size_t>& variables,
    const eckit::Parametrisation& config ) const {
  return functionspace_->createField<DATATYPE>(name,variables,config);
}

template< typename DATATYPE >
Field EdgeColumns::createField(
    const std::string& name,
    size_t levels,
    const std::vector<size_t>& variables,
    const eckit::Parametrisation& config ) const {
  return functionspace_->createField<DATATYPE>(name,levels,variables,config);
}


} // namespace functionspace
} // namespace atlas

