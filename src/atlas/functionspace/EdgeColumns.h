/*
 * (C) Copyright 1996-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef atlas_functionspace_EdgeColumnsFunctionSpace_h
#define atlas_functionspace_EdgeColumnsFunctionSpace_h

#include "eckit/memory/SharedPtr.h"
#include "atlas/mesh/Halo.h"
#include "atlas/field/FieldSet.h"
#include "atlas/functionspace/FunctionSpace.h"
#include "atlas/util/Config.h"

// ----------------------------------------------------------------------------
// Forward declarations

namespace atlas {
namespace mesh {
    class Mesh;
    class HybridElements;
}
}

namespace atlas {
namespace field {
    class FieldSet;
}
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

// ----------------------------------------------------------------------------

class EdgeColumns : public FunctionSpace
{
public:

    EdgeColumns( mesh::Mesh& mesh, const mesh::Halo &, const eckit::Parametrisation & );
    EdgeColumns( mesh::Mesh& mesh, const mesh::Halo & );
    EdgeColumns( mesh::Mesh& mesh );

    virtual ~EdgeColumns();

    virtual std::string name() const { return "Edges"; }

    size_t nb_edges() const;
    size_t nb_edges_global() const; // Only on MPI rank 0, will this be different from 0
    std::vector<size_t> nb_edges_global_foreach_rank() const;

    const mesh::Mesh& mesh() const { return *mesh_.get(); }
          mesh::Mesh& mesh()       { return *mesh_.get(); }

    const mesh::HybridElements& edges() const { return edges_; }
          mesh::HybridElements& edges()       { return edges_; }



// -- Field creation methods

    /// @brief Create a named scalar field
    template< typename DATATYPE > field::Field* createField(
              const std::string& name,
              const eckit::Parametrisation& = util::NoConfig() ) const;

    template< typename DATATYPE > field::Field* createField(
              const std::string& name,
              size_t levels,
              const eckit::Parametrisation& = util::NoConfig() ) const;

    /// @brief Create a named scalar field
    field::Field* createField(
        const std::string& name,
        array::DataType,
        const eckit::Parametrisation& = util::NoConfig() ) const;

    field::Field* createField(
        const std::string& name,
        array::DataType,
        size_t levels,
        const eckit::Parametrisation& = util::NoConfig() ) const;

    /// @brief Create a named field with specified dimensions for the variables
    template< typename DATATYPE >  field::Field* createField(
        const std::string& name,
        const std::vector<size_t>& variables,
        const eckit::Parametrisation& = util::NoConfig() ) const;

    template< typename DATATYPE >  field::Field* createField(
        const std::string& name,
        size_t levels,
        const std::vector<size_t>& variables,
        const eckit::Parametrisation& = util::NoConfig() ) const;

    /// @brief Create a named field with specified dimensions for the variables
    field::Field* createField(const std::string& name,
                              array::DataType,
                              const std::vector<size_t>& variables,
                              const eckit::Parametrisation& = util::NoConfig() ) const;

    field::Field* createField(const std::string& name,
                              array::DataType,
                              size_t levels,
                              const std::vector<size_t>& variables,
                              const eckit::Parametrisation& = util::NoConfig() ) const;

    /// @brief Create a named field based on other field (datatype and dimensioning)
    field::Field* createField(const std::string& name,
                              const field::Field&,
                              const eckit::Parametrisation& = util::NoConfig() ) const;



// -- Parallelisation aware methods

    void haloExchange( field::FieldSet& ) const;
    void haloExchange( field::Field& ) const;
    const parallel::HaloExchange& halo_exchange() const;

    void gather( const field::FieldSet&, field::FieldSet& ) const;
    void gather( const field::Field&, field::Field& ) const;
    const parallel::GatherScatter& gather() const;

    void scatter( const field::FieldSet&, field::FieldSet& ) const;
    void scatter( const field::Field&, field::Field& ) const;
    const parallel::GatherScatter& scatter() const;

    std::string checksum( const field::FieldSet& ) const;
    std::string checksum( const field::Field& ) const;
    const parallel::Checksum& checksum() const;

private: // methods

    std::string gather_scatter_name() const;
    std::string checksum_name() const;
    void constructor();
    size_t config_size(const eckit::Parametrisation& config) const;

private: // data

    eckit::SharedPtr<mesh::Mesh> mesh_; // non-const because functionspace may modify mesh
    mesh::HybridElements& edges_; // non-const because functionspace may modify mesh
    size_t nb_edges_;
    size_t nb_edges_global_;

    eckit::SharedPtr<parallel::GatherScatter> gather_scatter_; // without ghost
    eckit::SharedPtr<parallel::HaloExchange>  halo_exchange_;
    eckit::SharedPtr<parallel::Checksum>      checksum_;
};

// -------------------------------------------------------------------

template< typename DATATYPE >
field::Field* EdgeColumns::createField(
    const std::string& name,
    const eckit::Parametrisation& options ) const
{
    return createField(name,array::DataType::create<DATATYPE>(),options);
}

template< typename DATATYPE >
field::Field* EdgeColumns::createField(
    const std::string& name,
    size_t levels,
    const eckit::Parametrisation& options ) const
{
    return createField(name,array::DataType::create<DATATYPE>(),levels,options);
}

template< typename DATATYPE >
field::Field* EdgeColumns::createField(
    const std::string& name,
    const std::vector<size_t>& variables,
    const eckit::Parametrisation& options ) const
{
    return createField(name,array::DataType::create<DATATYPE>(),variables,options);
}

template< typename DATATYPE >
field::Field* EdgeColumns::createField(
    const std::string& name,
    size_t levels,
    const std::vector<size_t>& variables,
    const eckit::Parametrisation& options ) const
{
    return createField(name,array::DataType::create<DATATYPE>(),levels,variables,options);
}

// -------------------------------------------------------------------------------
#define mesh_Edges mesh::HybridElements
#define Char char
#define GatherScatter parallel::GatherScatter
#define Checksum parallel::Checksum
#define HaloExchange parallel::HaloExchange
#define field_Field field::Field
#define field_FieldSet field::FieldSet
#define mesh_Mesh mesh::Mesh
#define Options eckit::Parametrisation

extern "C" {

EdgeColumns* atlas__functionspace__Edges__new (mesh_Mesh* mesh, int halo);
EdgeColumns* atlas__functionspace__Edges__new_mesh (mesh_Mesh* mesh);
void atlas__functionspace__Edges__delete (EdgeColumns* This);
int atlas__functionspace__Edges__nb_edges(const EdgeColumns* This);
mesh_Mesh* atlas__functionspace__Edges__mesh(EdgeColumns* This);
mesh_Edges* atlas__functionspace__Edges__edges(EdgeColumns* This);
field_Field* atlas__functionspace__Edges__create_field (const EdgeColumns* This, const char* name, int kind, const Options* options);
field_Field* atlas__functionspace__Edges__create_field_vars (const EdgeColumns* This, const char* name, int variables[], int variables_size, int fortran_ordering, int kind, const Options* options);

field_Field* atlas__functionspace__Edges__create_field_lev (const EdgeColumns* This, const char* name, int levels, int kind, const Options* options);
field_Field* atlas__functionspace__Edges__create_field_lev_vars (const EdgeColumns* This, const char* name, int levels, int variables[], int variables_size, int fortran_ordering, int kind, const Options* options);


field_Field* atlas__functionspace__Edges__create_field_template (const EdgeColumns* This, const char* name, const field_Field* field_template, const Options* options);

void atlas__functionspace__Edges__halo_exchange_fieldset(const EdgeColumns* This, field_FieldSet* fieldset);
void atlas__functionspace__Edges__halo_exchange_field(const EdgeColumns* This, field_Field* field);
const HaloExchange* atlas__functionspace__Edges__get_halo_exchange(const EdgeColumns* This);

void atlas__functionspace__Edges__gather_fieldset(const EdgeColumns* This, const field_FieldSet* local, field_FieldSet* global);
void atlas__functionspace__Edges__gather_field(const EdgeColumns* This, const field_Field* local, field_Field* global);
const GatherScatter* atlas__functionspace__Edges__get_gather(const EdgeColumns* This);

void atlas__functionspace__Edges__scatter_fieldset(const EdgeColumns* This, const field_FieldSet* global, field_FieldSet* local);
void atlas__functionspace__Edges__scatter_field(const EdgeColumns* This, const field_Field* global, field_Field* local);
const GatherScatter* atlas__functionspace__Edges__get_scatter(const EdgeColumns* This);

void atlas__functionspace__Edges__checksum_fieldset(const EdgeColumns* This, const field_FieldSet* fieldset, Char* &checksum, int &size, int &allocated);
void atlas__functionspace__Edges__checksum_field(const EdgeColumns* This, const field_Field* field, Char* &checksum, int &size, int &allocated);
const Checksum* atlas__functionspace__Edges__get_checksum(const EdgeColumns* This);
}

#undef mesh_Edges
#undef Char
#undef GatherScatter
#undef Checksum
#undef HaloExchange
#undef field_Field
#undef field_FieldSet
#undef mesh_Mesh
#undef Options

} // namespace functionspace
} // namespace atlas

#endif // atlas_functionspace_EdgeColumnsFunctionSpace_h
