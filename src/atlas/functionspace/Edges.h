/*
 * (C) Copyright 1996-2015 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef atlas_functionspace_EdgesFunctionSpace_h
#define atlas_functionspace_EdgesFunctionSpace_h


#include "eckit/memory/SharedPtr.h"
#include "atlas/FunctionSpace.h"
#include "atlas/FieldSet.h"
#include "atlas/functionspace/Nodes.h"

namespace atlas { class Mesh; }
namespace atlas { namespace mesh { class HybridElements; } }
namespace atlas { class FieldSet; }

namespace atlas {
namespace functionspace {

// -------------------------------------------------------------------

class Edges : public FunctionSpace
{
public:

    Edges( Mesh& mesh, const Halo &, const eckit::Parametrisation & );
    Edges( Mesh& mesh, const Halo & );
    Edges( Mesh& mesh );

    virtual ~Edges();

    virtual std::string name() const { return "Edges"; }

    size_t nb_edges() const;
    size_t nb_edges_global() const; // Only on MPI rank 0, will this be different from 0
    std::vector<size_t> nb_edges_global_foreach_rank() const;

    const atlas::Mesh& mesh() const { return *mesh_.get(); }
          atlas::Mesh& mesh()       { return *mesh_.get(); }

    const mesh::HybridElements& edges() const { return edges_; }
          mesh::HybridElements& edges()       { return edges_; }



// -- Local Field creation methods

    /// @brief Create a named scalar field
    template< typename DATATYPE > Field* createField(const std::string& name) const;
    template< typename DATATYPE > Field* createField(const std::string& name, size_t levels) const;

    /// @brief Create a named scalar field
    Field* createField(const std::string& name, DataType) const;
    Field* createField(const std::string& name, DataType, size_t levels) const;

    /// @brief Create a named field with specified dimensions for the variables
    template< typename DATATYPE >  Field* createField(const std::string& name, const std::vector<size_t>& variables) const;
    template< typename DATATYPE >  Field* createField(const std::string& name, size_t levels, const std::vector<size_t>& variables) const;

    /// @brief Create a named field with specified dimensions for the variables
    Field* createField(const std::string& name, DataType, const std::vector<size_t>& variables) const;
    Field* createField(const std::string& name, DataType, size_t levels, const std::vector<size_t>& variables) const;

    /// @brief Create a named field based on other field (datatype and dimensioning)
    Field* createField(const std::string& name, const Field&) const;



// -- Global Field creation methods

    /// @brief Create a named global scalar field
    template< typename DATATYPE >  Field* createGlobalField(const std::string& name) const;
    template< typename DATATYPE >  Field* createGlobalField(const std::string& name,size_t levels) const;

    /// @brief Create a named global scalar field
    Field* createGlobalField(const std::string& name, DataType) const;
    Field* createGlobalField(const std::string& name, DataType, size_t levels) const;

    /// @brief Create a named global field with specified dimensions for the variables
    template< typename DATATYPE >  Field* createGlobalField(const std::string& name, const std::vector<size_t>& variables) const;
    template< typename DATATYPE >  Field* createGlobalField(const std::string& name, size_t levels, const std::vector<size_t>& variables) const;

    /// @brief Create a named field with specified dimensions for the variables
    Field* createGlobalField(const std::string& name, DataType, const std::vector<size_t>& variables) const;
    Field* createGlobalField(const std::string& name, DataType, size_t levels, const std::vector<size_t>& variables) const;

    /// @brief Create a named global field based on other field (datatype and dimensioning)
    Field* createGlobalField(const std::string& name, const Field&) const;



// -- Parallelisation aware methods

    void haloExchange( FieldSet& ) const;
    void haloExchange( Field& ) const;
    const mpl::HaloExchange& halo_exchange() const;

    void gather( const FieldSet&, FieldSet& ) const;
    void gather( const Field&, Field& ) const;
    const mpl::GatherScatter& gather() const;

    void scatter( const FieldSet&, FieldSet& ) const;
    void scatter( const Field&, Field& ) const;
    const mpl::GatherScatter& scatter() const;

    std::string checksum( const FieldSet& ) const;
    std::string checksum( const Field& ) const;
    const mpl::Checksum& checksum() const;

private: // methods

    std::string gather_scatter_name() const;
    std::string checksum_name() const;
    void constructor();

private: // data

    eckit::SharedPtr<Mesh> mesh_; // non-const because functionspace may modify mesh
    mesh::HybridElements& edges_; // non-const because functionspace may modify mesh
    size_t nb_edges_;
    size_t nb_edges_global_;

    eckit::SharedPtr<mpl::GatherScatter> gather_scatter_; // without ghost
    eckit::SharedPtr<mpl::HaloExchange>  halo_exchange_;
    eckit::SharedPtr<mpl::Checksum>      checksum_;
};

// -------------------------------------------------------------------

template< typename DATATYPE >
Field* Edges::createField(const std::string& name) const
{
    return createField(name,DataType::create<DATATYPE>());
}

template< typename DATATYPE >
Field* Edges::createField(const std::string& name, size_t levels) const
{
    return createField(name,DataType::create<DATATYPE>(),levels);
}

template< typename DATATYPE >
Field* Edges::createField(const std::string& name,const std::vector<size_t>& variables) const
{
    return createField(name,DataType::create<DATATYPE>(),variables);
}

template< typename DATATYPE >
Field* Edges::createField(const std::string& name, size_t levels, const std::vector<size_t>& variables) const
{
    return createField(name,DataType::create<DATATYPE>(),levels,variables);
}

template< typename DATATYPE >
Field* Edges::createGlobalField(const std::string& name) const
{
    return createGlobalField(name,DataType::create<DATATYPE>());
}

template< typename DATATYPE >
Field* Edges::createGlobalField(const std::string& name,size_t levels) const
{
    return createGlobalField(name,DataType::create<DATATYPE>(),levels);
}

template< typename DATATYPE >
Field* Edges::createGlobalField(const std::string& name, const std::vector<size_t>& variables) const
{
    return createGlobalField(name,DataType::create<DATATYPE>(),variables);
}

template< typename DATATYPE >
Field* Edges::createGlobalField(const std::string& name, size_t levels, const std::vector<size_t>& variables) const
{
    return createGlobalField(name,DataType::create<DATATYPE>(),levels,variables);
}

// -------------------------------------------------------------------------------
#define mesh_Edges mesh::HybridElements
#define Char char
#define GatherScatter mpl::GatherScatter
#define Checksum mpl::Checksum
#define HaloExchange mpl::HaloExchange

extern "C" {

Edges* atlas__functionspace__Edges__new (Mesh* mesh, int halo);
Edges* atlas__functionspace__Edges__new_mesh (Mesh* mesh);
void atlas__functionspace__Edges__delete (Edges* This);
int atlas__functionspace__Edges__nb_edges(const Edges* This);
Mesh* atlas__functionspace__Edges__mesh(Edges* This);
mesh_Edges* atlas__functionspace__Edges__edges(Edges* This);
Field* atlas__functionspace__Edges__create_field (const Edges* This, const char* name, int kind);
Field* atlas__functionspace__Edges__create_field_vars (const Edges* This, const char* name, int variables[], int variables_size, int fortran_ordering, int kind);

Field* atlas__functionspace__Edges__create_field_lev (const Edges* This, const char* name, int levels, int kind);
Field* atlas__functionspace__Edges__create_field_lev_vars (const Edges* This, const char* name, int levels, int variables[], int variables_size, int fortran_ordering, int kind);


Field* atlas__functionspace__Edges__create_field_template (const Edges* This, const char* name, const Field* field_template);
Field* atlas__functionspace__Edges__create_global_field (const Edges* This, const char* name, int kind);
Field* atlas__functionspace__Edges__create_global_field_vars (const Edges* This, const char* name, int variables[], int variables_size, int fortran_ordering, int kind);

Field* atlas__functionspace__Edges__create_global_field_lev (const Edges* This, const char* name, int levels, int kind);
Field* atlas__functionspace__Edges__create_global_field_lev_vars (const Edges* This, const char* name, int levels, int variables[], int variables_size, int fortran_ordering, int kind);

Field* atlas__functionspace__Edges__create_global_field_template (const Edges* This, const char* name, const Field* field_template);

void atlas__functionspace__Edges__halo_exchange_fieldset(const Edges* This, FieldSet* fieldset);
void atlas__functionspace__Edges__halo_exchange_field(const Edges* This, Field* field);
const HaloExchange* atlas__functionspace__Edges__get_halo_exchange(const Edges* This);

void atlas__functionspace__Edges__gather_fieldset(const Edges* This, const FieldSet* local, FieldSet* global);
void atlas__functionspace__Edges__gather_field(const Edges* This, const Field* local, Field* global);
const GatherScatter* atlas__functionspace__Edges__get_gather(const Edges* This);

void atlas__functionspace__Edges__scatter_fieldset(const Edges* This, const FieldSet* global, FieldSet* local);
void atlas__functionspace__Edges__scatter_field(const Edges* This, const Field* global, Field* local);
const GatherScatter* atlas__functionspace__Edges__get_scatter(const Edges* This);

void atlas__functionspace__Edges__checksum_fieldset(const Edges* This, const FieldSet* fieldset, Char* &checksum, int &size, int &allocated);
void atlas__functionspace__Edges__checksum_field(const Edges* This, const Field* field, Char* &checksum, int &size, int &allocated);
const Checksum* atlas__functionspace__Edges__get_checksum(const Edges* This);
}

#undef mesh_Edges
#undef Char
#undef GatherScatter
#undef Checksum
#undef HaloExchange

} // namespace functionspace
} // namespace atlas

#endif // atlas_functionspace_EdgesFunctionSpace_h
