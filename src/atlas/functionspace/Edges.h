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

    const atlas::Mesh& mesh() const { return mesh_; }
          atlas::Mesh& mesh()       { return mesh_; }

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

    Mesh& mesh_; // non-const because functionspace may modify mesh
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

} // namespace functionspace
} // namespace atlas

#endif // atlas_functionspace_EdgesFunctionSpace_h
