/*
 * (C) Copyright 1996-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef atlas_functionspace_NodeColumnsFunctionSpace_h
#define atlas_functionspace_NodeColumnsFunctionSpace_h

#include "eckit/memory/SharedPtr.h"
#include "atlas/internals/atlas_config.h"
#include "atlas/mesh/Halo.h"
#include "atlas/field/FieldSet.h"
#include "atlas/field/Options.h"
#include "atlas/functionspace/FunctionSpace.h"

// ----------------------------------------------------------------------------
// Forward declarations

namespace atlas {
namespace mesh {
    class Mesh;
    class Nodes;
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

class NodeColumns : public FunctionSpace
{
public:

    NodeColumns( mesh::Mesh& mesh, const mesh::Halo &, const eckit::Parametrisation & );
    NodeColumns( mesh::Mesh& mesh, const mesh::Halo & );
    NodeColumns( mesh::Mesh& mesh );

    virtual ~NodeColumns();

    virtual std::string name() const { return "Nodes"; }

    size_t nb_nodes() const;
    size_t nb_nodes_global() const; // All MPI ranks will have same output
    std::vector<size_t> nb_nodes_global_foreach_rank() const;

    const mesh::Mesh& mesh() const { return mesh_; }
          mesh::Mesh& mesh()       { return mesh_; }

    const mesh::Nodes& nodes() const { return nodes_; }
          mesh::Nodes& nodes()       { return nodes_; }



// -- Local Field creation methods

    /// @brief Create a named scalar field
    template< typename DATATYPE > field::Field* createField(
              const std::string& name,
              const eckit::Parametrisation& = util::NoConfig()) const;

    template< typename DATATYPE > field::Field* createField(
              const std::string& name,
              size_t levels,
              const eckit::Parametrisation& = util::NoConfig()) const;

    /// @brief Create a named scalar field
    field::Field* createField(
        const std::string& name,
        array::DataType,
        const eckit::Parametrisation& = util::NoConfig() ) const;

    field::Field* createField(
        const std::string& name,
        array::DataType, size_t levels,
        const eckit::Parametrisation& = util::NoConfig()) const;

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
    field::Field* createField(
        const std::string& name,
        array::DataType,
        const std::vector<size_t>& variables,
        const eckit::Parametrisation& = util::NoConfig()) const;

    field::Field* createField(
        const std::string& name,
        array::DataType, size_t levels,
        const std::vector<size_t>& variables,
        const eckit::Parametrisation& = util::NoConfig())  const;

    /// @brief Create a named field based on other field (datatype and dimensioning)
    field::Field* createField(
        const std::string& name,
        const field::Field&,
        const eckit::Parametrisation& = util::NoConfig()) const;

    field::Field* createField(const eckit::Parametrisation&) const;

// -- Global field::Field creation methods

    /// @brief Create a named global scalar field
    template< typename DATATYPE >  field::Field* createGlobalField(const std::string& name) const;
    template< typename DATATYPE >  field::Field* createGlobalField(const std::string& name,size_t levels) const;

    /// @brief Create a named global scalar field
    field::Field* createGlobalField(const std::string& name, array::DataType) const;
    field::Field* createGlobalField(const std::string& name, array::DataType, size_t levels) const;

    /// @brief Create a named global field with specified dimensions for the variables
    template< typename DATATYPE >  field::Field* createGlobalField(const std::string& name, const std::vector<size_t>& variables) const;
    template< typename DATATYPE >  field::Field* createGlobalField(const std::string& name, size_t levels, const std::vector<size_t>& variables) const;

    /// @brief Create a named field with specified dimensions for the variables
    field::Field* createGlobalField(const std::string& name, array::DataType, const std::vector<size_t>& variables) const;
    field::Field* createGlobalField(const std::string& name, array::DataType, size_t levels, const std::vector<size_t>& variables) const;

    /// @brief Create a named global field based on other field (datatype and dimensioning)
    field::Field* createGlobalField(const std::string& name, const field::Field&) const;

// -- Parallelisation aware methods

    const mesh::Halo& halo() const { return halo_; }

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

    /// @brief Compute sum of scalar field
    /// @param [out] sum    Scalar value containing the sum of the full 3D field
    /// @param [out] N      Number of values that are contained in the sum (nodes*levels)
    template< typename DATATYPE >
    void sum( const field::Field&, DATATYPE& sum, size_t& N ) const;

    /// @brief Compute sum of field for each variable
    /// @param [out] sum    For each field-variable, the sum of the full 3D field
    /// @param [out] N      Number of values that are contained in the sum (nodes*levels)
    template< typename DATATYPE >
    void sum( const field::Field&, std::vector<DATATYPE>& sum, size_t& N ) const;

    /// @brief Compute sum of field for each vertical level separately
    /// @param [out] sum    field::Field of dimension of input without the nodes index
    /// @param [out] N      Number of nodes used to sum each level
    void sumPerLevel( const field::Field&, field::Field& sum, size_t& N ) const;

    /// @brief Compute order independent sum of scalar field
    /// @param [out] sum    Scalar value containing the sum of the full 3D field
    /// @param [out] N      Number of values that are contained in the sum (nodes*levels)
    template< typename DATATYPE >
    void orderIndependentSum( const field::Field&, DATATYPE& sum, size_t& N ) const;

    /// @brief Compute order independent sum of field for each variable
    /// @param [out] sum    For each field-variable, the sum of the full 3D field
    /// @param [out] N      Number of values that are contained in the sum (nodes*levels)
    template< typename DATATYPE >
    void orderIndependentSum( const field::Field&, std::vector<DATATYPE>&, size_t& N ) const;

    /// @brief Compute order independent sum of field for each vertical level separately
    /// @param [out] sum    field::Field of dimension of input without the nodes index
    /// @param [out] N      Number of nodes used to sum each level
    void orderIndependentSumPerLevel( const field::Field&, field::Field& sum, size_t& N ) const;

    /// @brief Compute minimum of scalar field
    template< typename DATATYPE >
    void minimum( const field::Field&, DATATYPE& minimum ) const;

    /// @brief Compute maximum of scalar field
    template< typename DATATYPE >
    void maximum( const field::Field&, DATATYPE& maximum ) const;

    /// @brief Compute minimum of field for each field-variable
    template< typename DATATYPE >
    void minimum( const field::Field&, std::vector<DATATYPE>& ) const;

    /// @brief Compute maximum of field for each field-variable
    template< typename DATATYPE >
    void maximum( const field::Field&, std::vector<DATATYPE>& ) const;

    /// @brief Compute minimum of field for each vertical level separately
    /// @param [out] min    field::Field of dimension of input without the nodes index
    void minimumPerLevel( const field::Field&, field::Field& min ) const;

    /// @brief Compute maximum of field for each vertical level separately
    /// @param [out] max    field::Field of dimension of input without the nodes index
    void maximumPerLevel( const field::Field&, field::Field& max ) const;

    /// @brief Compute minimum of scalar field, as well as the global index and level.
    template< typename DATATYPE >
    void minimumAndLocation( const field::Field&, DATATYPE& minimum, gidx_t& glb_idx ) const;

    /// @brief Compute maximum of scalar field, as well as the global index and level.
    template< typename DATATYPE >
    void maximumAndLocation( const field::Field&, DATATYPE& maximum, gidx_t& glb_idx ) const;

    /// @brief Compute minimum of scalar field, as well as the global index and level.
    template< typename DATATYPE >
    void minimumAndLocation( const field::Field&, DATATYPE& minimum, gidx_t& glb_idx, size_t& level ) const;

    /// @brief Compute maximum of scalar field, as well as the global index and level.
    template< typename DATATYPE >
    void maximumAndLocation( const field::Field&, DATATYPE& maximum, gidx_t& glb_idx, size_t& level ) const;

    /// @brief Compute minimum of field for each field-variable, as well as the global indices and levels.
    template< typename DATATYPE >
    void minimumAndLocation( const field::Field&, std::vector<DATATYPE>& minimum, std::vector<gidx_t>& glb_idx ) const;

    /// @brief Compute maximum of field for each field-variable, as well as the global indices and levels.
    template< typename DATATYPE >
    void maximumAndLocation( const field::Field&, std::vector<DATATYPE>& maximum, std::vector<gidx_t>& glb_idx ) const;

    /// @brief Compute minimum of field for each field-variable, as well as the global indices and levels.
    template< typename DATATYPE >
    void minimumAndLocation( const field::Field&, std::vector<DATATYPE>& minimum, std::vector<gidx_t>& glb_idx, std::vector<size_t>& level ) const;

    /// @brief Compute maximum of field for each field-variable, as well as the global indices and levels.
    template< typename DATATYPE >
    void maximumAndLocation( const field::Field&, std::vector<DATATYPE>& maximum, std::vector<gidx_t>& glb_idx, std::vector<size_t>& level ) const;

    /// @brief Compute minimum and its location of a field for each vertical level separately
    void minimumAndLocationPerLevel( const field::Field&, field::Field& column, field::Field& glb_idx ) const;

    /// @brief Compute maximum and its location of a field for each vertical level separately
    void maximumAndLocationPerLevel( const field::Field&, field::Field& column, field::Field& glb_idx ) const;

    /// @brief Compute mean value of scalar field
    /// @param [out] mean    Mean value
    /// @param [out] N       Number of value used to create the mean
    template< typename DATATYPE >
    void mean( const field::Field&, DATATYPE& mean, size_t& N ) const;

    /// @brief Compute mean value of field for each field-variable
    /// @param [out] mean    Mean values for each variable
    /// @param [out] N       Number of values used to create the means
    template< typename DATATYPE >
    void mean( const field::Field&, std::vector<DATATYPE>& mean, size_t& N ) const;

    /// @brief Compute mean values of field for vertical level separately
    /// @param [out] mean    field::Field of dimension of input without the nodes index
    /// @param [out] N       Number of values used to create the means
    void meanPerLevel( const field::Field&, field::Field& mean, size_t& N ) const;

    /// @brief Compute mean value and standard deviation of scalar field
    /// @param [out] mean      Mean value
    /// @param [out] stddev    Standard deviation
    /// @param [out] N         Number of value used to create the mean
    template< typename DATATYPE >
    void meanAndStandardDeviation( const field::Field&, DATATYPE& mean, DATATYPE& stddev, size_t& N ) const;

    /// @brief Compute mean values and standard deviations of scalar field for each field-variable
    /// @param [out] mean      Mean values for each field-variable
    /// @param [out] stddev    Standard deviation for each field-variable
    /// @param [out] N         Number of value used to create the means
    template< typename DATATYPE >
    void meanAndstandardDeviation( const field::Field&, std::vector<DATATYPE>& mean, std::vector<DATATYPE>& stddev, size_t& N ) const;

    /// @brief Compute mean values and standard deviations of field for vertical level separately
    /// @param [out] mean      field::Field of dimension of input without the nodes index
    /// @param [out] stddev    field::Field of dimension of input without the nodes index
    /// @param [out] N         Number of values used to create the means
    void meanAndStandardDeviationPerLevel( const field::Field&, field::Field& mean, field::Field& stddev, size_t& N ) const;

private: // methods

    std::string halo_name() const;
    std::string gather_scatter_name() const;
    std::string checksum_name() const;
    void constructor();

    size_t config_nb_nodes(const eckit::Parametrisation&) const;
    array::DataType config_datatype(const eckit::Parametrisation&) const;
    std::string config_name(const eckit::Parametrisation&) const;
    size_t config_levels(const eckit::Parametrisation&) const;

private: // data

    mesh::Mesh& mesh_; // non-const because functionspace may modify mesh
    mesh::Nodes& nodes_; // non-const because functionspace may modify mesh
    mesh::Halo halo_;
    size_t nb_nodes_;
    size_t nb_nodes_global_;
    std::vector<size_t> nb_nodes_global_foreach_rank_;

    eckit::SharedPtr<parallel::GatherScatter> gather_scatter_; // without ghost
    eckit::SharedPtr<parallel::HaloExchange>  halo_exchange_;
    eckit::SharedPtr<parallel::Checksum>      checksum_;

};

// -------------------------------------------------------------------

template< typename DATATYPE >
field::Field* NodeColumns::createField(
    const std::string& name,
    const eckit::Parametrisation& options) const
{
    return createField(name,array::DataType::create<DATATYPE>(),options);
}

template< typename DATATYPE >
field::Field* NodeColumns::createField(
    const std::string& name,
    size_t levels,
    const eckit::Parametrisation& options) const
{
    return createField(name,array::DataType::create<DATATYPE>(),levels,options);
}

template< typename DATATYPE >
field::Field* NodeColumns::createField(
    const std::string& name,
    const std::vector<size_t>& variables,
    const eckit::Parametrisation& options) const
{
    return createField(name,array::DataType::create<DATATYPE>(),variables,options);
}

template< typename DATATYPE >
field::Field* NodeColumns::createField(
    const std::string& name,
    size_t levels,
    const std::vector<size_t>& variables,
    const eckit::Parametrisation& options) const
{
    return createField(name,array::DataType::create<DATATYPE>(),levels,variables,options);
}

template< typename DATATYPE >
field::Field* NodeColumns::createGlobalField(const std::string& name) const
{
    return createGlobalField(name,array::DataType::create<DATATYPE>());
}

template< typename DATATYPE >
field::Field* NodeColumns::createGlobalField(const std::string& name,size_t levels) const
{
    return createGlobalField(name,array::DataType::create<DATATYPE>(),levels);
}

template< typename DATATYPE >
field::Field* NodeColumns::createGlobalField(const std::string& name, const std::vector<size_t>& variables) const
{
    return createGlobalField(name,array::DataType::create<DATATYPE>(),variables);
}

template< typename DATATYPE >
field::Field* NodeColumns::createGlobalField(const std::string& name, size_t levels, const std::vector<size_t>& variables) const
{
    return createGlobalField(name,array::DataType::create<DATATYPE>(),levels,variables);
}

} // namespace functionspace
} // namespace atlas

#endif // atlas_functionspace_NodeColumnsFunctionSpace_h
