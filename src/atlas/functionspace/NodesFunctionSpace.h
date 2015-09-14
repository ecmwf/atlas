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
  Halo(const Mesh& mesh);
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
    std::vector<size_t> nb_nodes_global_foreach_rank() const;

    const Mesh& mesh() const { return mesh_; }
          Mesh& mesh()       { return mesh_; }

    const Nodes& nodes() const { return nodes_; }
          Nodes& nodes()       { return nodes_; }



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

    /// @brief Compute sum of scalar field
    /// @param [out] sum    Scalar value containing the sum of the full 3D field
    /// @param [out] N      Number of values that are contained in the sum (nodes*levels)
    template< typename DATATYPE >
    void sum( const Field&, DATATYPE& sum, size_t& N ) const;

    /// @brief Compute sum of field for each variable
    /// @param [out] sum    For each field-variable, the sum of the full 3D field
    /// @param [out] N      Number of values that are contained in the sum (nodes*levels)
    template< typename DATATYPE >
    void sum( const Field&, std::vector<DATATYPE>& sum, size_t& N ) const;

    /// @brief Compute sum of field for each vertical level separately
    /// @param [out] sum    Field of dimension of input without the nodes index
    /// @param [out] N      Number of nodes used to sum each level
    void sumPerLevel( const Field&, Field& sum, size_t& N ) const;

    /// @brief Compute order independent sum of scalar field
    /// @param [out] sum    Scalar value containing the sum of the full 3D field
    /// @param [out] N      Number of values that are contained in the sum (nodes*levels)
    template< typename DATATYPE >
    void orderIndependentSum( const Field&, DATATYPE& sum, size_t& N ) const;

    /// @brief Compute order independent sum of field for each variable
    /// @param [out] sum    For each field-variable, the sum of the full 3D field
    /// @param [out] N      Number of values that are contained in the sum (nodes*levels)
    template< typename DATATYPE >
    void orderIndependentSum( const Field&, std::vector<DATATYPE>&, size_t& N ) const;

    /// @brief Compute order independent sum of field for each vertical level separately
    /// @param [out] sum    Field of dimension of input without the nodes index
    /// @param [out] N      Number of nodes used to sum each level
    void orderIndependentSumPerLevel( const Field&, Field& sum, size_t& N ) const;

    /// @brief Compute minimum of scalar field
    template< typename DATATYPE >
    void minimum( const Field&, DATATYPE& minimum ) const;

    /// @brief Compute maximum of scalar field
    template< typename DATATYPE >
    void maximum( const Field&, DATATYPE& maximum ) const;

    /// @brief Compute minimum of field for each field-variable
    template< typename DATATYPE >
    void minimum( const Field&, std::vector<DATATYPE>& ) const;

    /// @brief Compute maximum of field for each field-variable
    template< typename DATATYPE >
    void maximum( const Field&, std::vector<DATATYPE>& ) const;

    /// @brief Compute minimum of field for each vertical level separately
    /// @param [out] min    Field of dimension of input without the nodes index
    void minimumPerLevel( const Field&, Field& min ) const;

    /// @brief Compute maximum of field for each vertical level separately
    /// @param [out] max    Field of dimension of input without the nodes index
    void maximumPerLevel( const Field&, Field& max ) const;

    /// @brief Compute minimum of scalar field, as well as the global index and level.
    template< typename DATATYPE >
    void minimumAndLocation( const Field&, DATATYPE& minimum, gidx_t& glb_idx ) const;

    /// @brief Compute maximum of scalar field, as well as the global index and level.
    template< typename DATATYPE >
    void maximumAndLocation( const Field&, DATATYPE& maximum, gidx_t& glb_idx ) const;

    /// @brief Compute minimum of scalar field, as well as the global index and level.
    template< typename DATATYPE >
    void minimumAndLocation( const Field&, DATATYPE& minimum, gidx_t& glb_idx, size_t& level ) const;

    /// @brief Compute maximum of scalar field, as well as the global index and level.
    template< typename DATATYPE >
    void maximumAndLocation( const Field&, DATATYPE& maximum, gidx_t& glb_idx, size_t& level ) const;

    /// @brief Compute minimum of field for each field-variable, as well as the global indices and levels.
    template< typename DATATYPE >
    void minimumAndLocation( const Field&, std::vector<DATATYPE>& minimum, std::vector<gidx_t>& glb_idx ) const;

    /// @brief Compute maximum of field for each field-variable, as well as the global indices and levels.
    template< typename DATATYPE >
    void maximumAndLocation( const Field&, std::vector<DATATYPE>& maximum, std::vector<gidx_t>& glb_idx ) const;

    /// @brief Compute minimum of field for each field-variable, as well as the global indices and levels.
    template< typename DATATYPE >
    void minimumAndLocation( const Field&, std::vector<DATATYPE>& minimum, std::vector<gidx_t>& glb_idx, std::vector<size_t>& level ) const;

    /// @brief Compute maximum of field for each field-variable, as well as the global indices and levels.
    template< typename DATATYPE >
    void maximumAndLocation( const Field&, std::vector<DATATYPE>& maximum, std::vector<gidx_t>& glb_idx, std::vector<size_t>& level ) const;

    /// @brief Compute minimum and its location of a field for each vertical level separately
    void minimumAndLocationPerLevel( const Field&, Field& column, Field& glb_idx ) const;

    /// @brief Compute maximum and its location of a field for each vertical level separately
    void maximumAndLocationPerLevel( const Field&, Field& column, Field& glb_idx ) const;

    /// @brief Compute mean value of scalar field
    /// @param [out] mean    Mean value
    /// @param [out] N       Number of value used to create the mean
    template< typename DATATYPE >
    void mean( const Field&, DATATYPE& mean, size_t& N ) const;

    /// @brief Compute mean value of field for each field-variable
    /// @param [out] mean    Mean values for each variable
    /// @param [out] N       Number of values used to create the means
    template< typename DATATYPE >
    void mean( const Field&, std::vector<DATATYPE>& mean, size_t& N ) const;

    /// @brief Compute mean values of field for vertical level separately
    /// @param [out] mean    Field of dimension of input without the nodes index
    /// @param [out] N       Number of values used to create the means
    void meanPerLevel( const Field&, Field& mean, size_t& N ) const;

    /// @brief Compute mean value and standard deviation of scalar field
    /// @param [out] mean      Mean value
    /// @param [out] stddev    Standard deviation
    /// @param [out] N         Number of value used to create the mean
    template< typename DATATYPE >
    void meanAndStandardDeviation( const Field&, DATATYPE& mean, DATATYPE& stddev, size_t& N ) const;

    /// @brief Compute mean values and standard deviations of scalar field for each field-variable
    /// @param [out] mean      Mean values for each field-variable
    /// @param [out] stddev    Standard deviation for each field-variable
    /// @param [out] N         Number of value used to create the means
    template< typename DATATYPE >
    void meanAndstandardDeviation( const Field&, std::vector<DATATYPE>& mean, std::vector<DATATYPE>& stddev, size_t& N ) const;

    /// @brief Compute mean values and standard deviations of field for vertical level separately
    /// @param [out] mean      Field of dimension of input without the nodes index
    /// @param [out] stddev    Field of dimension of input without the nodes index
    /// @param [out] N         Number of values used to create the means
    void meanAndStandardDeviationPerLevel( const Field&, Field& mean, Field& stddev, size_t& N ) const;

private: // methods

    std::string halo_name() const;
    std::string gather_scatter_name() const;
    std::string checksum_name() const;

private: // data

    Mesh& mesh_; // non-const because functionspace may modify mesh
    Nodes& nodes_; // non-const because functionspace may modify mesh
    Halo halo_;
    size_t nb_nodes_;
    size_t nb_nodes_global_;
    std::vector<size_t> nb_nodes_global_foreach_rank_;
};

// -------------------------------------------------------------------

typedef NodesFunctionSpace NodesColumnFunctionSpace;

// -------------------------------------------------------------------

template< typename DATATYPE >
Field* NodesFunctionSpace::createField(const std::string& name) const
{
    return createField(name,DataType::create<DATATYPE>());
}

template< typename DATATYPE >
Field* NodesFunctionSpace::createField(const std::string& name, size_t levels) const
{
    return createField(name,DataType::create<DATATYPE>(),levels);
}

template< typename DATATYPE >
Field* NodesFunctionSpace::createField(const std::string& name,const std::vector<size_t>& variables) const
{
    return createField(name,DataType::create<DATATYPE>(),variables);
}

template< typename DATATYPE >
Field* NodesFunctionSpace::createField(const std::string& name, size_t levels, const std::vector<size_t>& variables) const
{
    return createField(name,DataType::create<DATATYPE>(),levels,variables);
}

template< typename DATATYPE >
Field* NodesFunctionSpace::createGlobalField(const std::string& name) const
{
    return createGlobalField(name,DataType::create<DATATYPE>());
}

template< typename DATATYPE >
Field* NodesFunctionSpace::createGlobalField(const std::string& name,size_t levels) const
{
    return createGlobalField(name,DataType::create<DATATYPE>(),levels);
}

template< typename DATATYPE >
Field* NodesFunctionSpace::createGlobalField(const std::string& name, const std::vector<size_t>& variables) const
{
    return createGlobalField(name,DataType::create<DATATYPE>(),variables);
}

template< typename DATATYPE >
Field* NodesFunctionSpace::createGlobalField(const std::string& name, size_t levels, const std::vector<size_t>& variables) const
{
    return createGlobalField(name,DataType::create<DATATYPE>(),levels,variables);
}

} // namespace functionspace
} // namespace atlas

#endif // atlas_functionspace_NodesFunctionSpace_h
