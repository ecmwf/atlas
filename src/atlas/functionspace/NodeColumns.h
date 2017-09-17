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
#include "atlas/library/config.h"
#include "atlas/mesh/Halo.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/field/FieldSet.h"
#include "atlas/option.h"
#include "atlas/functionspace/FunctionSpace.h"

// ----------------------------------------------------------------------------
// Forward declarations

namespace atlas {
namespace mesh {
    class Nodes;
}
}

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

class NodeColumns : public FunctionSpaceImpl
{
public:

    NodeColumns( Mesh& mesh, const mesh::Halo & );
    NodeColumns( Mesh& mesh, const eckit::Configuration & );
    NodeColumns( Mesh& mesh );

    virtual ~NodeColumns();

    virtual std::string name() const { return "Nodes"; }

    size_t nb_nodes() const;
    size_t nb_nodes_global() const; // All MPI ranks will have same output

    const Mesh& mesh() const { return mesh_; }
    
    size_t levels() const { return nb_levels_; }

    mesh::Nodes& nodes() const { return nodes_; }

// -- Field creation methods
    
    using FunctionSpaceImpl::createField;

    /// @brief Create a field
    virtual Field createField(const eckit::Configuration&) const;

    virtual Field createField(const Field&, const eckit::Configuration& ) const;

// -- Parallelisation aware methods

    const mesh::Halo& halo() const { return halo_; }

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


    /// @brief Compute sum of scalar field
    /// @param [out] sum    Scalar value containing the sum of the full 3D field
    /// @param [out] N      Number of values that are contained in the sum (nodes*levels)
    template< typename Value >
    void sum( const Field&, Value& sum, size_t& N ) const;

//    /// @brief Compute sum of field for each variable
//    /// @param [out] sum    For each field-variable, the sum of the full 3D field
//    /// @param [out] N      Number of values that are contained in the sum (nodes*levels)
//    template< typename Value >
//    void sum( const Field&, std::vector<Value>& sum, size_t& N ) const;

    /// @brief Compute sum of field for each vertical level separately
    /// @param [out] sum    Field of dimension of input without the nodes index
    /// @param [out] N      Number of nodes used to sum each level
    void sumPerLevel( const Field&, Field& sum, size_t& N ) const;

    /// @brief Compute order independent sum of scalar field
    /// @param [out] sum    Scalar value containing the sum of the full 3D field
    /// @param [out] N      Number of values that are contained in the sum (nodes*levels)
    template< typename Value >
    void orderIndependentSum( const Field&, Value& sum, size_t& N ) const;

//    /// @brief Compute order independent sum of field for each variable
//    /// @param [out] sum    For each field-variable, the sum of the full 3D field
//    /// @param [out] N      Number of values that are contained in the sum (nodes*levels)
//    template< typename Value >
//    void orderIndependentSum( const Field&, std::vector<Value>&, size_t& N ) const;

    /// @brief Compute order independent sum of field for each vertical level separately
    /// @param [out] sum    Field of dimension of input without the nodes index
    /// @param [out] N      Number of nodes used to sum each level
    void orderIndependentSumPerLevel( const Field&, Field& sum, size_t& N ) const;

    /// @brief Compute minimum of scalar field
    template< typename Value >
    void minimum( const Field&, Value& minimum ) const;

    /// @brief Compute maximum of scalar field
    template< typename Value >
    void maximum( const Field&, Value& maximum ) const;

//    /// @brief Compute minimum of field for each field-variable
//    template< typename Value >
//    void minimum( const Field&, std::vector<Value>& ) const;

//    /// @brief Compute maximum of field for each field-variable
//    template< typename Value >
//    void maximum( const Field&, std::vector<Value>& ) const;

    /// @brief Compute minimum of field for each vertical level separately
    /// @param [out] min    Field of dimension of input without the nodes index
    void minimumPerLevel( const Field&, Field& min ) const;

    /// @brief Compute maximum of field for each vertical level separately
    /// @param [out] max    Field of dimension of input without the nodes index
    void maximumPerLevel( const Field&, Field& max ) const;

    /// @brief Compute minimum of scalar field, as well as the global index and level.
    template< typename Value >
    void minimumAndLocation( const Field&, Value& minimum, gidx_t& glb_idx ) const;

    /// @brief Compute maximum of scalar field, as well as the global index and level.
    template< typename Value >
    void maximumAndLocation( const Field&, Value& maximum, gidx_t& glb_idx ) const;

    /// @brief Compute minimum of scalar field, as well as the global index and level.
    template< typename Value >
    void minimumAndLocation( const Field&, Value& minimum, gidx_t& glb_idx, size_t& level ) const;

    /// @brief Compute maximum of scalar field, as well as the global index and level.
    template< typename Value >
    void maximumAndLocation( const Field&, Value& maximum, gidx_t& glb_idx, size_t& level ) const;

    /// @brief Compute minimum of field for each field-variable, as well as the global indices and levels.
    template< typename Vector >
    void minimumAndLocation( const Field&, Vector& minimum, std::vector<gidx_t>& glb_idx ) const;

    /// @brief Compute maximum of field for each field-variable, as well as the global indices and levels.
    template< typename Vector >
    void maximumAndLocation( const Field&, Vector& maximum, std::vector<gidx_t>& glb_idx ) const;

    /// @brief Compute minimum of field for each field-variable, as well as the global indices and levels.
    template< typename Vector >
    void minimumAndLocation( const Field&, Vector& minimum, std::vector<gidx_t>& glb_idx, std::vector<size_t>& level ) const;

    /// @brief Compute maximum of field for each field-variable, as well as the global indices and levels.
    template< typename Vector >
    void maximumAndLocation( const Field&, Vector& maximum, std::vector<gidx_t>& glb_idx, std::vector<size_t>& level ) const;

    /// @brief Compute minimum and its location of a field for each vertical level separately
    void minimumAndLocationPerLevel( const Field&, Field& column, Field& glb_idx ) const;

    /// @brief Compute maximum and its location of a field for each vertical level separately
    void maximumAndLocationPerLevel( const Field&, Field& column, Field& glb_idx ) const;

    /// @brief Compute mean value of scalar field
    /// @param [out] mean    Mean value
    /// @param [out] N       Number of value used to create the mean
    template< typename Value >
    void mean( const Field&, Value& mean, size_t& N ) const;

//    /// @brief Compute mean value of field for each field-variable
//    /// @param [out] mean    Mean values for each variable
//    /// @param [out] N       Number of values used to create the means
//    template< typename Value >
//    void mean( const Field&, std::vector<Value>& mean, size_t& N ) const;

    /// @brief Compute mean values of field for vertical level separately
    /// @param [out] mean    Field of dimension of input without the nodes index
    /// @param [out] N       Number of values used to create the means
    void meanPerLevel( const Field&, Field& mean, size_t& N ) const;

    /// @brief Compute mean value and standard deviation of scalar field
    /// @param [out] mean      Mean value
    /// @param [out] stddev    Standard deviation
    /// @param [out] N         Number of value used to create the mean
    template< typename Value >
    void meanAndStandardDeviation( const Field&, Value& mean, Value& stddev, size_t& N ) const;

//    /// @brief Compute mean values and standard deviations of scalar field for each field-variable
//    /// @param [out] mean      Mean values for each field-variable
//    /// @param [out] stddev    Standard deviation for each field-variable
//    /// @param [out] N         Number of value used to create the means
//    template< typename Value >
//    void meanAndStandardDeviation( const Field&, std::vector<Value>& mean, std::vector<Value>& stddev, size_t& N ) const;

    /// @brief Compute mean values and standard deviations of field for vertical level separately
    /// @param [out] mean      Field of dimension of input without the nodes index
    /// @param [out] stddev    Field of dimension of input without the nodes index
    /// @param [out] N         Number of values used to create the means
    void meanAndStandardDeviationPerLevel( const Field&, Field& mean, Field& stddev, size_t& N ) const;

private: // methods

    std::string halo_name() const;
    std::string gather_scatter_name() const;
    std::string checksum_name() const;
    void constructor();

    size_t config_nb_nodes(const eckit::Configuration&) const;
    array::DataType config_datatype(const eckit::Configuration&) const;
    std::string config_name(const eckit::Configuration&) const;
    size_t config_levels(const eckit::Configuration&) const;
    array::ArrayShape config_shape(const eckit::Configuration&) const;
    void set_field_metadata(const eckit::Configuration&, Field& ) const;
    
    size_t footprint() const;

private: // data

    Mesh mesh_;    // non-const because functionspace may modify mesh
    mesh::Nodes& nodes_; // non-const because functionspace may modify mesh
    mesh::Halo halo_;
    size_t nb_nodes_;
    size_t nb_nodes_global_;
    size_t nb_levels_;

    eckit::SharedPtr<parallel::GatherScatter> gather_scatter_; // without ghost
    eckit::SharedPtr<parallel::HaloExchange>  halo_exchange_;
    eckit::SharedPtr<parallel::Checksum>      checksum_;

private:

    template< typename Value >
    struct FieldStatisticsT {
      FieldStatisticsT(const NodeColumns*);
      void sum( const Field&, Value& sum, size_t& N ) const;
      void orderIndependentSum( const Field&, Value& sum, size_t& N ) const;
      void minimum( const Field&, Value& minimum ) const;
      void maximum( const Field&, Value& maximum ) const;
      void minimumAndLocation( const Field&, Value& minimum, gidx_t& glb_idx ) const;
      void maximumAndLocation( const Field&, Value& maximum, gidx_t& glb_idx ) const;
      void minimumAndLocation( const Field&, Value& minimum, gidx_t& glb_idx, size_t& level ) const;
      void maximumAndLocation( const Field&, Value& maximum, gidx_t& glb_idx, size_t& level ) const;
      void mean( const Field&, Value& mean, size_t& N ) const;
      void meanAndStandardDeviation( const Field&, Value& mean, Value& stddev, size_t& N ) const;
      const NodeColumns& functionspace;
    };

    template< typename Vector >
    struct FieldStatisticsVectorT {
      FieldStatisticsVectorT(const NodeColumns*);
      void sum( const Field&, Vector& sum, size_t& N ) const;
      void orderIndependentSum( const Field&, Vector&, size_t& N ) const;
      void minimum( const Field&, Vector& ) const;
      void maximum( const Field&, Vector& ) const;
      void minimumAndLocation( const Field&, Vector& minimum, std::vector<gidx_t>& glb_idx ) const;
      void maximumAndLocation( const Field&, Vector& maximum, std::vector<gidx_t>& glb_idx ) const;
      void minimumAndLocation( const Field&, Vector& minimum, std::vector<gidx_t>& glb_idx, std::vector<size_t>& level ) const;
      void maximumAndLocation( const Field&, Vector& maximum, std::vector<gidx_t>& glb_idx, std::vector<size_t>& level ) const;
      void mean( const Field&, Vector& mean, size_t& N ) const;
      void meanAndStandardDeviation( const Field&, Vector& mean, Vector& stddev, size_t& N ) const;
      const NodeColumns& functionspace;
    };

    struct FieldStatistics {
      FieldStatistics(const NodeColumns*);
      void sumPerLevel( const Field&, Field& sum, size_t& N ) const;
      void orderIndependentSumPerLevel( const Field&, Field& sum, size_t& N ) const;
      void minimumPerLevel( const Field&, Field& min ) const;
      void maximumPerLevel( const Field&, Field& max ) const;
      void minimumAndLocationPerLevel( const Field&, Field& column, Field& glb_idx ) const;
      void maximumAndLocationPerLevel( const Field&, Field& column, Field& glb_idx ) const;
      void meanPerLevel( const Field&, Field& mean, size_t& N ) const;
      void meanAndStandardDeviationPerLevel( const Field&, Field& mean, Field& stddev, size_t& N ) const;
      const NodeColumns& functionspace;
    };

    template< typename T >
    struct FieldStatisticsSelector {
      using type = typename std::conditional< std::is_pod<T>::value,
          FieldStatisticsT<T>,
          FieldStatisticsVectorT<T>
      >::type;
    };
};

// -------------------------------------------------------------------

template< typename Value >
void NodeColumns::sum( const Field& field, Value& sum, size_t& N ) const {
  typename FieldStatisticsSelector<Value>::type(this).sum(field,sum,N);
}

inline void NodeColumns::sumPerLevel( const Field& field, Field& sum, size_t& N ) const {
  FieldStatistics(this).sumPerLevel(field,sum,N);
}

template< typename Value >
void NodeColumns::orderIndependentSum( const Field& field, Value& sum, size_t& N ) const {
  typename FieldStatisticsSelector<Value>::type(this).orderIndependentSum(field,sum,N);
}

inline void NodeColumns::orderIndependentSumPerLevel( const Field& field, Field& sum, size_t& N ) const {
  FieldStatistics(this).orderIndependentSumPerLevel(field,sum,N);
}

template< typename Value >
void NodeColumns::minimum( const Field& field, Value& minimum ) const {
  typename FieldStatisticsSelector<Value>::type(this).minimum(field,minimum);
}

template< typename Value >
void NodeColumns::maximum( const Field& field, Value& maximum ) const {
  typename FieldStatisticsSelector<Value>::type(this).maximum(field,maximum);
}

inline void NodeColumns::minimumPerLevel( const Field& field, Field& minimum ) const {
  return FieldStatistics(this).minimumPerLevel(field,minimum);
}

inline void NodeColumns::maximumPerLevel( const Field& field, Field& maximum ) const {
  FieldStatistics(this).maximumPerLevel(field,maximum);
}

template< typename Value >
void NodeColumns::minimumAndLocation( const Field& field, Value& minimum, gidx_t& glb_idx ) const {
  FieldStatisticsT<Value>(this).minimumAndLocation(field,minimum,glb_idx);
}

template< typename Value >
void NodeColumns::maximumAndLocation( const Field& field, Value& maximum, gidx_t& glb_idx ) const {
  FieldStatisticsT<Value>(this).maximumAndLocation(field,maximum,glb_idx);
}

template< typename Value >
void NodeColumns::minimumAndLocation( const Field& field, Value& minimum, gidx_t& glb_idx, size_t& level ) const {
  FieldStatisticsT<Value>(this).minimumAndLocation(field,minimum,glb_idx,level);
}

template< typename Value >
void NodeColumns::maximumAndLocation( const Field& field, Value& maximum, gidx_t& glb_idx, size_t& level ) const {
  FieldStatisticsT<Value>(this).maximumAndLocation(field,maximum,glb_idx,level);
}

template< typename Vector >
void NodeColumns::minimumAndLocation( const Field& field, Vector& minimum, std::vector<gidx_t>& glb_idx ) const {
  FieldStatisticsVectorT<Vector>(this).minimumAndLocation(field,minimum,glb_idx);
}

template< typename Vector >
void NodeColumns::maximumAndLocation( const Field& field, Vector& maximum, std::vector<gidx_t>& glb_idx ) const {
  FieldStatisticsVectorT<Vector>(this).maximumAndLocation(field,maximum,glb_idx);
}

template< typename Vector >
void NodeColumns::minimumAndLocation( const Field& field, Vector& minimum, std::vector<gidx_t>& glb_idx, std::vector<size_t>& level ) const {
  FieldStatisticsVectorT<Vector>(this).minimumAndLocation(field,minimum,glb_idx,level);
}

template< typename Vector >
void NodeColumns::maximumAndLocation( const Field& field, Vector& maximum, std::vector<gidx_t>& glb_idx, std::vector<size_t>& level ) const {
  FieldStatisticsVectorT<Vector>(this).maximumAndLocation(field,maximum,glb_idx,level);
}

inline void NodeColumns::minimumAndLocationPerLevel( const Field& field, Field& column, Field& glb_idx ) const {
  FieldStatistics(this).minimumAndLocationPerLevel(field,column,glb_idx);
}

inline void NodeColumns::maximumAndLocationPerLevel( const Field& field, Field& column, Field& glb_idx ) const {
  FieldStatistics(this).maximumAndLocationPerLevel(field,column,glb_idx);
}

template< typename Value >
void NodeColumns::mean( const Field& field, Value& mean, size_t& N ) const {
  typename FieldStatisticsSelector<Value>::type(this).mean(field,mean,N);
}

inline void NodeColumns::meanPerLevel( const Field& field, Field& mean, size_t& N ) const {
  FieldStatistics(this).meanPerLevel(field,mean,N);
}

template< typename Value >
void NodeColumns::meanAndStandardDeviation( const Field& field, Value& mean, Value& stddev, size_t& N ) const {
  typename FieldStatisticsSelector<Value>::type(this).meanAndStandardDeviation(field,mean,stddev,N);
}

inline void NodeColumns::meanAndStandardDeviationPerLevel( const Field& field, Field& mean, Field& stddev, size_t& N ) const {
  FieldStatistics(this).meanAndStandardDeviationPerLevel(field,mean,stddev,N);
}

} // namespace detail

// -------------------------------------------------------------------

class NodeColumns : public FunctionSpace
{
public:


    NodeColumns();
    NodeColumns( const FunctionSpace& );

    NodeColumns( Mesh& mesh, const mesh::Halo & );
    NodeColumns( Mesh& mesh );
    NodeColumns( Mesh& mesh, const eckit::Configuration & );

    operator bool() const { return valid(); }
    bool valid() const { return functionspace_; }

    size_t nb_nodes() const;
    size_t nb_nodes_global() const; // All MPI ranks will have same output

    const Mesh& mesh() const;

    size_t levels() const;

    mesh::Nodes& nodes() const;

// -- Parallelisation aware methods

    const mesh::Halo& halo() const;

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


    /// @brief Compute sum of scalar field
    /// @param [out] sum    Scalar value containing the sum of the full 3D field
    /// @param [out] N      Number of values that are contained in the sum (nodes*levels)
    template< typename Value >
    void sum( const Field&, Value& sum, size_t& N ) const;

//    /// @brief Compute sum of field for each variable
//    /// @param [out] sum    For each field-variable, the sum of the full 3D field
//    /// @param [out] N      Number of values that are contained in the sum (nodes*levels)
//    template< typename Value >
//    void sum( const Field&, std::vector<Value>& sum, size_t& N ) const;

    /// @brief Compute sum of field for each vertical level separately
    /// @param [out] sum    Field of dimension of input without the nodes index
    /// @param [out] N      Number of nodes used to sum each level
    void sumPerLevel( const Field&, Field& sum, size_t& N ) const;

    /// @brief Compute order independent sum of scalar field
    /// @param [out] sum    Scalar value containing the sum of the full 3D field
    /// @param [out] N      Number of values that are contained in the sum (nodes*levels)
    template< typename Value >
    void orderIndependentSum( const Field&, Value& sum, size_t& N ) const;

//    /// @brief Compute order independent sum of field for each variable
//    /// @param [out] sum    For each field-variable, the sum of the full 3D field
//    /// @param [out] N      Number of values that are contained in the sum (nodes*levels)
//    template< typename Value >
//    void orderIndependentSum( const Field&, std::vector<Value>&, size_t& N ) const;

    /// @brief Compute order independent sum of field for each vertical level separately
    /// @param [out] sum    Field of dimension of input without the nodes index
    /// @param [out] N      Number of nodes used to sum each level
    void orderIndependentSumPerLevel( const Field&, Field& sum, size_t& N ) const;

    /// @brief Compute minimum of scalar field
    template< typename Value >
    void minimum( const Field&, Value& minimum ) const;

    /// @brief Compute maximum of scalar field
    template< typename Value >
    void maximum( const Field&, Value& maximum ) const;

//    /// @brief Compute minimum of field for each field-variable
//    template< typename Value >
//    void minimum( const Field&, std::vector<Value>& ) const;

//    /// @brief Compute maximum of field for each field-variable
//    template< typename Value >
//    void maximum( const Field&, std::vector<Value>& ) const;

    /// @brief Compute minimum of field for each vertical level separately
    /// @param [out] min    Field of dimension of input without the nodes index
    void minimumPerLevel( const Field&, Field& min ) const;

    /// @brief Compute maximum of field for each vertical level separately
    /// @param [out] max    Field of dimension of input without the nodes index
    void maximumPerLevel( const Field&, Field& max ) const;

    /// @brief Compute minimum of scalar field, as well as the global index and level.
    template< typename Value >
    void minimumAndLocation( const Field&, Value& minimum, gidx_t& glb_idx ) const;

    /// @brief Compute maximum of scalar field, as well as the global index and level.
    template< typename Value >
    void maximumAndLocation( const Field&, Value& maximum, gidx_t& glb_idx ) const;

    /// @brief Compute minimum of scalar field, as well as the global index and level.
    template< typename Value >
    void minimumAndLocation( const Field&, Value& minimum, gidx_t& glb_idx, size_t& level ) const;

    /// @brief Compute maximum of scalar field, as well as the global index and level.
    template< typename Value >
    void maximumAndLocation( const Field&, Value& maximum, gidx_t& glb_idx, size_t& level ) const;

    /// @brief Compute minimum of field for each field-variable, as well as the global indices and levels.
    template< typename Value >
    void minimumAndLocation( const Field&, std::vector<Value>& minimum, std::vector<gidx_t>& glb_idx ) const;

    /// @brief Compute maximum of field for each field-variable, as well as the global indices and levels.
    template< typename Value >
    void maximumAndLocation( const Field&, std::vector<Value>& maximum, std::vector<gidx_t>& glb_idx ) const;

    /// @brief Compute minimum of field for each field-variable, as well as the global indices and levels.
    template< typename Value >
    void minimumAndLocation( const Field&, std::vector<Value>& minimum, std::vector<gidx_t>& glb_idx, std::vector<size_t>& level ) const;

    /// @brief Compute maximum of field for each field-variable, as well as the global indices and levels.
    template< typename Value >
    void maximumAndLocation( const Field&, std::vector<Value>& maximum, std::vector<gidx_t>& glb_idx, std::vector<size_t>& level ) const;

    /// @brief Compute minimum and its location of a field for each vertical level separately
    void minimumAndLocationPerLevel( const Field&, Field& column, Field& glb_idx ) const;

    /// @brief Compute maximum and its location of a field for each vertical level separately
    void maximumAndLocationPerLevel( const Field&, Field& column, Field& glb_idx ) const;

    /// @brief Compute mean value of scalar field
    /// @param [out] mean    Mean value
    /// @param [out] N       Number of value used to create the mean
    template< typename Value >
    void mean( const Field&, Value& mean, size_t& N ) const;

//    /// @brief Compute mean value of field for each field-variable
//    /// @param [out] mean    Mean values for each variable
//    /// @param [out] N       Number of values used to create the means
//    template< typename Value >
//    void mean( const Field&, std::vector<Value>& mean, size_t& N ) const;

    /// @brief Compute mean values of field for vertical level separately
    /// @param [out] mean    Field of dimension of input without the nodes index
    /// @param [out] N       Number of values used to create the means
    void meanPerLevel( const Field&, Field& mean, size_t& N ) const;

    /// @brief Compute mean value and standard deviation of scalar field
    /// @param [out] mean      Mean value
    /// @param [out] stddev    Standard deviation
    /// @param [out] N         Number of value used to create the mean
    template< typename Value >
    void meanAndStandardDeviation( const Field&, Value& mean, Value& stddev, size_t& N ) const;

//    /// @brief Compute mean values and standard deviations of scalar field for each field-variable
//    /// @param [out] mean      Mean values for each field-variable
//    /// @param [out] stddev    Standard deviation for each field-variable
//    /// @param [out] N         Number of value used to create the means
//    template< typename Value >
//    void meanAndStandardDeviation( const Field&, std::vector<Value>& mean, std::vector<Value>& stddev, size_t& N ) const;

    /// @brief Compute mean values and standard deviations of field for vertical level separately
    /// @param [out] mean      Field of dimension of input without the nodes index
    /// @param [out] stddev    Field of dimension of input without the nodes index
    /// @param [out] N         Number of values used to create the means
    void meanAndStandardDeviationPerLevel( const Field&, Field& mean, Field& stddev, size_t& N ) const;

private:

    const detail::NodeColumns* functionspace_;
};

// -------------------------------------------------------------------

inline size_t NodeColumns::levels() const {
    return functionspace_->levels();
}

template< typename Value >
void NodeColumns::sum( const Field& field, Value& sum, size_t& N ) const {
  functionspace_->sum(field,sum,N);
}

//template< typename Value >
//void NodeColumns::sum( const Field& field, std::vector<Value>& sum, size_t& N ) const {
//  functionspace_->sum(field,sum,N);
//}

inline void NodeColumns::sumPerLevel( const Field& field, Field& sum, size_t& N ) const {
  functionspace_->sumPerLevel(field,sum,N);
}

template< typename Value >
void NodeColumns::orderIndependentSum( const Field& field, Value& sum, size_t& N ) const {
  functionspace_->orderIndependentSum(field,sum,N);
}

//template< typename Value >
//void NodeColumns::orderIndependentSum( const Field& field, std::vector<Value>& sum, size_t& N ) const {
//  functionspace_->orderIndependentSum(field,sum,N);
//}

inline void NodeColumns::orderIndependentSumPerLevel( const Field& field, Field& sum, size_t& N ) const {
  functionspace_->orderIndependentSumPerLevel(field,sum,N);
}

template< typename Value >
void NodeColumns::minimum( const Field& field, Value& minimum ) const {
  functionspace_->minimum(field,minimum);
}

template< typename Value >
void NodeColumns::maximum( const Field& field, Value& maximum ) const {
  functionspace_->maximum(field,maximum);
}

//template< typename Value >
//void NodeColumns::minimum( const Field& field, std::vector<Value>& minimum ) const {
//  functionspace_->minimum(field,minimum);
//}

//template< typename Value >
//void NodeColumns::maximum( const Field& field, std::vector<Value>& maximum) const {
//  functionspace_->maximum(field,maximum);
//}

inline void NodeColumns::minimumPerLevel( const Field& field, Field& minimum ) const {
  return functionspace_->minimumPerLevel(field,minimum);
}

inline void NodeColumns::maximumPerLevel( const Field& field, Field& maximum ) const {
  functionspace_->maximumPerLevel(field,maximum);
}

template< typename Value >
void NodeColumns::minimumAndLocation( const Field& field, Value& minimum, gidx_t& glb_idx ) const {
  functionspace_->minimumAndLocation(field,minimum,glb_idx);
}

template< typename Value >
void NodeColumns::maximumAndLocation( const Field& field, Value& maximum, gidx_t& glb_idx ) const {
  functionspace_->maximumAndLocation(field,maximum,glb_idx);
}

template< typename Value >
void NodeColumns::minimumAndLocation( const Field& field, Value& minimum, gidx_t& glb_idx, size_t& level ) const {
  functionspace_->minimumAndLocation(field,minimum,glb_idx,level);
}

template< typename Value >
void NodeColumns::maximumAndLocation( const Field& field, Value& maximum, gidx_t& glb_idx, size_t& level ) const {
  functionspace_->maximumAndLocation(field,maximum,glb_idx,level);
}

template< typename Value >
void NodeColumns::minimumAndLocation( const Field& field, std::vector<Value>& minimum, std::vector<gidx_t>& glb_idx ) const {
  functionspace_->minimumAndLocation(field,minimum,glb_idx);
}

template< typename Value >
void NodeColumns::maximumAndLocation( const Field& field, std::vector<Value>& maximum, std::vector<gidx_t>& glb_idx ) const {
  functionspace_->maximumAndLocation(field,maximum,glb_idx);
}

template< typename Value >
void NodeColumns::minimumAndLocation( const Field& field, std::vector<Value>& minimum, std::vector<gidx_t>& glb_idx, std::vector<size_t>& level ) const {
  functionspace_->minimumAndLocation(field,minimum,glb_idx,level);
}

template< typename Value >
void NodeColumns::maximumAndLocation( const Field& field, std::vector<Value>& maximum, std::vector<gidx_t>& glb_idx, std::vector<size_t>& level ) const {
  functionspace_->maximumAndLocation(field,maximum,glb_idx,level);
}

inline void NodeColumns::minimumAndLocationPerLevel( const Field& field, Field& column, Field& glb_idx ) const {
  functionspace_->minimumAndLocationPerLevel(field,column,glb_idx);
}

inline void NodeColumns::maximumAndLocationPerLevel( const Field& field, Field& column, Field& glb_idx ) const {
  functionspace_->maximumAndLocationPerLevel(field,column,glb_idx);
}

template< typename Value >
void NodeColumns::mean( const Field& field, Value& mean, size_t& N ) const {
  functionspace_->mean(field,mean,N);
}

//template< typename Value >
//void NodeColumns::mean( const Field& field, std::vector<Value>& mean, size_t& N ) const {
//  functionspace_->mean(field,mean,N);
//}

inline void NodeColumns::meanPerLevel( const Field& field, Field& mean, size_t& N ) const {
  functionspace_->meanPerLevel(field,mean,N);
}

template< typename Value >
void NodeColumns::meanAndStandardDeviation( const Field& field, Value& mean, Value& stddev, size_t& N ) const {
  functionspace_->meanAndStandardDeviation(field,mean,stddev,N);
}

//template< typename Value >
//void NodeColumns::meanAndStandardDeviation( const Field& field, std::vector<Value>& mean, std::vector<Value>& stddev, size_t& N ) const {
//  functionspace_->meanAndStandardDeviation(field,mean,stddev,N);
//}

inline void NodeColumns::meanAndStandardDeviationPerLevel( const Field& field, Field& mean, Field& stddev, size_t& N ) const {
  functionspace_->meanAndStandardDeviationPerLevel(field,mean,stddev,N);
}

// -------------------------------------------------------------------

} // namespace functionspace
} // namespace atlas
