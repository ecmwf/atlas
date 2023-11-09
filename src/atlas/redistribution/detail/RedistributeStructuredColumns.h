/*
 * (C) Crown Copyright 2021 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <vector>

#include "atlas/redistribution/detail/RedistributionImpl.h"
#include "atlas/redistribution/detail/RedistributionImplFactory.h"


namespace atlas {

class Field;
class FieldSet;
class FunctionSpace;

namespace functionspace {
class StructuredColumns;
}  // namespace functionspace
}  // namespace atlas

namespace atlas {
namespace redistribution {
namespace detail {

// Forward declarations.
class RedistributeStructuredColumns;
class StructuredIndexRange;

// Type aliases.
class StructuredIndexRange;
using idxPair                    = std::pair<idx_t, idx_t>;
using idxPairVector              = std::vector<idxPair>;
using StructuredIndexRangeVector = std::vector<StructuredIndexRange>;

/// \brief    Concrete redistributor class for StructuredColumns to
///           StructuredColumns.
///
/// \details  Class to map two function spaces with the same grid but
///           different partitioners.
class RedistributeStructuredColumns : public RedistributionImpl {
public:
    RedistributeStructuredColumns() = default;

    /// \brief    Initialises the redistributor.
    void do_setup() override;

    static std::string static_type() { return "RedistributeStructuredColumns"; }

    std::string type() const override { return static_type(); }

    /// \brief    Redistributes source field to target field.
    ///
    /// \details  Transfers source field to target field via an
    ///           MPI_Alltoallv. Function space of source field must match
    ///           sourceFunctionSpace supplied to the constructor. Same
    ///           applies to target field.
    ///
    /// \param[in]  source  input field matching sourceFunctionSpace.
    /// \param[out] target  output field matching targetFunctionSpace.
    void execute(const Field& source, Field& target) const override;

    /// \brief    Redistributes source field set to target fields set.
    ///
    /// \details  Transfers source field set to target field set via
    ///           multiple invocations of execute(sourceField, targetField).
    ///
    /// \param[in]  source  input field set.
    /// \param[out] target  output field set.
    void execute(const FieldSet& source, FieldSet& target) const override;

private:
    // Generic execute call to handle different field types.
    template <typename fieldType>
    void do_execute(const Field& source, Field& target) const;

    // FunctionSpaces recast to StructuredColumns.
    functionspace::StructuredColumns source_;
    functionspace::StructuredColumns target_;

    // Vectors of index range intersection objects.
    StructuredIndexRangeVector sendIntersections_{};
    StructuredIndexRangeVector recvIntersections_{};

    // Counts and displacements for MPI communications.
    std::vector<int> sendCounts_{};
    std::vector<int> sendDisplacements_{};
    std::vector<int> recvCounts_{};
    std::vector<int> recvDisplacements_{};

    std::string mpi_comm_;
};

/// \brief    Helper class for function space intersections.
class StructuredIndexRange {
public:
    /// \brief    Default Constructor.
    StructuredIndexRange() = default;

    /// \brief    Constructor.
    StructuredIndexRange(const functionspace::StructuredColumns&);

    /// \brief    Get index ranges from all PEs.
    StructuredIndexRangeVector getStructuredIndexRanges() const;

    /// \brief    Count number of elements.
    idx_t getElemCount() const;

    /// \brief    Intersection operator.
    StructuredIndexRange operator&(const StructuredIndexRange&) const;

    /// \brief    Iterate over all indices and do something with functor.
    template <typename functorType>
    void forEach(const functorType&) const;

    const std::string& mpi_comm() const { return mpi_comm_; }

private:
    // Begin and end of j range.
    idxPair jBeginEnd_{};

    // Begin and end of i range for each j.
    idxPairVector iBeginEnd_{};

    std::string mpi_comm_;
};

}  // namespace detail
}  // namespace redistribution
}  // namespace atlas
