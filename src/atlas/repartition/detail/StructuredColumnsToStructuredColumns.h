/*
 * (C) Crown Copyright 2021 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <mpi.h>
#include <vector>

#include "atlas/functionspace/StructuredColumns.h"
#include "atlas/repartition/detail/RepartitionImpl.h"

namespace atlas {

  class Field;
  class FieldSet;
  class FunctionSpace;

  namespace functionspace {
    namespace detail {
      class StructuredColumns;
    }
  }
}

namespace atlas {
  namespace repartition {

    // Forward declarations.
    class StructuredColumnsToStructuredColumns;
    class FuncSpaceRange;

    // Type aliases.
    struct FuncSpaceRange;
    using idxPair = std::pair<idx_t, idx_t>;
    using idxPairVector = std::vector<idxPair>;
    using FuncSpaceRangeVector = std::vector<FuncSpaceRange>;
    using BlockData = std::pair<std::vector<int>, std::vector<int>>;
    using BlockDataVector = std::vector<BlockData>;

    // Atlas to MPI type conversion.
    template <typename T>
    struct typeTraits {};

    template <>
    struct typeTraits<double> {
      static const auto atlasType = array::DataType::KIND_REAL64;
      static const auto mpiType = MPI_DOUBLE;
    };

    template <>
    struct typeTraits<float> {
      static const auto atlasType = array::DataType::KIND_REAL32;
      static const auto mpiType = MPI_FLOAT;
    };

    template <>
    struct typeTraits<int> {
      static const auto atlasType = array::DataType::KIND_INT32;
      static const auto mpiType = MPI_INT;
    };

    template <>
    struct typeTraits<long> {
      static const auto atlasType = array::DataType::KIND_INT64;
      static const auto mpiType = MPI_LONG;
    };

    template <>
    struct typeTraits<unsigned long> {
      static const auto atlasType = array::DataType::KIND_UINT64;
      static const auto mpiType = MPI_UNSIGNED_LONG;
    };

    using functionspace::detail::StructuredColumns;

    /// \brief    Concrete repartitioning class for StructuredColumns to
    ///           StructuredColumns.
    ///
    /// \details  Class to map two function spaces with the same grid but
    ///           different partitioners. Creates an MPI Graph communicator
    ///           to perform an efficient, scalable Neighbor_Alltoallw
    ///           communication.
    class StructuredColumnsToStructuredColumns : public RepartitionImpl {

    public:

      /// \brief    Constructs and initialises the repartitioner.
      ///
      /// \details  Performs MPI_Allgatherv to determine the (i, j, k) ranges of
      ///           each source and target function space on each PE. Then
      ///           calculates range intersections to determine what needs to be
      ///           sent where via an MPI_Alltoallw. The grids of source and
      ///           target function space must match.
      ///
      /// \param[in]  sourceFunctionSpace  Function space of source fields.
      /// \param[in]  targetFunctionSpace  Function space of target fields.
      StructuredColumnsToStructuredColumns(
        const FunctionSpace& sourceFunctionSpace,
        const FunctionSpace& targetFunctionSpace);

      /// \brief    Destructor.
      ~StructuredColumnsToStructuredColumns() override;

      /// \brief    Repartitions source field to target field.
      ///
      /// \details  Transfers source field to target field via an
      ///           MPI_Neighbor_Alltoallw. Function space of source field
      ///           must match sourceFunctionSpace supplied to the constructor.
      ///           Same applies to target field.
      ///
      /// \param[in]  sourceField  input field matching sourceFunctionSpace.
      /// \param[out] targetField  output field matching targetFunctionSpace.
      void execute(const Field& sourceField, Field& targetField) const override;

      /// \brief    Repartitions source field set to target fields set.
      ///
      /// \details  Transfers source field set to target field set via multiple
      ///           invocations of execute(sourceField, targetField).
      ///
      /// \param[in]  sourceFieldSet  input field set.
      /// \param[out] targetFieldSet  output field set.
      void execute(const FieldSet& sourceFieldSet,
        FieldSet& targetFieldSet) const override;

      // Disable copying.
      StructuredColumnsToStructuredColumns(
        const StructuredColumnsToStructuredColumns&) = delete;
      StructuredColumnsToStructuredColumns operator=(
      const StructuredColumnsToStructuredColumns&) = delete;


    private:

      // Generic execute call to handle different field types.
      template <typename fieldType>
      void doExecute(const Field& sourceField, Field& targetField) const;

      // FunctionSpaces recast to StructuredColumns.
      const StructuredColumns* sourceStructuredColumnsPtr_{};
      const StructuredColumns* targetStructuredColumnsPtr_{};

      // MPI distributed graph communicator.
      MPI_Comm graphComm_{};

      // Block lenghts and displacements for MPI indexed datatype.
      BlockDataVector sendBlockDataVector_{};
      BlockDataVector recvBlockDataVector_{};

      // Data counts and displacements.
      std::vector<int> sendCounts_{};
      std::vector<MPI_Aint> sendDisplacements_{};
      std::vector<int> recvCounts_{};
      std::vector<MPI_Aint> recvDisplacements_{};

    };

    // Helper class for function space intersections.
    class FuncSpaceRange {

    public:

      // Default Constructor.
      FuncSpaceRange() = default;

      // Constructor.
      FuncSpaceRange(const StructuredColumns* const structuredColumnsPtr);

      // Get index ranges from all PEs.
      FuncSpaceRangeVector getFuncSpaceRanges() const;

      // Count number of elements.
      idx_t getElemCount() const;

      // Intersection operator.
      FuncSpaceRange operator&(const FuncSpaceRange& indexRange) const;

      // Get block lengths and displacements for index range.
      BlockData
        getBlockData(
        const StructuredColumns* const structuredColumnsPtr) const;

      // Return copy of rank.
      int getRank() const {return rank_;}

    private:

      // MPI rank of range.
      int rank_{};

      // Begin and end of j range.
      idxPair jBeginEnd_{};

      // Begin and end of i range for each j.
      idxPairVector iBeginEnd_{};

    };
  }
}
