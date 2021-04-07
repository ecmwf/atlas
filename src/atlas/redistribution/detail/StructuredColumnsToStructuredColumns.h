/*
 * (C) Crown Copyright 2021 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <vector>

#include "atlas/functionspace/StructuredColumns.h"
#include "atlas/redistribution/detail/RedistributionImpl.h"


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
  namespace redistribution {
    namespace detail {

      // Forward declarations.
      class StructuredColumnsToStructuredColumns;
      class StructuredIndexRange;

      // Type aliases.
      class StructuredIndexRange;
      using idxPair = std::pair<idx_t, idx_t>;
      using idxPairVector = std::vector<idxPair>;
      using StructuredIndexRangeVector = std::vector<StructuredIndexRange>;

      using functionspace::detail::StructuredColumns;

      /// \brief    Concrete redistributor class for StructuredColumns to
      ///           StructuredColumns.
      ///
      /// \details  Class to map two function spaces with the same grid but
      ///           different partitioners.
      class StructuredColumnsToStructuredColumns : public RedistributionImpl {

      public:

        /// \brief    Constructs and initialises the redistributor.
        ///
        /// \details  Performs MPI_Allgatherv to determine the (i, j, k) ranges
        ///           of each source and target function space on each PE.
        ///           The grids of source and target function space must match.
        ///
        /// \param[in]  sourceFunctionSpace  Function space of source fields.
        /// \param[in]  targetFunctionSpace  Function space of target fields.
        StructuredColumnsToStructuredColumns(
          const FunctionSpace& sourceFunctionSpace,
          const FunctionSpace& targetFunctionSpace);

        /// \brief    Redistributes source field to target field.
        ///
        /// \details  Transfers source field to target field via an
        ///           MPI_Alltoallv. Function space of source field must match
        ///           sourceFunctionSpace supplied to the constructor. Same
        ///           applies to target field.
        ///
        /// \param[in]  sourceField  input field matching sourceFunctionSpace.
        /// \param[out] targetField  output field matching targetFunctionSpace.
        void execute(
          const Field& sourceField, Field& targetField) const override;

        /// \brief    Redistributes source field set to target fields set.
        ///
        /// \details  Transfers source field set to target field set via
        ///           multiple invocations of execute(sourceField, targetField).
        ///
        /// \param[in]  sourceFieldSet  input field set.
        /// \param[out] targetFieldSet  output field set.
        void execute(const FieldSet& sourceFieldSet,
          FieldSet& targetFieldSet) const override;

      private:

        // Generic execute call to handle different field types.
        template <typename fieldType>
        void doExecute(const Field& sourceField, Field& targetField) const;

        // FunctionSpaces recast to StructuredColumns.
        const StructuredColumns* sourceStructuredColumnsPtr_{};
        const StructuredColumns* targetStructuredColumnsPtr_{};

        // Vectors of index range intersection objects.
        StructuredIndexRangeVector sendIntersections_{};
        StructuredIndexRangeVector recvIntersections_{};

        // Counts and displacements for MPI communications.
        std::vector<int> sendCounts_{};
        std::vector<int> sendDisplacements_{};
        std::vector<int> recvCounts_{};
        std::vector<int> recvDisplacements_{};

      };

      /// \brief    Helper class for function space intersections.
      class StructuredIndexRange {

      public:

        /// \brief    Default Constructor.
        StructuredIndexRange() = default;

        /// \brief    Constructor.
        StructuredIndexRange(const StructuredColumns* const structuredColumnsPtr);

        /// \brief    Get index ranges from all PEs.
        StructuredIndexRangeVector getStructuredIndexRanges() const;

        /// \brief    Count number of elements.
        idx_t getElemCount() const;

        /// \brief    Intersection operator.
        StructuredIndexRange operator&(const StructuredIndexRange& indexRange) const;

        /// \brief    Iterate over all indices and do something with functor.
        template <typename functorType>
        void forEach(const functorType& functor) const;

      private:

        // Begin and end of j range.
        idxPair jBeginEnd_{};

        // Begin and end of i range for each j.
        idxPairVector iBeginEnd_{};

      };
    }
  }
}
