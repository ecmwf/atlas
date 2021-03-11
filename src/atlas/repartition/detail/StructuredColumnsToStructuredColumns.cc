/*
 * (C) Crown Copyright 2021 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <algorithm>
#include <numeric>
#include <string>

#include "atlas/array/MakeView.h"
#include "atlas/field.h"
#include "atlas/field/FieldSet.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/repartition/detail/StructuredColumnsToStructuredColumns.h"
#include "atlas/repartition/detail/RepartitionUtils.h"

namespace atlas {
  namespace repartition {

    namespace  {

      //========================================================================
      // Type traits structs.
      //========================================================================

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

      //========================================================================
      // Helper function forward declarations.
      //========================================================================

      idx_t countElements(const IndexRange& indexRange);

      IndexRangeVector
        getIndexRanges(const StructuredColumns* const structuredColumnsPtr);

      IndexRange getIntersection(const IndexRange& indexRangeA,
        const IndexRange& indexRangeB);

      IndexRangeVector getIntersections(const IndexRange& indexRangeA,
        const IndexRangeVector& indexRangesB);

      MPI_Comm getGraphComm(const IndexRangeVector& sendIntersections,
        const IndexRangeVector& recvIntersections);

      std::pair<std::vector<std::vector<int>>, std::vector<std::vector<int>>>
        getBlockLengthsDisplacements(const IndexRangeVector& intersections,
        const StructuredColumns* const structuredColumnsPtr);

      std::vector<MPI_Datatype> getDataTypes(
        const std::vector<std::vector<int>>& blockLengths,
        const std::vector<std::vector<int>>& blockDisplacements,
        const MPI_Datatype elemDataType);

      void freeDataTypes(std::vector<MPI_Datatype>& dataTypes);

    }

    //==========================================================================
    // Class public methods implementation.
    //==========================================================================

    // Constructor.
    StructuredColumnsToStructuredColumns::StructuredColumnsToStructuredColumns(
      const FunctionSpace& sourceFunctionSpace,
      const FunctionSpace& targetFunctionSpace) :
      RepartitionImpl(sourceFunctionSpace, targetFunctionSpace),
      sourceStructuredColumnsPtr_(
        getSourceFunctionSpace()->cast<StructuredColumns>()),
      targetStructuredColumnsPtr_(
        getTargetFunctionSpace()->cast<StructuredColumns>()) {

      // Check that grids match (will also check for bad casts).
      checkGrids<StructuredColumns>(
        sourceStructuredColumnsPtr_, targetStructuredColumnsPtr_,
        "sourceStructuredColumnsPtr_", "targetStructuredColumnsPtr_");

      // Get MPI rank.
      const auto mpiRank = static_cast<size_t>(mpi::comm().rank());

      // Get ranges of source function spaces.
      const auto sourceRanges = getIndexRanges(sourceStructuredColumnsPtr_);

      // Get ranges of target function spaces.
      const auto targetRanges = getIndexRanges(targetStructuredColumnsPtr_);

      // Calculate intersection between this source range and all target ranges.
      const auto sendIntersections =
        getIntersections(sourceRanges[mpiRank], targetRanges);

      // Calculate intersection between this target range and all source ranges.
      const auto recvIntersections =
        getIntersections(targetRanges[mpiRank], sourceRanges);

      // Get distrubted graph communicator.
      graphComm_ = getGraphComm(sendIntersections, recvIntersections);

      // Get send block lengths and displacements.
      std::tie(sendBlockLengths_, sendBlockDisplacements_) =
        getBlockLengthsDisplacements(sendIntersections,
          sourceStructuredColumnsPtr_);

      // Get receive block lengths and displacements.
      std::tie(recvBlockLengths_, recvBlockDisplacements_) =
        getBlockLengthsDisplacements(recvIntersections,
          targetStructuredColumnsPtr_);

      // Set counts and displacement vectors.
      sendCounts_ = std::vector<int>(sendBlockLengths_.size(), 1);
      sendDisplacements_ = std::vector<MPI_Aint>(sendBlockLengths_.size(), 0);
      recvCounts_ = std::vector<int>(recvBlockLengths_.size(), 1);
      recvDisplacements_ = std::vector<MPI_Aint>(recvBlockLengths_.size(), 0);

      return;
    }

    // Destructor.
    StructuredColumnsToStructuredColumns::
      ~StructuredColumnsToStructuredColumns() {

      // Free graph
      MPI_Comm_free(&graphComm_);

      return;
    }

    void StructuredColumnsToStructuredColumns::execute(
      const Field& sourceField, Field& targetField) const {

      // Check source grids match.
      checkGrids<StructuredColumns>(
        sourceField.functionspace().get(), sourceStructuredColumnsPtr_,
        "sourceField.functionspace()", "sourceStructuredColumnsPtr_");

      // Check target grids match.
      checkGrids<StructuredColumns>(
        targetField.functionspace().get(), targetStructuredColumnsPtr_,
        "targetField.functionspace()", "targetStructuredColumnsPtr_");

      // Check data types match.
      checkFieldDataType(sourceField, targetField,
        "sourceField", "targetField");

      // Determine data type of field and execute.
      switch(sourceField.datatype().kind()){

        case typeTraits<double>::atlasType:
          doExecute<double>(sourceField, targetField);
          break;

        case typeTraits<float>::atlasType:
          doExecute<float>(sourceField, targetField);
          break;

        case typeTraits<int>::atlasType:
          doExecute<int>(sourceField, targetField);
          break;

        case typeTraits<long>::atlasType:
          doExecute<long>(sourceField, targetField);
          break;

        case typeTraits<unsigned long>::atlasType:
          doExecute<unsigned long>(sourceField, targetField);
          break;

        default:
          throw eckit::NotImplemented(
            "No implementation for data type " +
            sourceField.datatype().str(), Here());
      }

      return;
    }

    void StructuredColumnsToStructuredColumns::execute(
      const FieldSet& sourceFieldSet, FieldSet& targetFieldSet) const {

      // Check that both FieldSets are the same size.
      checkFieldSetSize(sourceFieldSet, targetFieldSet,
        "sourceFieldSet", "targetFieldSet");

      std::for_each(sourceFieldSet.cbegin(), sourceFieldSet.cend(),
        [&, targetFieldSetIt = targetFieldSet.begin()]
        (const Field& sourceField) mutable {
          execute(sourceField, *targetFieldSetIt++);
          return;
        });

      return;
    }

    //==========================================================================
    // Class private methods implementation.
    //==========================================================================

    template<typename fieldType>
    void StructuredColumnsToStructuredColumns::doExecute(
      const Field& sourceField, Field& targetField) const {

      // Make Atlas view objects.
      const auto sourceView = array::make_view<fieldType, 2>(sourceField);
      auto targetView = array::make_view<fieldType, 2>(targetField);

      // Make vector of MPI send data types.
      auto sendTypes = getDataTypes(sendBlockLengths_, sendBlockDisplacements_,
        typeTraits<fieldType>::mpiType);

      // Make vector of MPI recv data types.
      auto recvTypes = getDataTypes(recvBlockLengths_, recvBlockDisplacements_,
        typeTraits<fieldType>::mpiType);

      // Call MPI_Alltoallw.
      MPI_Neighbor_alltoallw(sourceView.data(), sendCounts_.data(),
        sendDisplacements_.data(), sendTypes.data(), targetView.data(),
        recvCounts_.data(), recvDisplacements_.data(), recvTypes.data(),
        graphComm_);

      // Free data types.
      freeDataTypes(sendTypes);
      freeDataTypes(recvTypes);

      return;
    }

    //==========================================================================
    // Helper function definitions.
    //==========================================================================

    namespace  {

      // Count the number of elements in an indexRange
      idx_t countElements(const IndexRange& indexRange) {

        // Accumulate size of positive i range.
        const idx_t count = std::accumulate(indexRange.iBeginEnd.cbegin(),
          indexRange.iBeginEnd.cend(), 0,
          [](const int& cumulant, const idxPair& iElem) -> idx_t {

            // Only count positive differences.
            return cumulant + std::max(iElem.second - iElem.first, 0);

          }) * indexRange.levels;

        return count;
      }

      // Communicate index range of function spaces to all PEs
      IndexRangeVector getIndexRanges(
        const StructuredColumns* const structuredColumnsPtr) {

        // Get MPI communicator size.
        const auto mpiSize = static_cast<size_t>(mpi::comm().size());

        // Set send buffer for j range.
        const auto jSendBuffer = std::make_pair(
          structuredColumnsPtr->j_begin(), structuredColumnsPtr->j_end());

        // Set recv buffer for j range.
        auto jRecvBuffer = idxPairVector(mpiSize);

        // Perform all gather.
        mpi::comm().allGather(
          jSendBuffer, jRecvBuffer.begin(), jRecvBuffer.end());

        // Set send buffer for i range.
        auto iSendBuffer = idxPairVector{};
        std::generate_n(std::back_inserter(iSendBuffer),
          jSendBuffer.second - jSendBuffer.first,
          [&, j = jSendBuffer.first]() mutable {
            const auto iElem = std::make_pair(
              structuredColumnsPtr->i_begin(j), structuredColumnsPtr->i_end(j));
            ++j;
            return iElem;
          });

        // Set recv counts for i range.
        auto irecvCounts = std::vector<int>{};
        std::transform(jRecvBuffer.cbegin(), jRecvBuffer.cend(),
          std::back_inserter(irecvCounts),
          [](const idxPair& jElem){
            return static_cast<int>(jElem.second - jElem.first);
          });

        // Set recv displacements for i range.
        auto irecvDisplacements = std::vector<int>{0};
        std::partial_sum(irecvCounts.cbegin(), irecvCounts.cend() - 1,
          std::back_inserter(irecvDisplacements));

        // Set recv buffer for i range.
        auto irecvBuffer = idxPairVector(static_cast<size_t>(
          irecvDisplacements.back() + irecvCounts.back()));

        // Perform all gather.
        mpi::comm().allGatherv(
          iSendBuffer.cbegin(), iSendBuffer.cend(), irecvBuffer.begin(),
          irecvCounts.data(), irecvDisplacements.data());

        // Make vector of indexRange structs.
        auto indexRanges = IndexRangeVector{};
        std::generate_n(std::back_inserter(indexRanges), mpiSize,
          [&, i = size_t{0}]() mutable {

            auto rangeElem = IndexRange{};
            rangeElem.levels = structuredColumnsPtr->levels();
            rangeElem.jBeginEnd = jRecvBuffer[i];

            const auto iBegin =
              irecvBuffer.cbegin() + irecvDisplacements[i];
            const auto iEnd = iBegin + irecvCounts[i];
            std::copy(iBegin, iEnd, std::back_inserter(rangeElem.iBeginEnd));

            // Count number of elements.
            rangeElem.nElem = countElements(rangeElem);

            ++i;
            return rangeElem;
          });

        return indexRanges;
      }

      // Calculate the intersection between two index ranges.
      IndexRange getIntersection(const IndexRange& indexRangeA,
        const IndexRange& indexRangeB) {

        // Declare result.
        auto intersection = IndexRange{};

        // set number of levels.
        intersection.levels = indexRangeA.levels;

        // get j intersection range.
        intersection.jBeginEnd = std::make_pair(
          std::max(indexRangeA.jBeginEnd.first, indexRangeB.jBeginEnd.first),
          std::min(indexRangeA.jBeginEnd.second, indexRangeB.jBeginEnd.second));

        // get i intersection range.
        if (intersection.jBeginEnd.first < intersection.jBeginEnd.second) {

          // get iterators.
          const auto iBeginA = indexRangeA.iBeginEnd.cbegin()
            + intersection.jBeginEnd.first - indexRangeA.jBeginEnd.first;
          const auto iBeginB = indexRangeB.iBeginEnd.cbegin()
            + intersection.jBeginEnd.first - indexRangeB.jBeginEnd.first;
          const auto iEndA = indexRangeA.iBeginEnd.cend()
            + intersection.jBeginEnd.second - indexRangeA.jBeginEnd.second;

          std::transform(
            iBeginA, iEndA, iBeginB, std::back_inserter(intersection.iBeginEnd),
            [](const idxPair& iElemA, const idxPair& iElemB){
              return std::make_pair(std::max(iElemA.first, iElemB.first),
                std::min(iElemA.second, iElemB.second));
            });
        }

        // Count number of elements.
        intersection.nElem = countElements(intersection);

        return intersection;
      }

      // Calculate the intersections between an index range and a vector of
      // index ranges.
      IndexRangeVector getIntersections(const IndexRange& indexRangeA,
        const IndexRangeVector& indexRangesB) {

        // Declare result.
        auto intersections = IndexRangeVector{};

        // Calculate intersection between index ranges.
        std::transform(indexRangesB.cbegin(), indexRangesB.cend(),
          std::back_inserter(intersections), [&](const IndexRange& indexRangeB){
            return getIntersection(indexRangeA, indexRangeB);
          });

        return intersections;
      }

      // Returns a handle to an MPI directed graph for neighbourhood
      // communications.
      MPI_Comm getGraphComm(const IndexRangeVector& sendIntersections,
        const IndexRangeVector& recvIntersections) {

        // Declare result
        MPI_Comm graphComm{};

        // Get parent communicator.
        const auto oldComm = mpi::comm().communicator();

        // Get ranks of valid receive PEs.
        std::vector<int> sources{};
        std::for_each(recvIntersections.cbegin(), recvIntersections.cend(),
          [&, i = int{0}](const IndexRange& intersection) mutable {
            if (intersection.nElem) sources.push_back(i);
            ++i;
            return;
          });

        const auto inDegree = static_cast<int>(sources.size());

        // Get ranks of valid send PEs.
        std::vector<int> destinations{};
        std::for_each(sendIntersections.cbegin(), sendIntersections.cend(),
          [&, i = int{0}](const IndexRange& intersection) mutable {
            if (intersection.nElem) destinations.push_back(i);
            ++i;
            return;
          });

        const auto outDegree = static_cast<int>(destinations.size());

        // Create distributed graph communicator.
        MPI_Dist_graph_create_adjacent(oldComm, inDegree, sources.data(),
          MPI_UNWEIGHTED, outDegree, destinations.data(), MPI_UNWEIGHTED,
          MPI_INFO_NULL, false, &graphComm);

        return graphComm;
      }

      std::pair<std::vector<std::vector<int>>, std::vector<std::vector<int>>>
        getBlockLengthsDisplacements(const IndexRangeVector& intersections,
        const StructuredColumns* const structuredColumnsPtr) {

        // Declare results.
        auto blockLengths = std::vector<std::vector<int>>{};
        auto blockDisplacements = std::vector<std::vector<int>>{};

        // Loop over all intersections.
        std::for_each(intersections.cbegin(), intersections.cend(),
          [&](const IndexRange& intersection){

            // Add block information if intersection is valid.
            if (intersection.nElem) {

              // Add length and displacement vectors.
              blockLengths.push_back(std::vector<int>{});
              blockDisplacements.push_back(std::vector<int>{});

              // Get reference to newly added vectors.
              auto& blockLength = blockLengths.back();
              auto& blockDisplacement = blockDisplacements.back();

              // Loop over all i range pairs.
              std::for_each(intersection.iBeginEnd.cbegin(),
                intersection.iBeginEnd.cend(),
                [&, j = intersection.jBeginEnd.first]
                (const idxPair& iElem) mutable {

                  auto iSize = iElem.second - iElem.first;
                  if (iSize > 0) {

                    // Add block size.
                    blockLength.push_back(iSize * intersection.levels);

                    // Add displacement.
                    blockDisplacement.push_back(
                      structuredColumnsPtr->index(iElem.first, j)
                      * intersection.levels);
                  }

                  ++j;
                  return;
                });

            }

            return;
          });

        return std::make_pair(blockLengths, blockDisplacements);
      }

      // Get a vector of data types from arrays of block lengths, displacements
      // and an element data type.
      std::vector<MPI_Datatype> getDataTypes(
        const std::vector<std::vector<int>>& blockLengths,
        const std::vector<std::vector<int>>& blockDisplacements,
        const MPI_Datatype elemDataType) {

        // Declare result.
        auto dataTypes = std::vector<MPI_Datatype>(blockLengths.size());

        // Generate data types.
        std::for_each(dataTypes.begin(), dataTypes.end(),
          [&, blockLengthIt = blockLengths.begin(),
          blockDisplacementIt = blockDisplacements.begin()]
          (MPI_Datatype& sendType) mutable {

            // Create and commit MPI derived type.
            auto blockCount = static_cast<int>(blockLengthIt->size());
            MPI_Type_indexed(blockCount, blockLengthIt->data(),
              blockDisplacementIt->data(), elemDataType,
              &sendType);

            MPI_Type_commit(&sendType);
            ++blockLengthIt;
            ++blockDisplacementIt;

            return;
          });

        return dataTypes;
      }

      // Free a vector of data types.
      void freeDataTypes(std::vector<MPI_Datatype>& dataTypes) {

        std::for_each(dataTypes.begin(), dataTypes.end(),
          [](MPI_Datatype& dataType){
            MPI_Type_free(&dataType);
          });

        return;
      }

    }
  }
}
