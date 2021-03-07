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

    //==========================================================================
    // Helper function forward declarations.
    //==========================================================================

    namespace  {
      IndexRangeVector
        getIndexRanges(const StructuredColumns* const structuredColumnsPtr);

      IndexRange getIntersection(const IndexRange& indexRangeA,
        const IndexRange& indexRangeB);

      IndexRangeVector getIntersections(const IndexRange& indexRangeA,
        const IndexRangeVector& indexRangesB);

      std::pair<std::vector<int>, std::vector<int>>
        getCountsDisplacements(const IndexRangeVector& indexRanges);

      void trimIntersections(IndexRangeVector& intersections,
        const std::vector<int>& counts);

      template <typename Functor>
      void indexRangeForEach(const IndexRange& indexRange, Functor functor);

    }

    //==========================================================================
    // Class public methods implementation.
    //==========================================================================

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
      sendIntersections_ =
        getIntersections(sourceRanges[mpiRank], targetRanges);

      // Calculate intersection between this target range and all source ranges.
      recvIntersections_ =
        getIntersections(targetRanges[mpiRank], sourceRanges);

      // Get send counts and displacements.
      std::tie(sendCounts_, sendDisplacements_) =
        getCountsDisplacements(sendIntersections_);

      // Get receive counts and displacements.
      std::tie(recvCounts_, recvDisplacements_) =
        getCountsDisplacements(recvIntersections_);

      // Trim off empty send intersections.
      trimIntersections(sendIntersections_, sendCounts_);

      // Trim off empty recv intersections.
      trimIntersections(recvIntersections_, recvCounts_);

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

        case array::DataType::KIND_REAL64:
          doExecute<double>(sourceField, targetField);
          break;

        case array::DataType::KIND_REAL32:
          doExecute<float>(sourceField, targetField);
          break;

        case array::DataType::KIND_INT32:
          doExecute<int>(sourceField, targetField);
          break;

        case array::DataType::KIND_INT64:
          doExecute<long>(sourceField, targetField);
          break;

        case array::DataType::KIND_UINT64:
          doExecute<unsigned long>(sourceField, targetField);
          break;

        default:
          throw eckit::NotImplemented(
            "No implementation for data type " +
            sourceField.datatype().str(), Here());
          break;
      }

      return;
    }

    void StructuredColumnsToStructuredColumns::execute(
      const FieldSet& sourceFieldSet, FieldSet& targetFieldSet) const {

      // Check that both FieldSets are the same size.
      checkFieldSetSize(sourceFieldSet, targetFieldSet,
        "sourceFieldSet", "targetFieldSet");

      std::for_each(sourceFieldSet.begin(), sourceFieldSet.end(),
        [&, i = idx_t{0}](const Field& sourceField) mutable {
          execute(sourceField, targetFieldSet[i++]);
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


      // Allocate send buffer.
      const auto nSend = static_cast<size_t>(
        sendCounts_.back() + sendDisplacements_.back());
      auto sendBuffer = std::vector<fieldType>(nSend);

      // Allocate receive buffer
      const auto nRecv = static_cast<size_t>(
        recvCounts_.back() + recvDisplacements_.back());
      auto recvBuffer = std::vector<fieldType>(nRecv);

      // The constructor should have made sure that the intersections match up
      // with counts and displacements vectors.

      //Loop over send intersections and write to send buffer.
      std::for_each(sendIntersections_.cbegin(), sendIntersections_.cend(),
        [&, sendBufferIt = sendBuffer.begin()]
        (const IndexRange& intersection) mutable {

        indexRangeForEach(intersection,
          [&](const idx_t i, const idx_t j, const idx_t k){
            const auto iNode = sourceStructuredColumnsPtr_->index(i, j);
            *sendBufferIt++ = sourceView(iNode, k);
            return;
          });

        return;

      });

      // Perform all to all communication.
      mpi::comm().allToAllv(sendBuffer.data(), sendCounts_.data(),
        sendDisplacements_.data(), recvBuffer.data(), recvCounts_.data(),
        recvDisplacements_.data());


      //Loop over recv intersections and read from recv buffer.
      std::for_each(recvIntersections_.cbegin(), recvIntersections_.cend(),
        [&, recvBufferIt = recvBuffer.cbegin()]
        (const IndexRange& intersection) mutable {

        indexRangeForEach(intersection,
          [&](const idx_t i, const idx_t j, const idx_t k){
            const auto iNode = targetStructuredColumnsPtr_->index(i, j);
            targetView(iNode, k) = *recvBufferIt++;
            return;
          });

        return;

      });

      // By this point, everything should have worked.

      return;
    }

    //==========================================================================
    // Helper function definitions.
    //==========================================================================

    namespace  {

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

      // Get the count and displacement arrays needed for an allToAll MPI
      // communication.
      std::pair<std::vector<int>, std::vector<int>>
        getCountsDisplacements(const IndexRangeVector& indexRanges) {

        // Declare result.
        auto counts = std::vector<int>{};
        auto displacements = std::vector<int>{};

        // Count the number of elements in each index range.
        std::transform(indexRanges.cbegin(), indexRanges.cend(),
          std::back_inserter(counts),
          [](const IndexRange& rangeElem){

            // Accumulate size of positive i range.
            const int count = std::accumulate(rangeElem.iBeginEnd.cbegin(),
              rangeElem.iBeginEnd.cend(), 0,
              [](const int& cumulant, const idxPair& iElem) -> int {
                return cumulant + std::max(iElem.second - iElem.first, 0);
              }) * rangeElem.levels;

            return count;
          });

        // Calculate displacements from counts.
        displacements.push_back(0);
        std::partial_sum(counts.cbegin(), counts.cend() - 1,
          std::back_inserter(displacements));

        return std::make_pair(counts, displacements);
      }

      // Trim off index ranges where counts equals zero.
      void trimIntersections(IndexRangeVector& intersections,
        const std::vector<int>& counts) {

        intersections.erase(
          std::remove_if(intersections.begin(), intersections.end(),
          [&](const IndexRange& intersection) -> bool {

            // Some pointer arithmetic to get element of counts vector.
            const auto distance = &intersection - intersections.data();

            return !(*(counts.data() + distance));
          }), intersections.end());
        return;
      }

      // Visit all indices and call f(idx_t i, idx_t j, idx_t k).
      template <typename functor>
      void indexRangeForEach(const IndexRange& indexRange, functor f) {

        idx_t jBegin, jEnd;
        std::tie(jBegin, jEnd) = indexRange.jBeginEnd;

        for (idx_t j = jBegin; j < jEnd; ++j) {

          idx_t iBegin, iEnd;
          std::tie(iBegin, iEnd) =
            indexRange.iBeginEnd[static_cast<size_t>(j - jBegin)];

          for (idx_t i = iBegin; i < iEnd; ++i){

            for (idx_t k = 0; k < indexRange.levels; ++k) {

              // Call functor.
              f(i, j, k);

            }
          }
        }

        return;
      }

    }
  }
}
