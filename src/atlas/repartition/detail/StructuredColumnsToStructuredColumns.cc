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

      // Get source range of this function space.
      const auto sourceRange = FuncSpaceRange(sourceStructuredColumnsPtr_);

      // Get target range of this function space.
      const auto targetRange = FuncSpaceRange(targetStructuredColumnsPtr_);

      // Get source ranges over all PEs.
      const auto sourceRanges = sourceRange.getFuncSpaceRanges();

      // Get get target ranges of all PEs.
      const auto targetRanges = targetRange.getFuncSpaceRanges();

      // Get intersections between sourceRange and targetRanges.
      auto getIntersections = [](const FuncSpaceRange& range,
        const FuncSpaceRangeVector& ranges){

          auto intersections = FuncSpaceRangeVector{};
          std::transform(ranges.cbegin(), ranges.cend(),
            std::back_inserter(intersections),
            [&](const FuncSpaceRange& rangesElem){
              return range & rangesElem;
            });

          return intersections;
        };

      auto sendIntersections = getIntersections(sourceRange, targetRanges);
      auto recvIntersections = getIntersections(targetRange, sourceRanges);

      // Trim off invalid intersections.
      auto trimIntersections = [](FuncSpaceRangeVector& intersections){

        intersections.erase(
          std::remove_if(intersections.begin(), intersections.end(),
            [](const FuncSpaceRange& intersection) -> bool {
              return !(intersection.getElemCount());
            }), intersections.end());
        };

      trimIntersections(sendIntersections);
      trimIntersections(recvIntersections);

      // Get block lengths and displacements for MPI data types.
      auto getBlockDataVector = [](const FuncSpaceRangeVector& intersections,
        const StructuredColumns* const StructuredColumnsPtr){

          auto blockDataVector = BlockDataVector{};
          std::transform(intersections.cbegin(), intersections.cend(),
          std::back_inserter(blockDataVector),
          [&](const FuncSpaceRange& intersection){
              return intersection.getBlockData(StructuredColumnsPtr);
            });

          return blockDataVector;
        };

      sendBlockDataVector_ =
        getBlockDataVector(sendIntersections, sourceStructuredColumnsPtr_);
      recvBlockDataVector_ =
        getBlockDataVector(recvIntersections, targetStructuredColumnsPtr_);

      // Make distributed graph communicator.
      auto getRanks = [](FuncSpaceRangeVector& intersections){

          auto ranks = std::vector<int>{};
          std::transform(intersections.cbegin(), intersections.cend(),
          std::back_inserter(ranks), [](const FuncSpaceRange& intersection){
            return intersection.getRank();
          });
          return ranks;
        };

      const auto oldComm = mpi::comm().communicator();
      const auto sources = getRanks(recvIntersections);
      const auto destinations = getRanks(sendIntersections);
      const auto inDegree = static_cast<int>(sources.size());
      const auto outDegree = static_cast<int>(destinations.size());

      checkMPI(MPI_Dist_graph_create_adjacent(oldComm, inDegree, sources.data(),
        MPI_UNWEIGHTED, outDegree, destinations.data(), MPI_UNWEIGHTED,
        MPI_INFO_NULL, false, &graphComm_), "MPI_Dist_graph_create_adjacent");

      // Set counts and displacement vectors.
      sendCounts_ = std::vector<int>(sendIntersections.size(), 1);
      sendDisplacements_ = std::vector<MPI_Aint>(sendIntersections.size(), 0);
      recvCounts_ = std::vector<int>(recvIntersections.size(), 1);
      recvDisplacements_ = std::vector<MPI_Aint>(recvIntersections.size(), 0);

      return;
    }

    // Destructor.
    StructuredColumnsToStructuredColumns::
      ~StructuredColumnsToStructuredColumns() {

      // Free graph
      checkMPI(MPI_Comm_free(&graphComm_), "MPI_Comm_free");

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

      // Get MPI indexed data types.
      auto getDataTypes = [&](const BlockDataVector& blockDataVector){

          auto dataTypes = std::vector<MPI_Datatype>(blockDataVector.size());
          std::transform(blockDataVector.cbegin(), blockDataVector.cend(),
            dataTypes.begin(), [&](const BlockData& blockData){
              MPI_Datatype dataType{};

              const auto& blockLengths = blockData.first;
              const auto& blockDisplacements = blockData.second;
              const auto blockCount = static_cast<int>(blockLengths.size());

              checkMPI(MPI_Type_indexed(
                blockCount, blockLengths.data(), blockDisplacements.data(),
                typeTraits<fieldType>::mpiType, &dataType), "MPI_Type_indexed");
              checkMPI(MPI_Type_commit(&dataType), "MPI_Type_commit");

              return dataType;
          });
          return dataTypes;
        };

      auto sendTypes = getDataTypes(sendBlockDataVector_);
      auto recvTypes = getDataTypes(recvBlockDataVector_);

      // Call MPI Neighbour all-to-all (no extra data copy).
      checkMPI(MPI_Neighbor_alltoallw(sourceView.data(), sendCounts_.data(),
        sendDisplacements_.data(), sendTypes.data(), targetView.data(),
        recvCounts_.data(), recvDisplacements_.data(), recvTypes.data(),
        graphComm_), "MPI_Neighbor_alltoallw");

      // Free data types.
      auto freeDataTypes = [](std::vector<MPI_Datatype>& dataTypes){

        std::for_each(dataTypes.begin(), dataTypes.end(),
          [](MPI_Datatype& dataType){
            checkMPI(MPI_Type_free(&dataType), "MPI_Type_free");
            return;
          });
        return;
      };

      freeDataTypes(sendTypes);
      freeDataTypes(recvTypes);

      return;
    }

    //==========================================================================
    // Index range methods.
    //==========================================================================

    // Constructor.
    FuncSpaceRange::FuncSpaceRange(
        const StructuredColumns* const structuredColumnsPtr) {

      rank_ = static_cast<int>(mpi::comm().rank());

      jBeginEnd_ = std::make_pair(
        structuredColumnsPtr->j_begin(), structuredColumnsPtr->j_end());

      for (auto j = jBeginEnd_.first; j < jBeginEnd_.second; ++j) {
        iBeginEnd_.push_back(std::make_pair(
          structuredColumnsPtr->i_begin(j), structuredColumnsPtr->i_end(j)));
      }
      return;
    }

    // Get index ranges from all PEs.
    FuncSpaceRangeVector FuncSpaceRange::getFuncSpaceRanges () const {

      // Get MPI communicator size.
      const auto mpiSize = static_cast<size_t>(mpi::comm().size());

      // Set recv buffer for j range.
      auto jRecvBuffer = idxPairVector(mpiSize);

      // Perform all gather.
      mpi::comm().allGather(
        jBeginEnd_, jRecvBuffer.begin(), jRecvBuffer.end());

      // Set i receive counts.
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
        iBeginEnd_.cbegin(), iBeginEnd_.cend(), irecvBuffer.begin(),
        irecvCounts.data(), irecvDisplacements.data());

      // Make vector of indexRange structs.
      auto indexRanges = FuncSpaceRangeVector{};
      for (size_t i = 0; i < mpiSize; ++i) {
        auto indexRange = FuncSpaceRange{};

        indexRange.rank_ = static_cast<int>(i);
        indexRange.jBeginEnd_ = jRecvBuffer[i];

        const auto iBegin = irecvBuffer.cbegin() + irecvDisplacements[i];
        const auto iEnd = iBegin + irecvCounts[i];
        std::copy(iBegin, iEnd, std::back_inserter(indexRange.iBeginEnd_));

        indexRanges.push_back(indexRange);
      }

      return indexRanges;
    }

    // Count number of elements in index range.
    idx_t FuncSpaceRange::getElemCount() const {

      // Accumulate size of positive i range.
      const idx_t count =
        std::accumulate(iBeginEnd_.cbegin(), iBeginEnd_.cend(),
        0, [](const int& cumulant, const idxPair& iElem) -> idx_t {

          // Only count positive differences.
          return cumulant + std::max(iElem.second - iElem.first, 0);

        });

      return count;
    }

    // Return the intersection between two index ranges.
    // Result takes the rank member of the right-hand argument.
    FuncSpaceRange FuncSpaceRange::operator&(
      const FuncSpaceRange& indexRange) const {

      // Declare result.
      auto intersection = FuncSpaceRange{};

      // get j intersection range.
      intersection.jBeginEnd_ = std::make_pair(
        std::max(jBeginEnd_.first, indexRange.jBeginEnd_.first),
        std::min(jBeginEnd_.second, indexRange.jBeginEnd_.second));

      // get i intersection range.
      if (intersection.jBeginEnd_.first < intersection.jBeginEnd_.second) {

        // get iterators.
        const auto iBeginA = iBeginEnd_.cbegin()
          + intersection.jBeginEnd_.first - jBeginEnd_.first;
        const auto iBeginB = indexRange.iBeginEnd_.cbegin()
          + intersection.jBeginEnd_.first - indexRange.jBeginEnd_.first;
        const auto iEndA = iBeginEnd_.cend()
          + intersection.jBeginEnd_.second - jBeginEnd_.second;

        std::transform(
          iBeginA, iEndA, iBeginB, std::back_inserter(intersection.iBeginEnd_),
          [](const idxPair& iElemA, const idxPair& iElemB){
            return std::make_pair(std::max(iElemA.first, iElemB.first),
              std::min(iElemA.second, iElemB.second));
          });
      }

      intersection.rank_ = indexRange.rank_;

      return intersection;
    }

    // Get MPI indexed type block lengths and displacements from an index
    // range.
    BlockData FuncSpaceRange::getBlockData(
      const StructuredColumns *const structuredColumnsPtr) const {

      // Declare results.
      auto blockLengths = std::vector<int>{};
      auto blockDisplacements = std::vector<int>{};

      // Get block data from index range.
      std::for_each(iBeginEnd_.cbegin(), iBeginEnd_.cend(),
        [&, j = jBeginEnd_.first](const idxPair& iElem) mutable {

          // Only add block if range is positive.
          if (iElem.second > iElem.first) {

            // Set block length.
            blockLengths.push_back(static_cast<int>(
              iElem.second - iElem.first) * structuredColumnsPtr->levels());

            // Set block displacement.
            blockDisplacements.push_back(static_cast<int>(
              structuredColumnsPtr->index(iElem.first, j)
              * structuredColumnsPtr->levels()));
          }

          ++j;
        });

      return std::make_pair(blockLengths, blockDisplacements);
    }

  }
}
