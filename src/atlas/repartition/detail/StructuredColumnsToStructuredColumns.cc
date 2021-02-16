#include <algorithm>
#include <numeric>

#include "atlas/array/MakeView.h"
#include "atlas/field.h"
#include "atlas/field/FieldSet.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/repartition/detail/StructuredColumnsToStructuredColumns.h"

namespace atlas {
  namespace repartition {

    // Helper function forward declarations.
    namespace  {
      std::vector<indexRange>
        getIndexRange(const StructuredColumns& functionSpace);

      indexRange getIntersection(const indexRange& indexRangeA,
        const indexRange& indexRangeB);

      std::vector<indexRange> getIntersections(const indexRange& indexRangeA,
        const std::vector<indexRange>& indexRangesB);

      std::pair<std::vector<int>, std::vector<int>>
        getCountsDisplacements(const std::vector<indexRange>& indexRanges);

      std::vector<indexRange> trimIntersections(
        const std::vector<indexRange>& intersections,
        const std::vector<int>& counts);

      template <typename fieldType>
      void performAllToAll(const Field& sourceField, const Field& targetField);
    }


    StructuredColumnsToStructuredColumns::StructuredColumnsToStructuredColumns(
      const FunctionSpace& sourceFunctionSpace,
      const FunctionSpace& targetFunctionSpace) :
      RepartitionImpl(sourceFunctionSpace, targetFunctionSpace),
      sourceStructuredColumns_(getSourceFunctionSpace()),
      targetStructuredColumns_(getTargetFunctionSpace()) {

      // Get MPI rank.
      auto mpiRank = static_cast<size_t>(mpi::comm().rank());

      // Get ranges of source function spaces.
      auto sourceRanges = getIndexRange(sourceStructuredColumns_);

      // Get ranges of target function spaces.
      auto targetRanges = getIndexRange(targetStructuredColumns_);

      // Calculate intersection between this source range and all target ranges.
      sendIntersections_ =
        getIntersections(sourceRanges[mpiRank], targetRanges);

      // Calculate intersection between this target range and all source ranges.
      receiveIntersections_ =
        getIntersections(targetRanges[mpiRank], sourceRanges);

      // Get send counts and displacements.
      std::tie(sendCounts_, sendDisplacements_) =
        getCountsDisplacements(sendIntersections_);

      // Get receive counts and displacements.
      std::tie(receiveCounts_, receiveDisplacements_) =
        getCountsDisplacements(receiveIntersections_);

      // Trim off empty send intersections.
      sendIntersections_ =
        trimIntersections(sendIntersections_, sendCounts_);

      // Trim off empty receive intersections.
      receiveIntersections_ =
        trimIntersections(receiveIntersections_, receiveCounts_);

      if (mpiRank == 0) {

          std::for_each(sendIntersections_.begin(), sendIntersections_.end(),
            [&](indexRange& rangeElem){
              std::cout << std::endl << "MPI rank " << mpiRank << std::endl;
              std::cout << rangeElem.jBeginEnd.first << " " << rangeElem.jBeginEnd.second << std::endl;
              std::cout << std::endl;
              std::for_each(rangeElem.iBeginEnd.begin(),rangeElem.iBeginEnd.end(),
                [](idxPair elem){
                std::cout << elem.first << " " << elem.second << std::endl;
              });
            });
        std::cout << sendCounts_ << std::endl;
      }

      return;
    }

    void StructuredColumnsToStructuredColumns::execute(
      const Field& sourceField, Field& targetField ) const {

      std::cout<< "Field execute on MPI rank " << mpi::rank()<< std::endl;

      return;
    }

    void StructuredColumnsToStructuredColumns::execute(
      const FieldSet& sourceField, FieldSet& targetField) const {

      std::cout<< "FieldSet execute on MPI rank " << mpi::rank()<< std::endl;

      return;
    }


    // Helper function definitions.
    namespace  {

      // Communicate index range of function spaces to all PEs
      std::vector<indexRange> getIndexRange(
        const StructuredColumns& functionSpace) {

        // Get MPI communicator size.
        auto mpiSize = static_cast<size_t>(mpi::comm().size());

        // Set send buffer for j range.
        auto jSendBuffer = idxPair(
          functionSpace.j_begin(), functionSpace.j_end());

        // Set receive buffer for j range.
        auto jReceiveBuffer = idxPairVector(mpiSize);

        // Perform all gather.
        mpi::comm().allGather(
          jSendBuffer, jReceiveBuffer.begin(), jReceiveBuffer.end());

        // Set send buffer for i range.
        auto iSendBuffer = idxPairVector();
        std::generate_n(std::back_inserter(iSendBuffer),
          jSendBuffer.second - jSendBuffer.first,
          [&, j = jSendBuffer.first]() mutable
          {
            auto iElem = idxPair(
              functionSpace.i_begin(j), functionSpace.i_end(j));
            ++j;
            return iElem;
          });

        // Set receive counts for i range.
        auto iReceiveCounts = std::vector<int>{};
        std::transform(jReceiveBuffer.begin(), jReceiveBuffer.end(),
          std::back_inserter(iReceiveCounts),
          [](const idxPair& jElem){
            return static_cast<int>(jElem.second - jElem.first);
          });

        // Set receive displacements for i range.
        auto iReceiveDisplacements = std::vector<int>{1};
        std::partial_sum(iReceiveCounts.begin(), iReceiveCounts.end() - 1,
          std::back_inserter(iReceiveDisplacements));

        // Set receive buffer for i range.
        auto iReceiveBuffer = idxPairVector(static_cast<size_t>(
          iReceiveDisplacements.back() + iReceiveCounts.back()));

        // Perform all gather.
        mpi::comm().allGatherv(
          iSendBuffer.begin(), iSendBuffer.end(), iReceiveBuffer.begin(),
          iReceiveCounts.data(), iReceiveDisplacements.data());

        // Make vector of indexRange structs.
        auto indexRanges = std::vector<indexRange>{};
        std::generate_n(std::back_inserter(indexRanges), mpiSize,
          [&, i = size_t{0}]() mutable {

            auto rangeElem = indexRange{};
            rangeElem.jBeginEnd = jReceiveBuffer[i];

            auto iBegin = iReceiveBuffer.begin() + iReceiveDisplacements[i];
            auto iEnd = iBegin + iReceiveCounts[i];
            std::copy(iBegin, iEnd, std::back_inserter(rangeElem.iBeginEnd));

            ++i;
            return rangeElem;
          });

        return indexRanges;
      }

      // Calculate the intersection between two index ranges.
      indexRange getIntersection(const indexRange& indexRangeA,
        const indexRange& indexRangeB) {

        // Declare result.
        auto intersection = indexRange{};

        // get j intersection range.
        intersection.jBeginEnd = idxPair(
          std::max(indexRangeA.jBeginEnd.first, indexRangeB.jBeginEnd.first),
          std::min(indexRangeA.jBeginEnd.second, indexRangeB.jBeginEnd.second));

        // get i intersection range.
        if (intersection.jBeginEnd.first < intersection.jBeginEnd.second) {

          // get iterators.
          auto iBeginA = indexRangeA.iBeginEnd.begin()
            + intersection.jBeginEnd.first - indexRangeA.jBeginEnd.first;
          auto iBeginB = indexRangeB.iBeginEnd.begin()
            + intersection.jBeginEnd.first - indexRangeB.jBeginEnd.first;
          auto iEndA = indexRangeA.iBeginEnd.end()
            + intersection.jBeginEnd.second - indexRangeA.jBeginEnd.second;

          std::transform(
            iBeginA, iEndA, iBeginB, std::back_inserter(intersection.iBeginEnd),
            [](const idxPair& iElemA, const idxPair& iElemB){
              return idxPair(std::max(iElemA.first, iElemB.first),
                std::min(iElemA.second, iElemB.second));
            });
        }
        return intersection;
      }

      // Calculate the intersections between an index range and a vector of
      // index ranges.
      std::vector<indexRange> getIntersections(const indexRange& indexRangeA,
        const std::vector<indexRange>& indexRangesB) {

        // Declare result.
        auto intersections = std::vector<indexRange>{};

        // Get MPI comm size.
        auto mpiSize = mpi::comm().size();

        // Calculate intersection between index ranges.
        std::generate_n(std::back_inserter(intersections), mpiSize,
          [&, i = size_t{0}]() mutable {
            return getIntersection(indexRangeA, indexRangesB[i++]);
          });

        return intersections;
      }

      // Get the count and displacement arrays needed for an allToAll MPI
      // communication.
      std::pair<std::vector<int>, std::vector<int>>
        getCountsDisplacements(const std::vector<indexRange>& indexRanges) {

        // Declare result.
        auto countsDisplacements =
          std::pair<std::vector<int>, std::vector<int>>{};

        // Count the number of elements in each index range.
        std::transform(indexRanges.begin(), indexRanges.end(),
          std::back_inserter(countsDisplacements.first),
          [](const indexRange& rangeElem){

            // Accumulate size of positive i range.
            int count = std::accumulate(rangeElem.iBeginEnd.begin(),
              rangeElem.iBeginEnd.end(), 0,
              [](const int& cumulant, const idxPair& iElem) -> int {
                return cumulant + std::max(iElem.second - iElem.first, 0);
              });

            return count;
          });

        // Calculate displacements from counts.
        countsDisplacements.second.push_back(0);
        std::partial_sum(countsDisplacements.first.begin(),
          countsDisplacements.first.end() - 1,
          std::back_inserter(countsDisplacements.second));

        return countsDisplacements;
      }

      // Produce a trimmed vector of index ranges with no zero-sized ranges.
      std::vector<indexRange> trimIntersections(
        const std::vector<indexRange>& intersections,
        const std::vector<int>& counts) {

        // Declare result.
        auto trimmedIntersections = std::vector<indexRange>{};

        // Remove intersections with no elements.
        std::for_each(intersections.begin(), intersections.end(),
          [&, i = size_t{0}](const indexRange& intersection) mutable {
            if (counts[i++]) trimmedIntersections.push_back(intersection);
            return;
          });

        return trimmedIntersections;
      }

    }
  }
}
