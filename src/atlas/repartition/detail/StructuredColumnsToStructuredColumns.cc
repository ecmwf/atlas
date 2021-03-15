/*
 * (C) Crown Copyright 2021 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <algorithm>
#include <numeric>

#include "atlas/array/MakeView.h"
#include "atlas/field.h"
#include "atlas/field/FieldSet.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/repartition/detail/StructuredColumnsToStructuredColumns.h"
#include "atlas/repartition/detail/RepartitionUtils.h"

namespace atlas {
  namespace repartition {
    namespace detail {

      // Anonymous name space for helper functions.
      namespace {

        // Transform a vector element-by-element with a functor.
        template <typename outType, typename inType, typename functorType>
        std::vector<outType> transformVector(
          const std::vector<inType>& inVector, functorType functor) {

          // Declare result.
          auto outVector = std::vector<outType>{};
          outVector.reserve(inVector.size());

          // Transform vector.
          std::transform(inVector.begin(), inVector.end(),
            std::back_inserter(outVector), functor);

          return outVector;
        }

        // Visit each index in a FuncSpaceRangeVector.
        template <typename functorType>
        void forEachIndex(const FuncSpaceRangeVector& ranges,
          functorType functor) {

          // Loop over all ranges.
          std::for_each(ranges.cbegin(), ranges.cend(),
            [&](const FuncSpaceRange& range) {

            range.forEach(functor);
            return;
          });
          return;
        }

      }

      //========================================================================
      // Class public methods implementation.
      //========================================================================

      // Constructor.
      StructuredColumnsToStructuredColumns::
        StructuredColumnsToStructuredColumns(
          const FunctionSpace& sourceFunctionSpace,
          const FunctionSpace& targetFunctionSpace) :
          RepartitionImpl(sourceFunctionSpace, targetFunctionSpace),
          sourceStructuredColumnsPtr_(
            getSourceFunctionSpace()->cast<StructuredColumns>()),
          targetStructuredColumnsPtr_(
            getTargetFunctionSpace()->cast<StructuredColumns>()) {

        // Check casts.
        TRY_CAST(StructuredColumns, sourceStructuredColumnsPtr_);
        TRY_CAST(StructuredColumns, targetStructuredColumnsPtr_);

        // Check that grids match.
        CHECK_GRIDS(StructuredColumns,
          sourceStructuredColumnsPtr_, targetStructuredColumnsPtr_);


        // Get source and target range of this function space.
        const auto sourceRange = FuncSpaceRange(sourceStructuredColumnsPtr_);
        const auto targetRange = FuncSpaceRange(targetStructuredColumnsPtr_);


        // Get source and target ranges over all PEs.
        const auto sourceRanges = sourceRange.getFuncSpaceRanges();
        const auto targetRanges = targetRange.getFuncSpaceRanges();


        // Get intersections between sourceRange and targetRanges.
        auto getIntersections = [](const FuncSpaceRange& range,
          const FuncSpaceRangeVector& ranges) {

          return transformVector<FuncSpaceRange>(ranges,
            [&](const FuncSpaceRange& rangesElem) {

            return range & rangesElem;
          });
        };

        sendIntersections_ = getIntersections(sourceRange, targetRanges);
        recvIntersections_ = getIntersections(targetRange, sourceRanges);


        // Get counts and displacements for MPI communication.
        auto getCountsDisplacements = [](
          const FuncSpaceRangeVector& intersections, const idx_t levels) {

          const auto counts = transformVector<int>(intersections,
            [&](const FuncSpaceRange& intersection) {

            return static_cast<int>(intersection.getElemCount() * levels);
          });

          auto displacements = std::vector<int>{0};
          std::partial_sum(counts.cbegin(), counts.cend() - 1,
            std::back_inserter(displacements));

          return std::make_pair(counts, displacements);
        };

        std::tie(sendCounts_, sendDisplacements_) =
          getCountsDisplacements(sendIntersections_,
          sourceStructuredColumnsPtr_->levels());

        std::tie(recvCounts_, recvDisplacements_) =
          getCountsDisplacements(recvIntersections_,
          targetStructuredColumnsPtr_->levels());


        // Trim off invalid intersections.
        auto trimIntersections = [](FuncSpaceRangeVector& intersections) {

          intersections.erase(
            std::remove_if(intersections.begin(), intersections.end(),
              [](const FuncSpaceRange& intersection) {

              return !(intersection.getElemCount());
            }), intersections.end());

          return;
        };

        trimIntersections(sendIntersections_);
        trimIntersections(recvIntersections_);

        return;
      }

      void StructuredColumnsToStructuredColumns::execute(
        const Field& sourceField, Field& targetField) const {

        // Check functionspace casts.
        TRY_CAST(StructuredColumns, sourceField.functionspace().get());
        TRY_CAST(StructuredColumns, targetField.functionspace().get());

        // Check source grids match.
        CHECK_GRIDS(StructuredColumns, sourceField.functionspace().get(),
          sourceStructuredColumnsPtr_);

        // Check target grids match.
        CHECK_GRIDS(StructuredColumns, targetField.functionspace().get(),
          targetStructuredColumnsPtr_);

        // Check data types match.
        CHECK_FIELD_DATA_TYPE(sourceField, targetField);

        // Determine data type of field and execute.
        switch(sourceField.datatype().kind()){

          case array::DataType::KIND_REAL64 :
            doExecute<double>(sourceField, targetField);
            break;

          case array::DataType::KIND_REAL32 :
            doExecute<float>(sourceField, targetField);
            break;

          case array::DataType::KIND_INT32 :
            doExecute<int>(sourceField, targetField);
            break;

          case array::DataType::KIND_INT64 :
            doExecute<long>(sourceField, targetField);
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
        CHECK_FIELD_SET_SIZE(sourceFieldSet, targetFieldSet);

        auto targetFieldSetIt = targetFieldSet.begin();
        std::for_each(sourceFieldSet.cbegin(), sourceFieldSet.cend(),
          [&](const Field& sourceField) {

          execute(sourceField, *targetFieldSetIt++);
          return;
        });

        return;
      }

      //========================================================================
      // Class private methods implementation.
      //========================================================================

      template<typename fieldType>
      void StructuredColumnsToStructuredColumns::doExecute(
        const Field& sourceField, Field& targetField) const {


        // Make Atlas view objects.
        const auto sourceView = array::make_view<fieldType, 2>(sourceField);
        auto targetView = array::make_view<fieldType, 2>(targetField);


        // Get buffer sizes.
        const auto nSend = (sendCounts_.back() + sendDisplacements_.back()) *
          sourceStructuredColumnsPtr_->levels();
        const auto nRecv = (recvCounts_.back() + recvDisplacements_.back()) *
          targetStructuredColumnsPtr_->levels();


        // Allocate send and receive buffers.
        auto sendBuffer = std::vector<fieldType>(static_cast<size_t>(nSend));
        auto recvBuffer = std::vector<fieldType>(static_cast<size_t>(nRecv));


        // Set send functor.
        auto sendBufferIt = sendBuffer.begin();
        auto sendFunctor = [&] (const idx_t i, const idx_t j) {

          // Loop over levels
          auto iNode = sourceStructuredColumnsPtr_->index(i, j);
          auto kEnd = sourceStructuredColumnsPtr_->levels();
          for (idx_t k = 0; k < kEnd; ++k) {

            *sendBufferIt++ = sourceView(iNode, k);
          }
          return;
        };


        // Set receive functor.
        auto recvBufferIt = recvBuffer.cbegin();
        auto recvFunctor = [&](const idx_t i, const idx_t j) {

          // Loop over levels
          auto iNode = targetStructuredColumnsPtr_->index(i, j);
          auto kEnd = targetStructuredColumnsPtr_->levels();
          for (idx_t k = 0; k < kEnd; ++k) {

            targetView(iNode, k) = *recvBufferIt++;
          }
          return;
        };


        // Write data to buffer.
        forEachIndex(sendIntersections_, sendFunctor);

        // Perform allToAllv.
        mpi::comm().allToAllv(
          sendBuffer.data(), sendCounts_.data(), sendDisplacements_.data(),
          recvBuffer.data(), recvCounts_.data(), recvDisplacements_.data());

        // Read data from buffer.
        forEachIndex(recvIntersections_, recvFunctor);

        return;
      }

      //========================================================================
      // Index range methods.
      //========================================================================

      // Constructor.
      FuncSpaceRange::FuncSpaceRange(
          const StructuredColumns* const structuredColumnsPtr) {

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
        const auto mpiSize = static_cast<size_t>(atlas::mpi::comm().size());

        // Set recv buffer for j range.
        auto jRecvBuffer = idxPairVector(mpiSize);

        // Perform all gather.
        atlas::mpi::comm().allGather(
          jBeginEnd_, jRecvBuffer.begin(), jRecvBuffer.end());

        // Set i receive counts.
        auto iRecvCounts = transformVector<int>(jRecvBuffer,
          [](const idxPair& jElem) {

          return static_cast<int>(jElem.second - jElem.first);
        });

        // Set recv displacements for i range.
        auto iRecvDisplacements = std::vector<int>{0};
        std::partial_sum(iRecvCounts.cbegin(), iRecvCounts.cend() - 1,
          std::back_inserter(iRecvDisplacements));

        // Set recv buffer for i range.
        auto irecvBuffer = idxPairVector(static_cast<size_t>(
          iRecvDisplacements.back() + iRecvCounts.back()));

        // Perform all gather.
        atlas::mpi::comm().allGatherv(
          iBeginEnd_.cbegin(), iBeginEnd_.cend(), irecvBuffer.begin(),
          iRecvCounts.data(), iRecvDisplacements.data());

        // Make vector of indexRange structs.
        auto indexRanges = FuncSpaceRangeVector{};
        for (size_t i = 0; i < mpiSize; ++i) {

          auto indexRange = FuncSpaceRange{};
          indexRange.jBeginEnd_ = jRecvBuffer[i];
          const auto iBegin = irecvBuffer.cbegin() + iRecvDisplacements[i];
          const auto iEnd = iBegin + iRecvCounts[i];
          std::copy(iBegin, iEnd, std::back_inserter(indexRange.iBeginEnd_));

          indexRanges.push_back(indexRange);
        }

        return indexRanges;
      }

      // Count number of elements in index range.
      idx_t FuncSpaceRange::getElemCount() const {

        // Accumulate size of positive i range.
        const auto count =
          std::accumulate(iBeginEnd_.cbegin(), iBeginEnd_.cend(),
          0, [](const int cumulant, const idxPair iElem) {

          // Only count positive differences.
          return cumulant +
            static_cast<int>(std::max(iElem.second - iElem.first, 0));
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

          std::transform(iBeginA, iEndA, iBeginB,
            std::back_inserter(intersection.iBeginEnd_),
            [](const idxPair iElemA, const idxPair iElemB) {

            return std::make_pair(std::max(iElemA.first, iElemB.first),
              std::min(iElemA.second, iElemB.second));
          });
        }
        return intersection;
      }

      // Loop over all indices. Functor should have signature
      // functor(const idx_t i, const idx_t j).
      template <typename functorType>
      void FuncSpaceRange::forEach(functorType functor) const {

        auto iBeginEndIt = iBeginEnd_.begin();

        // Loop over j.
        for (auto j = jBeginEnd_.first; j < jBeginEnd_.second; ++j) {

          const auto iBeginEnd = *iBeginEndIt++;

          // Loop over i.
          for (auto i = iBeginEnd.first; i < iBeginEnd.second; ++i) {

            // Call functor.
            functor(i, j);

          }
        }

        return;
      }

    }
  }
}
