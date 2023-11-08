/*
 * (C) Crown Copyright 2021 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "RedistributeStructuredColumns.h"

#include <algorithm>
#include <numeric>

#include "atlas/array/MakeView.h"
#include "atlas/field.h"
#include "atlas/field/FieldSet.h"
#include "atlas/functionspace/StructuredColumns.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/runtime/Exception.h"

namespace atlas {
namespace redistribution {
namespace detail {

// Anonymous name space for helper functions.
namespace {

// Transform a vector element-by-element with a functor.
template <typename outType, typename inType, typename functorType>
std::vector<outType> transformVector(const std::vector<inType>& inVector, const functorType& functor) {
    // Declare result.
    auto outVector = std::vector<outType>{};
    outVector.reserve(inVector.size());

    // Transform vector.
    std::transform(inVector.begin(), inVector.end(), std::back_inserter(outVector), functor);

    return outVector;
}

// Visit each index in a StructuredIndexRangeVector.
template <typename functorType>
void forEachIndex(const StructuredIndexRangeVector& ranges, const functorType& functor) {
    // Loop over all ranges.
    std::for_each(ranges.cbegin(), ranges.cend(), [&](const StructuredIndexRange& range) {
        range.forEach(functor);
        return;
    });
    return;
}

}  // namespace

//========================================================================
// Class public methods implementation.
//========================================================================

// Constructor.
void RedistributeStructuredColumns::do_setup() {
    source_ = source();
    target_ = target();

    // Assert that functionspaces are StructuredColumns.
    ATLAS_ASSERT(source_, source().type() + " must be StructuredColumns");
    ATLAS_ASSERT(target_, target().type() + " must be StructuredColumns");

    // Check that grids match.
    ATLAS_ASSERT(source_.grid().name() == target_.grid().name());

    // Check levels match.
    ATLAS_ASSERT(source_.levels() == target_.levels());

    // Check that communicators match.
    ATLAS_ASSERT(source_.mpi_comm() == target_.mpi_comm());
    mpi_comm_ = source_.mpi_comm();


    // Get source and target range of this function space.
    const auto sourceRange = StructuredIndexRange(source_);
    const auto targetRange = StructuredIndexRange(target_);


    // Get source and target ranges over all PEs.
    const auto sourceRanges = sourceRange.getStructuredIndexRanges();
    const auto targetRanges = targetRange.getStructuredIndexRanges();


    // Get intersections between sourceRange and targetRanges.
    auto getIntersections = [](const StructuredIndexRange& range, const StructuredIndexRangeVector& ranges) {
        return transformVector<StructuredIndexRange>(
            ranges, [&](const StructuredIndexRange& rangesElem) { return range & rangesElem; });
    };

    sendIntersections_ = getIntersections(sourceRange, targetRanges);
    recvIntersections_ = getIntersections(targetRange, sourceRanges);


    // Get counts and displacements for MPI communication.
    auto getCountsDisplacements = [](const StructuredIndexRangeVector& intersections, const idx_t levels) {
        const auto counts = transformVector<int>(intersections, [&](const StructuredIndexRange& intersection) {
            return static_cast<int>(intersection.getElemCount() * levels);
        });

        auto displacements = std::vector<int>{0};
        std::partial_sum(counts.cbegin(), counts.cend() - 1, std::back_inserter(displacements));

        return std::make_pair(counts, displacements);
    };

    std::tie(sendCounts_, sendDisplacements_) = getCountsDisplacements(sendIntersections_, source_.levels());

    std::tie(recvCounts_, recvDisplacements_) = getCountsDisplacements(recvIntersections_, target_.levels());


    // Trim off invalid intersections.
    auto trimIntersections = [](StructuredIndexRangeVector& intersections) {
        intersections.erase(
            std::remove_if(intersections.begin(), intersections.end(),
                           [](const StructuredIndexRange& intersection) { return !(intersection.getElemCount()); }),
            intersections.end());

        return;
    };

    trimIntersections(sendIntersections_);
    trimIntersections(recvIntersections_);

    return;
}

void RedistributeStructuredColumns::execute(const field::FieldImpl* sf, field::FieldImpl* tf) const {
    const field::FieldImpl& sourceField = *sf;
    field::FieldImpl& targetField = *tf;
    // Assert that fields are defined on StructuredColumns.
    ATLAS_ASSERT(functionspace::StructuredColumns(sourceField.functionspace()));
    ATLAS_ASSERT(functionspace::StructuredColumns(targetField.functionspace()));

    // Check that grids match.
    ATLAS_ASSERT(functionspace::StructuredColumns(sourceField.functionspace()).grid().name() == source_.grid().name());
    ATLAS_ASSERT(functionspace::StructuredColumns(targetField.functionspace()).grid().name() == target_.grid().name());

    // Check levels match.
    ATLAS_ASSERT(sourceField.levels() == source_.levels());
    ATLAS_ASSERT(targetField.levels() == target_.levels());

    // Determine data type of field and execute.
    switch (sourceField.datatype().kind()) {
        case array::DataType::KIND_REAL64:
            do_execute<double>(sourceField, targetField);
            break;

        case array::DataType::KIND_REAL32:
            do_execute<float>(sourceField, targetField);
            break;

        case array::DataType::KIND_INT32:
            do_execute<int>(sourceField, targetField);
            break;

        case array::DataType::KIND_INT64:
            do_execute<long>(sourceField, targetField);
            break;

        default:
            throw_NotImplemented("No implementation for data type " + sourceField.datatype().str(), Here());
    }

    return;
}

void RedistributeStructuredColumns::execute(const field::FieldSetImpl* sourceFieldSet, field::FieldSetImpl* targetFieldSet) const {
    // Check that both FieldSets are the same size.
    ATLAS_ASSERT(sourceFieldSet->size() == targetFieldSet->size());

    auto targetFieldSetIt = targetFieldSet->begin();
    for (int i = 0; i < sourceFieldSet->size(); i++) {
        execute(&sourceFieldSet[i], &targetFieldSet[i]);
    }
    //std::for_each(sourceFieldSet->cbegin(), sourceFieldSet->cend(), [&](const field::FieldImpl* sourceField) {
    //    execute(sourceField, targetFieldSetIt++);
    //    return;
    //});

    return;
}

//========================================================================
// Class private methods implementation.
//========================================================================

template <typename fieldType>
void RedistributeStructuredColumns::do_execute(const field::FieldImpl& sourceField, field::FieldImpl& targetField) const {
    // Make Atlas view objects.
    const auto sourceView = array::make_view<fieldType, 2>(sourceField);
    auto targetView       = array::make_view<fieldType, 2>(targetField);


    // Get buffer sizes.
    const auto nSend = sendCounts_.back() + sendDisplacements_.back();
    const auto nRecv = recvCounts_.back() + recvDisplacements_.back();


    // Allocate send and receive buffers.
    auto sendBuffer = std::vector<fieldType>(static_cast<size_t>(nSend));
    auto recvBuffer = std::vector<fieldType>(static_cast<size_t>(nRecv));


    // Set send functor.
    auto sendBufferIt = sendBuffer.begin();
    auto sendFunctor  = [&](const idx_t i, const idx_t j) {
        // Loop over levels
        const auto iNode = source_.index(i, j);
        const auto kEnd  = source_.levels();
        for (idx_t k = 0; k < kEnd; ++k) {
            *sendBufferIt++ = sourceView(iNode, k);
        }
        return;
    };


    // Set receive functor.
    auto recvBufferIt = recvBuffer.cbegin();
    auto recvFunctor  = [&](const idx_t i, const idx_t j) {
        // Loop over levels
        const auto iNode = target_.index(i, j);
        const auto kEnd  = target_.levels();
        for (idx_t k = 0; k < kEnd; ++k) {
            targetView(iNode, k) = *recvBufferIt++;
        }
        return;
    };


    // Write data to buffer.
    forEachIndex(sendIntersections_, sendFunctor);

    // Communicate.
    mpi::comm(mpi_comm_).allToAllv(sendBuffer.data(), sendCounts_.data(), sendDisplacements_.data(), recvBuffer.data(),
                                   recvCounts_.data(), recvDisplacements_.data());

    // Read data from buffer.
    forEachIndex(recvIntersections_, recvFunctor);

    return;
}

//========================================================================
// Index range methods.
//========================================================================

// Constructor.
StructuredIndexRange::StructuredIndexRange(const functionspace::StructuredColumns& structuredColumns) {
    jBeginEnd_ = std::make_pair(structuredColumns.j_begin(), structuredColumns.j_end());

    for (auto j = jBeginEnd_.first; j < jBeginEnd_.second; ++j) {
        iBeginEnd_.push_back(std::make_pair(structuredColumns.i_begin(j), structuredColumns.i_end(j)));
    }

    mpi_comm_ = structuredColumns.mpi_comm();

    return;
}

// Get index ranges from all PEs.
StructuredIndexRangeVector StructuredIndexRange::getStructuredIndexRanges() const {
    auto& comm = mpi::comm(mpi_comm());

    // Get MPI communicator size.
    const auto mpiSize = static_cast<size_t>(comm.size());

    // Set recv buffer for j range.
    auto jRecvBuffer = idxPairVector(mpiSize);

    // Perform all gather.
    comm.allGather(jBeginEnd_, jRecvBuffer.begin(), jRecvBuffer.end());

    // Set i receive counts.
    auto iRecvCounts = transformVector<int>(
        jRecvBuffer, [](const idxPair& jElem) { return static_cast<int>(jElem.second - jElem.first); });

    // Set recv displacements for i range.
    auto iRecvDisplacements = std::vector<int>{0};
    std::partial_sum(iRecvCounts.cbegin(), iRecvCounts.cend() - 1, std::back_inserter(iRecvDisplacements));

    // Set recv buffer for i range.
    auto irecvBuffer = idxPairVector(static_cast<size_t>(iRecvDisplacements.back() + iRecvCounts.back()));

    // Perform all gather.
    comm.allGatherv(iBeginEnd_.cbegin(), iBeginEnd_.cend(), irecvBuffer.begin(), iRecvCounts.data(),
                    iRecvDisplacements.data());

    // Make vector of indexRange structs.
    auto indexRanges = StructuredIndexRangeVector{};
    for (size_t i = 0; i < mpiSize; ++i) {
        auto indexRange       = StructuredIndexRange{};
        indexRange.jBeginEnd_ = jRecvBuffer[i];
        const auto iBegin     = irecvBuffer.cbegin() + iRecvDisplacements[i];
        const auto iEnd       = iBegin + iRecvCounts[i];
        std::copy(iBegin, iEnd, std::back_inserter(indexRange.iBeginEnd_));

        indexRanges.push_back(indexRange);
    }

    return indexRanges;
}

// Count number of elements in index range.
idx_t StructuredIndexRange::getElemCount() const {
    // Accumulate size of positive i range.
    const auto count =
        std::accumulate(iBeginEnd_.cbegin(), iBeginEnd_.cend(), 0, [](const idx_t cumulant, const idxPair iElem) {
            // Only count positive differences.
            return cumulant + static_cast<idx_t>(std::max(iElem.second - iElem.first, idx_t(0)));
        });

    return count;
}

// Return the intersection between two index ranges.
StructuredIndexRange StructuredIndexRange::operator&(const StructuredIndexRange& indexRange) const {
    // Declare result.
    auto intersection = StructuredIndexRange{};

    // get j intersection range.
    intersection.jBeginEnd_ = std::make_pair(std::max(jBeginEnd_.first, indexRange.jBeginEnd_.first),
                                             std::min(jBeginEnd_.second, indexRange.jBeginEnd_.second));

    // get i intersection range.
    if (intersection.jBeginEnd_.first < intersection.jBeginEnd_.second) {
        // get iterators.
        const auto iBeginA = iBeginEnd_.cbegin() + intersection.jBeginEnd_.first - jBeginEnd_.first;
        const auto iBeginB =
            indexRange.iBeginEnd_.cbegin() + intersection.jBeginEnd_.first - indexRange.jBeginEnd_.first;
        const auto iEndA = iBeginEnd_.cend() + intersection.jBeginEnd_.second - jBeginEnd_.second;

        std::transform(iBeginA, iEndA, iBeginB, std::back_inserter(intersection.iBeginEnd_),
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
void StructuredIndexRange::forEach(const functorType& functor) const {
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

namespace {
static RedistributionImplBuilder<RedistributeStructuredColumns> register_builder(
    RedistributeStructuredColumns::static_type());
}

}  // namespace detail
}  // namespace redistribution
}  // namespace atlas
