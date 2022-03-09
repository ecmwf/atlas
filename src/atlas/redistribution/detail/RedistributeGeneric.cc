/*
 * (C) Crown Copyright 2021 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <numeric>
#include <vector>

#include "atlas/field/Field.h"
#include "atlas/functionspace/CellColumns.h"
#include "atlas/functionspace/EdgeColumns.h"
#include "atlas/functionspace/NodeColumns.h"
#include "atlas/functionspace/PointCloud.h"
#include "atlas/functionspace/StructuredColumns.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/redistribution/detail/RedistributeGeneric.h"
#include "atlas/redistribution/detail/RedistributionImplFactory.h"
#include "atlas/util/Unique.h"


namespace atlas {
namespace redistribution {
namespace detail {

using mesh::HybridElements;
using mesh::Nodes;

// Helper type definitions and functions for redistribution.
namespace {

// Define index-UID struct. (Needed to overload "<").
struct IdxUid : public std::pair<idx_t, uidx_t> {
    using std::pair<idx_t, uidx_t>::pair;
};

Field getGhostField(const FunctionSpace& functionspace) {
    if (functionspace::NodeColumns(functionspace) || functionspace::StructuredColumns(functionspace) ||
        functionspace::PointCloud(functionspace)) {
        return functionspace.ghost();
    }
    if (functionspace::EdgeColumns(functionspace) || functionspace::CellColumns(functionspace)) {
        // TODO: Move something like this into the functionspace::EdgeColumns and functionspace::CellColumns

        // Get mesh elements.
        const auto& elems = functionspace::EdgeColumns(functionspace)
                                ? functionspace::EdgeColumns(functionspace).edges()
                                : functionspace::CellColumns(functionspace).cells();


        // Make ghost field.
        auto ghost_field = Field("ghost", array::make_datatype<int>(), array::make_shape(functionspace.size()));

        // Make views.
        auto ghost        = array::make_view<int, 1>(ghost_field);
        auto remote_index = array::make_indexview<idx_t, 1>(elems.remote_index());
        auto partition    = array::make_view<int, 1>(elems.partition());

        // Set ghost field.
        const auto thisPart = static_cast<int>(mpi::comm().rank());
        for (idx_t i = 0; i < ghost.shape(0); ++i) {
            ghost(i) = partition(i) != thisPart || remote_index(i) != i;
        }
        return ghost_field;
    }
    return functionspace.ghost();
}


// Create UID field.
Field getUidField(const FunctionSpace& functionspace) {
    // Note there is a difference between uid_field and global_index in the type: uidx_t vs gidx_t
    // Currently this aliases the same type though...

    if (not functionspace::PointCloud(functionspace)) {
        return functionspace.global_index();
    }
    else {
        // PointCloud has no global_index()... TODO, put something in place there.

        // Make UID field.
        Field uid_field = Field("Unique ID", array::make_datatype<uidx_t>(), array::make_shape(functionspace.size()));

        // Make views.
        auto uid          = array::make_view<uidx_t, 1>(uid_field);
        const auto lonlat = array::make_view<double, 2>(functionspace.lonlat());

        // Set UIDs
        for (idx_t i = 0; i < uid_field.shape(0); ++i) {
            uid(i) = util::unique_lonlat(lonlat(i, LON), lonlat(i, LAT));
        }
        return uid_field;
    }
}

// Resolve functionspace implementation type and get unique ID vector.
std::vector<IdxUid> getUidVec(const FunctionSpace& functionspace) {
    // Get ghost and unique ID fields from functionspace.
    const auto ghost_field = getGhostField(functionspace);
    const auto uid         = getUidField(functionspace);

    // Make views to fields.
    const auto ghost   = array::make_view<int, 1>(ghost_field);
    const auto uidView = array::make_view<uidx_t, 1>(uid);

    auto uidVec = std::vector<IdxUid>{};
    uidVec.reserve(functionspace.size());

    // Get UIDs for non ghost elems.
    for (idx_t i = 0; i < functionspace.size(); ++i) {
        if (not ghost(i)) {
            uidVec.emplace_back(i, uidView(i));
        }
    }

    // Sort by UID value.
    std::sort(uidVec.begin(), uidVec.end(), [](const IdxUid& a, const IdxUid& b) { return a.second < b.second; });

    // Check UIDs are unique.
    if (ATLAS_BUILD_TYPE_DEBUG) {
        const size_t vecSize = uidVec.size();
        std::unique(uidVec.begin(), uidVec.end(),
                    [](const IdxUid& a, const IdxUid& b) { return a.second == b.second; });
        ATLAS_ASSERT(uidVec.size() == vecSize, "Unique ID set has duplicate members");
    }

    return uidVec;
}

// Copy UID index
std::vector<idx_t> getUidIdx(const std::vector<IdxUid>& uidVec) {
    auto idxVec = std::vector<idx_t>{};
    idxVec.reserve(uidVec.size());
    std::transform(uidVec.begin(), uidVec.end(), std::back_inserter(idxVec),
                   [](const IdxUid& uid) { return uid.first; });
    return idxVec;
}

// Copy UID value
std::vector<uidx_t> getUidVal(const std::vector<IdxUid>& uidVec) {
    auto valVec = std::vector<uidx_t>{};
    valVec.reserve(uidVec.size());
    std::transform(uidVec.begin(), uidVec.end(), std::back_inserter(valVec),
                   [](const IdxUid& uid) { return uid.second; });
    return valVec;
}

// Communicate UID values, return receive buffer and displacements.
std::pair<std::vector<uidx_t>, std::vector<int>> communicateUid(const std::vector<uidx_t>& sendBuffer) {
    auto counts = std::vector<int>(mpi::comm().size());
    mpi::comm().allGather(static_cast<int>(sendBuffer.size()), counts.begin(), counts.end());

    auto disps = std::vector<int>{};
    disps.reserve(mpi::comm().size() + 1);
    disps.push_back(0);
    std::partial_sum(counts.begin(), counts.end(), std::back_inserter(disps));


    auto recvBuffer = std::vector<uidx_t>(static_cast<size_t>(disps.back()));

    mpi::comm().allGatherv(sendBuffer.begin(), sendBuffer.end(), recvBuffer.begin(), counts.data(), disps.data());

    return std::make_pair(recvBuffer, disps);
}

// Need to overload "<" operator for comparison of asymmetric types in
// std::set_intersection.
bool operator<(const IdxUid& lhs, const uidx_t& rhs) {
    return lhs.second < rhs;
}
bool operator<(const uidx_t& lhs, const IdxUid& rhs) {
    return lhs < rhs.second;
}

// Find the intersection between local and global UIDs, then return local
// indices of incections and PE dispacements in vector.
std::pair<std::vector<idx_t>, std::vector<int>> getUidIntersection(const std::vector<IdxUid>& localUids,
                                                                   const std::vector<uidx_t>& globalUids,
                                                                   const std::vector<int>& globalDisps) {
    auto uidIntersection = std::vector<IdxUid>{};
    uidIntersection.reserve(localUids.size());

    auto disps = std::vector<int>{};
    disps.reserve(mpi::comm().size() + 1);
    disps.push_back(0);

    // Loop over all PE and find UID intersection.
    for (size_t i = 0; i < mpi::comm().size(); ++i) {
        // Get displaced iterators.
        auto globalUidsBegin = globalUids.begin() + globalDisps[i];
        auto globalUidsEnd   = globalUids.begin() + globalDisps[i + 1];

        // Get intersection.
        std::set_intersection(localUids.begin(), localUids.end(), globalUidsBegin, globalUidsEnd,
                              std::back_inserter(uidIntersection));

        // Set intersection displacement.
        disps.push_back(static_cast<int>(uidIntersection.size()));
    }

    // Check that the set of all intersections matches UIDs on local PE.
    if (ATLAS_BUILD_TYPE_DEBUG) {
        auto tempUids = uidIntersection;
        std::sort(tempUids.begin(), tempUids.end(),
                  [](const IdxUid& a, const IdxUid& b) { return a.second < b.second; });
        ATLAS_ASSERT(tempUids == localUids, "Set of all UID intersections does not match local UIDs.");
    }

    // Return local indices of intersection and displacements.
    return make_pair(getUidIdx(uidIntersection), disps);
}


// Iterate over a field, in the order of an index list, and apply a functor to
// each element.

// Recursive ForEach to visit all elements of field.
template <int Rank, int Dim = 0>
struct ForEach {
    template <typename Value, typename Functor, typename... Idxs>
    static void apply(const std::vector<idx_t>& idxList, array::ArrayView<Value, Rank>& fieldView, const Functor& f,
                      Idxs... idxs) {
        // Iterate over dimension Dim of array.
        for (idx_t idx = 0; idx < fieldView.shape(Dim); ++idx) {
            ForEach<Rank, Dim + 1>::apply(idxList, fieldView, f, idxs..., idx);
        }
    }
};

// Beginning of recursion when Dim == 0.
template <int Rank>
struct ForEach<Rank, 0> {
    template <typename Value, typename Functor, typename... Idxs>
    static void apply(const std::vector<idx_t>& idxList, array::ArrayView<Value, Rank>& fieldView, const Functor& f,
                      Idxs... idxs) {
        // Iterate over dimension 0 of array in order defined by idxList.
        for (idx_t idx : idxList) {
            ForEach<Rank, 1>::apply(idxList, fieldView, f, idxs..., idx);
        }
    }
};

// End of recursion when Dim == Rank.
template <int Rank>
struct ForEach<Rank, Rank> {
    template <typename Value, typename Functor, typename... Idxs>
    static void apply(const std::vector<idx_t>& idxList, array::ArrayView<Value, Rank>& fieldView, const Functor& f,
                      Idxs... idxs) {
        // Apply functor.
        f(fieldView(idxs...));
    }
};

}  // namespace

void RedistributeGeneric::do_setup() {
    // get a unique ID (UID) for each owned member of functionspace.
    const auto sourceUidVec = getUidVec(source());
    const auto targetUidVec = getUidVec(target());

    // Communicate UID vectors to all PEs.
    auto sourceGlobalUids                         = std::vector<uidx_t>{};
    auto sourceGlobalDisps                        = std::vector<int>{};
    std::tie(sourceGlobalUids, sourceGlobalDisps) = communicateUid(getUidVal(sourceUidVec));
    auto targetGlobalUids                         = std::vector<uidx_t>{};
    auto targetGlobalDisps                        = std::vector<int>{};
    std::tie(targetGlobalUids, targetGlobalDisps) = communicateUid(getUidVal(targetUidVec));

    // Get intersection of local UIDs and Global UIDs.
    std::tie(sourceLocalIdx_, sourceDisps_) = getUidIntersection(sourceUidVec, targetGlobalUids, targetGlobalDisps);
    std::tie(targetLocalIdx_, targetDisps_) = getUidIntersection(targetUidVec, sourceGlobalUids, sourceGlobalDisps);
}

void RedistributeGeneric::execute(const Field& sourceField, Field& targetField) const {
    //Check functionspaces match.
    ATLAS_ASSERT(sourceField.functionspace().type() == source().type());
    ATLAS_ASSERT(targetField.functionspace().type() == target().type());

    // Check Field datatypes match.
    ATLAS_ASSERT(sourceField.datatype() == targetField.datatype());

    // Check Field ranks match.
    ATLAS_ASSERT(sourceField.rank() == targetField.rank());

    // Check number of levels and variables match.
    for (idx_t i = 1; i < sourceField.rank(); ++i) {
        ATLAS_ASSERT(sourceField.shape(i) == targetField.shape(i));
    }

    // Perform redistribution.
    do_execute(sourceField, targetField);
}

void RedistributeGeneric::execute(const FieldSet& sourceFieldSet, FieldSet& targetFieldSet) const {
    // Check field set sizes match.
    ATLAS_ASSERT(sourceFieldSet.size() == targetFieldSet.size());

    // Redistribute fields.
    for (idx_t i = 0; i < sourceFieldSet.size(); ++i) {
        execute(sourceFieldSet[i], targetFieldSet[i]);
    }
}

// Determine datatype.
void RedistributeGeneric::do_execute(const Field& sourceField, Field& targetField) const {
    // Available datatypes defined in array/LocalView.cc
    switch (sourceField.datatype().kind()) {
        case array::DataType::KIND_REAL64: {
            return do_execute<double>(sourceField, targetField);
        }
        case array::DataType::KIND_REAL32: {
            return do_execute<float>(sourceField, targetField);
        }
        case array::DataType::KIND_INT64: {
            return do_execute<long>(sourceField, targetField);
        }
        case array::DataType::KIND_INT32: {
            return do_execute<int>(sourceField, targetField);
        }
        default: {
            ATLAS_THROW_EXCEPTION("No implementation for data type " + sourceField.datatype().str());
        }
    }
}

// Determine rank.
template <typename Value>
void RedistributeGeneric::do_execute(const Field& sourceField, Field& targetField) const {
    // Available ranks defined in array/LocalView.cc
    switch (sourceField.rank()) {
        case 1: {
            return do_execute<Value, 1>(sourceField, targetField);
        }
        case 2: {
            return do_execute<Value, 2>(sourceField, targetField);
        }
        case 3: {
            return do_execute<Value, 3>(sourceField, targetField);
        }
        case 4: {
            return do_execute<Value, 4>(sourceField, targetField);
        }
        case 5: {
            return do_execute<Value, 5>(sourceField, targetField);
        }
        case 6: {
            return do_execute<Value, 6>(sourceField, targetField);
        }
        case 7: {
            return do_execute<Value, 7>(sourceField, targetField);
        }
        case 8: {
            return do_execute<Value, 8>(sourceField, targetField);
        }
        case 9: {
            return do_execute<Value, 9>(sourceField, targetField);
        }
        default: {
            ATLAS_THROW_EXCEPTION("No implementation for rank " + std::to_string(sourceField.rank()));
        }
    }
}

// Perform redistribution.
template <typename Value, int Rank>
void RedistributeGeneric::do_execute(const Field& sourceField, Field& targetField) const {
    // Get array views.
    auto sourceView = array::make_view<Value, Rank>(sourceField);
    auto targetView = array::make_view<Value, Rank>(targetField);

    // Get number of elems per column.
    int elemsPerCol = 1;
    for (int i = 1; i < Rank; ++i) {
        elemsPerCol *= sourceView.shape(i);
    }

    // Set send displacement and counts vectors.
    auto sendDisps = std::vector<int>{};
    sendDisps.reserve(mpi::comm().size() + 1);
    auto sendCounts = std::vector<int>{};
    sendCounts.reserve(mpi::comm().size());
    std::transform(sourceDisps_.begin(), sourceDisps_.end(), std::back_inserter(sendDisps),
                   [&](const int& disp) { return disp * elemsPerCol; });
    std::adjacent_difference(sendDisps.begin() + 1, sendDisps.end(), std::back_inserter(sendCounts));

    // Set recv displacement and counts vectors.
    auto recvDisps = std::vector<int>{};
    recvDisps.reserve(mpi::comm().size() + 1);
    auto recvCounts = std::vector<int>{};
    recvCounts.reserve(mpi::comm().size());
    std::transform(targetDisps_.begin(), targetDisps_.end(), std::back_inserter(recvDisps),
                   [&](const int& disp) { return disp * elemsPerCol; });
    std::adjacent_difference(recvDisps.begin() + 1, recvDisps.end(), std::back_inserter(recvCounts));

    // Allocate send and recv buffers.
    auto sendBuffer   = std::vector<Value>(static_cast<size_t>(sendDisps.back()));
    auto recvBuffer   = std::vector<Value>(static_cast<size_t>(recvDisps.back()));
    auto sendBufferIt = sendBuffer.begin();
    auto recvBufferIt = recvBuffer.cbegin();

    // Copy sourceField to sendBuffer.
    ForEach<Rank>::apply(sourceLocalIdx_, sourceView, [&](const Value& elem) { *sendBufferIt++ = elem; });

    // Perform MPI communication.
    mpi::comm().allToAllv(sendBuffer.data(), sendCounts.data(), sendDisps.data(), recvBuffer.data(), recvCounts.data(),
                          recvDisps.data());

    // Copy recvBuffer to targetField.
    ForEach<Rank>::apply(targetLocalIdx_, targetView, [&](Value& elem) { elem = *recvBufferIt++; });
}

namespace {
static RedistributionImplBuilder<RedistributeGeneric> register_builder(RedistributeGeneric::static_type());
}

}  // namespace detail
}  // namespace redistribution
}  // namespace atlas
