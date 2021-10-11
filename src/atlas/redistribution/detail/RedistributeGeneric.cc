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
#include "atlas/redistribution/detail/RedistributionImplFactory.h"
#include "atlas/redistribution/detail/RedistributeGeneric.h"
#include "atlas/util/Unique.h"


namespace atlas {
namespace redistribution {
namespace detail {

using functionspace::detail::CellColumns;
using functionspace::detail::EdgeColumns;
using functionspace::detail::NodeColumns;
using functionspace::detail::PointCloud;
using functionspace::detail::StructuredColumns;

using mesh::HybridElements;
using mesh::Nodes;

// Helper type definitions and functions for redistribution.
namespace  {

// Helper types to distinguish mesh-based and meshless functionspaces.
template <typename FunctionSpaceType>
using MeshBased = typename std::enable_if<
    std::is_same<FunctionSpaceType, CellColumns>::value ||
    std::is_same<FunctionSpaceType, EdgeColumns>::value ||
    std::is_same<FunctionSpaceType, NodeColumns>::value>::type*;

template <typename FunctionSpaceType>
using Meshless = typename std::enable_if<
    std::is_same<FunctionSpaceType, StructuredColumns>::value ||
    std::is_same<FunctionSpaceType, PointCloud>::value>::type*;

// Define index-UID struct. (Needed to overload "<").
struct IdxUid : public std::pair<idx_t, uidx_t> {
    using std::pair<idx_t, uidx_t>::pair;
};

// Get elements from mesh.
const HybridElements& getElems( const CellColumns* functionSpaceImpl ) {
    return functionSpaceImpl->cells();
}
const HybridElements& getElems( const EdgeColumns* functionSpaceImpl ) {
    return functionSpaceImpl->edges();
}
const Nodes& getElems( const NodeColumns* functionSpaceImpl ) {
    return functionSpaceImpl->nodes();
}

// Create ghost field for mesh based functionspace.
template <typename FunctionSpaceType, MeshBased<FunctionSpaceType> = nullptr>
Field getGhostField( const FunctionSpaceType* functionSpaceImpl ) {

    // Get mesh elements.
    const auto& elems = getElems( functionSpaceImpl );

    // Make ghost field.
    auto ghost = Field( "ghost" , array::make_datatype<int>(),
                        array::make_shape( functionSpaceImpl->size() ) );

    // Make views.
    auto ghostView = array::make_view<int, 1>( ghost );
    auto remoteIndexView = array::make_indexview<idx_t, 1>( elems.remote_index() );
    auto partitionView = array::make_view<int, 1>( elems.partition() );

    // Set ghost field.
    const auto thisPart = static_cast<int>( mpi::comm().rank() );
    for ( idx_t i = 0; i < ghost.shape( 0 ); ++i ) {
        ghostView( i ) = partitionView( i ) != thisPart ||
                         remoteIndexView ( i ) != i;
    }
    return ghost;
}

// Get ghost field from meshless functionspace.
template <typename FunctionSpaceType, Meshless<FunctionSpaceType> = nullptr>
Field getGhostField( const FunctionSpaceType* functionSpaceImpl ) {
    return functionSpaceImpl->ghost();
}

// Get UID field from mesh based functionspace.
template <typename FunctionSpaceType, MeshBased<FunctionSpaceType> = nullptr>
Field getUidField( const FunctionSpaceType* functionSpaceImpl ) {
    return getElems( functionSpaceImpl ).global_index();
}

// Get UID field from StructuredColumns.
template <Meshless<StructuredColumns> = nullptr>
Field getUidField( const StructuredColumns* functionSpaceImpl ) {
    return functionSpaceImpl->global_index();
}

// Create UIDs from lonlat for PointCloud.
template <Meshless<PointCloud> = nullptr>
Field getUidField( const PointCloud* functionSpaceImpl ) {

    // Make UID field.
    Field uid = Field( " Uniquie ID ", array::make_datatype<uidx_t>(),
                       array::make_shape( functionSpaceImpl->size() ) );

    // Make views.
    auto uidView = array::make_view<uidx_t, 1>( uid );
    const auto lonLatView = array::make_view<double, 2>( functionSpaceImpl->lonlat() );

    // Set UIDs
    for ( idx_t i = 0; i < uid.shape(0); ++i ) {
        uidView( i ) = util::unique_lonlat( lonLatView (i, LON), lonLatView(i, LAT ) );
    }
    return uid;
}

// Get vector of (local index, UID) pairs.
template <typename FunctionSpaceType>
std::vector<IdxUid> getUidVec(const FunctionSpaceType* functionSpaceImpl) {

    // Get ghost and unique ID fields from functionspace.
    const auto ghost = getGhostField( functionSpaceImpl );
    const auto uid = getUidField( functionSpaceImpl );

    // Make views to fields.
    const auto ghostView = array::make_view<int, 1>( ghost );
    const auto uidView = array::make_view<uidx_t, 1>( uid );

    auto uidVec = std::vector<IdxUid>{};

    // Get UIDs for non ghost elems.
    for ( idx_t i = 0; i < functionSpaceImpl->size(); ++i ) {
        if ( ghostView( i ) ) {
            continue;
        }
       uidVec.push_back( IdxUid( i, uidView( i ) ) );
    }

    // Sort by UID value.
    std::sort( uidVec.begin(), uidVec.end(),
        []( const IdxUid& a, const IdxUid& b ){ return a.second < b.second; } );

    // Check UIDs are unique.
    if ( ATLAS_BUILD_TYPE_DEBUG ) {
        const size_t vecSize = uidVec.size();
        std::unique( uidVec.begin(), uidVec.end(),
            []( const IdxUid& a, const IdxUid& b ){ return a.second == b.second; });
        ATLAS_ASSERT( uidVec.size() == vecSize, "Unique ID set has duplicate members" );
    }

    return uidVec;
}

// Resolve functionspace implementation type and get unique ID vector.
std::vector<IdxUid> getUidVec(const FunctionSpace& functionSpace) {

    // Get implementation pointer.
    const auto* functionSpaceImpl = functionSpace.get();

    if ( functionSpaceImpl->cast<CellColumns>() ) {
        return getUidVec( functionSpaceImpl->cast<CellColumns>() );
    }
    if ( functionSpaceImpl->cast<EdgeColumns>() ) {
        return getUidVec( functionSpaceImpl->cast<EdgeColumns>() );
    }
    if ( functionSpaceImpl->cast<NodeColumns>() ) {
        return getUidVec( functionSpaceImpl->cast<NodeColumns>() );
    }
    if ( functionSpaceImpl->cast<PointCloud>() ) {
        return getUidVec( functionSpaceImpl->cast<PointCloud>() );
    }
    if ( functionSpaceImpl->cast<StructuredColumns>() ) {
        return getUidVec( functionSpaceImpl->cast<StructuredColumns>() );
    }

    const auto errorMessage =
        "Cannot redistribute functionspace type " + functionSpaceImpl->type();
    throw_NotImplemented( errorMessage, Here() );

}

// Copy UID index
std::vector<idx_t> getUidIdx( const std::vector<IdxUid>& uidVec ) {
    auto idxVec = std::vector<idx_t>{};
    std::transform( uidVec.begin(), uidVec.end(), std::back_inserter( idxVec ),
                   []( const IdxUid& uid ){ return uid.first; } );
    return idxVec;
}

// Copy UID value
std::vector<uidx_t> getUidVal( const std::vector<IdxUid>& uidVec ) {
    auto valVec = std::vector<uidx_t>{};
    std::transform( uidVec.begin(), uidVec.end(), std::back_inserter( valVec ),
                   []( const IdxUid& uid ){ return uid.second; } );
    return valVec;
}

// Communicate UID values, return receive buffer and displacements.
std::pair<std::vector<uidx_t>, std::vector<int>>
    communicateUid( const std::vector<uidx_t>& sendBuffer ) {

    auto counts = std::vector<int>( mpi::comm().size() );
    mpi::comm().allGather( static_cast<int>( sendBuffer.size() ),
                           counts.begin(), counts.end() );

    auto disps = std::vector<int>{ 0 };
    std::partial_sum( counts.begin(), counts.end(), std::back_inserter( disps ) );


    auto recvBuffer = std::vector<uidx_t>( static_cast<size_t>( disps.back() ) );

    mpi::comm().allGatherv( sendBuffer.begin(), sendBuffer.end(),
                            recvBuffer.begin(), counts.data(), disps.data() );

    return std::make_pair( recvBuffer, disps );
}

// Need to overload "<" operator as C++11 std::set_intersection does not support
// comparison lambda.
bool operator<( const IdxUid& lhs, const uidx_t& rhs ) {
    return lhs.second < rhs;
}
bool operator<( const uidx_t& lhs, const IdxUid& rhs ) {
    return lhs < rhs.second;
}

// Find the intersection between local and global UIDs, then return local
// indices of incections and PE dispacements in vector.
std::pair<std::vector<idx_t>, std::vector<int>>
    getUidIntersection( const std::vector<IdxUid>& localUids,
    const std::vector<uidx_t>& globalUids, const std::vector<int>& globalDisps ) {

    auto uidIntersection = std::vector<IdxUid>{};
    auto disps = std::vector<int>{ 0 };

    // Loop over all PE and find UID intersection.
    for ( size_t i = 0; i < mpi::comm().size(); ++i ) {

        // Get displaced iterators.
        auto globalUidsBegin = globalUids.begin() + globalDisps[i];
        auto globalUidsEnd = globalUids.begin() + globalDisps[i + 1];

        // Get intersection.
        std::set_intersection( localUids.begin(), localUids.end(),
            globalUidsBegin, globalUidsEnd, std::back_inserter( uidIntersection) );

        // Set intersection displacement.
        disps.push_back( static_cast<int>( uidIntersection.size() ) );

    }

    // Check that the set of all intersections matches UIDs on local PE.
    if ( ATLAS_BUILD_TYPE_DEBUG ) {
        auto tempUids = uidIntersection;
        std::sort( tempUids.begin(), tempUids.end(),
            []( const IdxUid& a, const IdxUid& b ){ return a.second < b.second; } );
        ATLAS_ASSERT( tempUids == localUids,
        "Set of all UID intersections does not match local UIDs." );
    }

    // Return local indices of intersection and displacements.
    return make_pair( getUidIdx( uidIntersection ), disps );
}


// Iterate over a field, in the order of an index list, and apply a functor to
// each element.

// Rank 1 overload.
template <typename Value, typename Functor>
void iterateField( const std::vector<idx_t>& idxList,
                   array::ArrayView<Value, 1>& fieldView, const Functor& f ) {

    for ( const idx_t i : idxList ) {
        f( fieldView( i ) );
    }
}

// Rank 2 overload.
template <typename Value, typename Functor>
void iterateField( const std::vector<idx_t>& idxList,
                   array::ArrayView<Value, 2>& fieldView, const Functor f ) {

    for ( const idx_t i : idxList ) {
        for ( idx_t j = 0; j < fieldView.shape(1); ++j ) {
            f( fieldView( i, j ) );
        }
    }
}

// Rank 3 overload.
template <typename Value, typename Functor>
void iterateField( const std::vector<idx_t>& idxList,
                   array::ArrayView<Value, 3>& fieldView, const Functor f ) {

    for ( const idx_t i : idxList ) {
        for ( idx_t j = 0; j < fieldView.shape(1); ++j ) {
            for ( idx_t k = 0; k < fieldView.shape(2); ++k )
                f( fieldView( i, j, k ) );
        }
    }
}

}

void RedistributeGeneric::setup(const FunctionSpace &sourceFunctionSpace,
                                const FunctionSpace &targetFunctionSpace) {

    // Assign function spaces.
    source() = sourceFunctionSpace;
    target() = targetFunctionSpace;

    // get a unique ID (UID) for each owned member of functionspace.
    const auto sourceUidVec = getUidVec( source() );
    const auto targetUidVec = getUidVec( target() );

    // Communicate UID vectors to all PEs.
    auto sourceGlobalUids = std::vector<uidx_t>{};
    auto sourceGlobalDisps = std::vector<int>{};
    std::tie( sourceGlobalUids, sourceGlobalDisps ) =
        communicateUid( getUidVal( sourceUidVec ) );
    auto targetGlobalUids = std::vector<uidx_t>{};
    auto targetGlobalDisps = std::vector<int>{};
    std::tie( targetGlobalUids, targetGlobalDisps ) =
        communicateUid( getUidVal( targetUidVec ) );

    // Get intersection of local UIDs and Global UIDs.
    std::tie( sourceLocalIdx_, sourceDisps_ ) =
        getUidIntersection( sourceUidVec, targetGlobalUids, targetGlobalDisps );
    std::tie( targetLocalIdx_, targetDisps_ ) =
        getUidIntersection( targetUidVec, sourceGlobalUids, sourceGlobalDisps );

}

void RedistributeGeneric::execute( const Field &sourceField, Field &targetField ) const {

    //Check functionspaces match.
    ATLAS_ASSERT( sourceField.functionspace().type() == source().type() );
    ATLAS_ASSERT( targetField.functionspace().type() == target().type() );

    // Check Field datatypes match.
    ATLAS_ASSERT( sourceField.datatype() == targetField.datatype() );

    // Check Field ranks match.
    ATLAS_ASSERT( sourceField.rank() == targetField.rank() );

    // Check number of levels and variables match.
    for ( idx_t i = 1; i < sourceField.rank(); ++i ) {
        ATLAS_ASSERT( sourceField.shape(i) == targetField.shape()[i] );
    }

    // Perform redistribution.
    doExecute( sourceField, targetField );
}

void RedistributeGeneric::execute( const FieldSet &sourceFieldSet, FieldSet &targetFieldSet ) const {

    // Check field set sizes match.
    ATLAS_ASSERT( sourceFieldSet.size() == targetFieldSet.size() );

    // Redistribute fields.
    for ( idx_t i = 0; i < sourceFieldSet.size(); ++i ) {
        execute( sourceFieldSet[i], targetFieldSet[i] );
    }
}

// Determine datatype.
void RedistributeGeneric::doExecute( const Field& sourceField, Field& targetField ) const {

    switch ( sourceField.datatype().kind() ) {
        case array::DataType::KIND_REAL64 : {
            doExecute<double>( sourceField, targetField );
            break;
        }
        case array::DataType::KIND_REAL32 : {
            doExecute<float>( sourceField, targetField );
            break;
        }
        case array::DataType::KIND_INT64 : {
            doExecute<long>( sourceField, targetField );
            break;
        }
        case array::DataType::KIND_INT32 : {
            doExecute<int>( sourceField, targetField );
            break;
        }
        default : {
            throw_NotImplemented( "No implementation for data type " +
                                  sourceField.datatype().str(), Here() );
        }
    }
}

// Determine rank.
template <typename Value>
void RedistributeGeneric::doExecute( const Field& sourceField, Field& targetField ) const {

    switch ( sourceField.rank() ) {
        case 1 : {
            doExecute<Value, 1>( sourceField, targetField );
            break;
        }
        case 2 : {
            doExecute<Value, 2>( sourceField, targetField );
            break;
        }
        case 3 : {
            doExecute<Value, 3>( sourceField, targetField );
            break;
        }
        default : {
            throw_NotImplemented( "No implementation for rank " +
                                  std::to_string ( sourceField.rank() ), Here() );
        }
    }
}

// Perform redistribution.
template <typename Value, int Rank>
void RedistributeGeneric::doExecute( const Field& sourceField, Field& targetField ) const {

    // Get array views.
    auto sourceView = array::make_view<Value, Rank>( sourceField );
    auto targetView = array::make_view<Value, Rank>( targetField );

    // Get number of elems per column.
    int elemsPerCol = 1;
    for ( size_t i = 1; i < Rank; ++i ) {
        elemsPerCol *= sourceView.shape()[i];
    }

    // Set send displacement and counts vectors.
    auto sendDisps = std::vector<int>{};
    auto sendCounts = std::vector<int>{};
    std::transform( sourceDisps_.begin(), sourceDisps_.end(),
                    std::back_inserter( sendDisps ),
                    [&]( const int& disp ){ return disp * elemsPerCol; } );
    std::adjacent_difference( sendDisps.begin() + 1, sendDisps.end(),
                              std::back_inserter( sendCounts ) );

    // Set recv displacement and counts vectors.
    auto recvDisps = std::vector<int>{};
    auto recvCounts = std::vector<int>{};
    std::transform( targetDisps_.begin(), targetDisps_.end(),
                    std::back_inserter( recvDisps ),
                    [&]( const int& disp ){ return disp * elemsPerCol; } );
    std::adjacent_difference( recvDisps.begin() + 1, recvDisps.end(),
                              std::back_inserter( recvCounts) );

    // Allocate send and recv buffers.
    auto sendBuffer = std::vector<Value>( static_cast<size_t>( sendDisps.back() ) );
    auto recvBuffer = std::vector<Value>( static_cast<size_t>( recvDisps.back() ) );
    auto sendBufferIt = sendBuffer.begin();
    auto recvBufferIt = recvBuffer.cbegin();

    // Copy sourceField to sendBuffer.
    iterateField( sourceLocalIdx_, sourceView,
                  [&]( const Value& elem ){ *sendBufferIt++ = elem; } );

    // Perform MPI communication.
    mpi::comm().allToAllv( sendBuffer.data(), sendCounts.data(), sendDisps.data(),
                           recvBuffer.data(), recvCounts.data(), recvDisps.data() );

    // Copy recvBuffer to targetField.
    iterateField( targetLocalIdx_, targetView,
                  [&]( Value& elem ){ elem = *recvBufferIt++; } );

}

namespace {
static RedistributionImplBuilder<RedistributeGeneric>
  register_builder( RedistributeGeneric::static_type() );
}

} // namespace detail
} // namespace redistribution
} // namespace atlas
