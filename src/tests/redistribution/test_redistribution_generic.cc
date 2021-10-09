/*
 * (C) Crown Copyright 2021 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <cmath>

#include "atlas/field/FieldSet.h"
#include "atlas/functionspace.h"
#include "atlas/grid.h"
#include "atlas/grid/Partitioner.h"
#include "atlas/meshgenerator.h"
#include "atlas/output/Gmsh.h"
#include "atlas/redistribution/Redistribution.h"
#include "atlas/util/Config.h"
#include "atlas/util/Constants.h"
#include "atlas/util/Unique.h"

#include "tests/AtlasTestEnvironment.h"

namespace atlas {
namespace test {

// Set floating point tolerance.
template<typename Value>
Value tolerance() { return  std::numeric_limits<double>::epsilon() * 8; }

// Set field config for different ranks.
template<int Rank>
util::Config fieldConfig();
// Rank 1 config.
template<>
util::Config fieldConfig<1>() {
    return util::Config();
}
// Rank 2 config.
template<>
util::Config fieldConfig<2>() {
    auto config = util::Config();
    config.set( "levels", 10 );
    return config;
}
// Rank 3 config.
template<>
util::Config fieldConfig<3>() {
    auto config = util::Config();
    config.set( "levels", 10 );
    config.set( "variables", 2 );
    return config;
}

// Define test pattern for grid.
template <typename Value>
Value testPattern( double lambda, double phi, idx_t level ) {
    return static_cast<Value>( std::cos( lambda * ( 1 + level ) * M_PI / 180. ) *
                           std::sin( phi * ( 1 + level ) * M_PI / 180. ) );
}

// Try and get a mesh from functionspace.
const Mesh* getMesh( const FunctionSpace& functionSpace) {

    // Try and create a pointer to one of the mesh based functionspaces.
    const auto* cellColumnsPtr = functionSpace.get()->cast<functionspace::detail::CellColumns>();
    const auto* edgeColumnsPtr = functionSpace.get()->cast<functionspace::detail::EdgeColumns>();
    const auto* nodeColumnsPtr = functionSpace.get()->cast<functionspace::detail::NodeColumns>();

    if ( cellColumnsPtr ) {
        return &( cellColumnsPtr->mesh() );
    }
    else if ( edgeColumnsPtr ) {
        return &( edgeColumnsPtr->mesh() );
    }
    else if ( nodeColumnsPtr ) {
        return &( nodeColumnsPtr->mesh() );
    }
    else {
        return nullptr;
    }
}

// Try and get cells or edges from functionspace.
const mesh::HybridElements* getHybridElems( const FunctionSpace& functionSpace ) {

    // Try and create a pointer to CellColumns or EdgeColumns.
    const auto* cellColumnsPtr = functionSpace.get()->cast<functionspace::detail::CellColumns>();
    const auto* edgeColumnsPtr = functionSpace.get()->cast<functionspace::detail::EdgeColumns>();

    if ( cellColumnsPtr ) {
        return &( cellColumnsPtr->cells() );
    }
    else if ( edgeColumnsPtr ) {
        return &( edgeColumnsPtr->edges() );
    }
    else {
        return nullptr;
    }

}

// Class to test functionspace redistribution.
template <typename Value, int Rank>
struct TestRedistribution {
public:

    TestRedistribution( const FunctionSpace& sourceFunctionSpace,
        const FunctionSpace& targetFunctionSpace ) :
        sourceFunctionSpace_( sourceFunctionSpace ), targetFunctionSpace_(targetFunctionSpace),
        redist_( sourceFunctionSpace, targetFunctionSpace ),
        sourceField_( sourceFunctionSpace_.createField<Value>( fieldConfig<Rank>() ) ),
        targetField_( targetFunctionSpace_.createField<Value>( fieldConfig<Rank>() ) ),
        sourceView_( array::make_view<Value, Rank>(sourceField_) ),
        targetView_( array::make_view<Value, Rank>(targetField_) ) {}

    void outputFields() {

        // Try and create a pointer to one of the mesh based functionspaces.
        const auto* sourceCellColumnsPtr = this->sourceFunctionSpace_.get()->template
                                     cast<functionspace::detail::CellColumns>();
        const auto* sourceEdgeColumnsPtr = this->sourceFunctionSpace_.get()->template
                                     cast<functionspace::detail::EdgeColumns>();
        const auto* sourceNodeColumnsPtr = this->sourceFunctionSpace_.get()->template
                                     cast<functionspace::detail::NodeColumns>();

        const auto* targetCellColumnsPtr = this->targetFunctionSpace_.get()->template
                                     cast<functionspace::detail::CellColumns>();
        const auto* targetEdgeColumnsPtr = this->targetFunctionSpace_.get()->template
                                     cast<functionspace::detail::EdgeColumns>();
        const auto* targetNodeColumnsPtr = this->targetFunctionSpace_.get()->template
                                     cast<functionspace::detail::NodeColumns>();



    }


    FunctionSpace sourceFunctionSpace_;
    FunctionSpace targetFunctionSpace_;

    Redistribution redist_;

    Field sourceField_;
    Field targetField_;

    array::ArrayView<Value, Rank> sourceView_;
    array::ArrayView<Value, Rank> targetView_;

};

// Test rank 1 fields with lonlat method.
template <typename Value>
struct TestRedistributionPoints1 : public TestRedistribution<Value, 1> {
    using TestRedistribution<Value, 1>::TestRedistribution;
    void execute() {

        auto sourceLonlatView = array::make_view<double, 2>( this->sourceFunctionSpace_.lonlat() );
        auto targetLonlatView = array::make_view<double, 2>( this->targetFunctionSpace_.lonlat() );

        // Set source field.
        for ( idx_t i = 0; i < this->sourceView_.shape( 0 ); ++i ) {
            this->sourceView_( i ) =
                testPattern<Value>( sourceLonlatView( i, 0 ),
                                    sourceLonlatView( i, 1 ), 0 );
        }

        // Perform redistribution.
        this->redist_.execute( this->sourceField_, this->targetField_ );

        // Perform halo excahnge;
        this->targetFunctionSpace_.haloExchange( this->targetField_ );

        // Check target field.
        for ( idx_t i = 0; i < this->targetView_.shape( 0 ); ++i ) {
            EXPECT_APPROX_EQ( this->targetView_( i ),
                testPattern<Value>( targetLonlatView( i, 0 ),
                                    targetLonlatView( i, 1 ), 0 ), tolerance<Value>() );
        }

    }
};

// Test rank 2 fields with lonlat method.
template <typename Value>
struct TestRedistributionPoints2 : public TestRedistribution<Value, 2> {
    using TestRedistribution<Value, 2>::TestRedistribution;
    void execute() {

        auto sourceLonlatView = array::make_view<double, 2>( this->sourceFunctionSpace_.lonlat() );
        auto targetLonlatView = array::make_view<double, 2>( this->targetFunctionSpace_.lonlat() );

        // Set source field.
        for ( idx_t i = 0; i < this->sourceView_.shape( 0 ); ++i ) {
            for ( idx_t j = 0; j < this->sourceView_.shape( 1 ); ++j ) {
                this->sourceView_( i, j ) =
                    testPattern<Value>( sourceLonlatView( i, 0 ),
                                        sourceLonlatView( i, 1 ), j );
            }
        }

        // Perform redistribution.
        this->redist_.execute( this->sourceField_, this->targetField_ );

        // Perform halo excahnge;
        this->targetFunctionSpace_.haloExchange( this->targetField_);

        // Check target field.
        for ( idx_t i = 0; i < this->targetView_.shape( 0 ); ++i ) {
            for ( idx_t j = 0; j < this->targetView_.shape( 1 ); ++j ) {
                EXPECT_APPROX_EQ( this->targetView_( i, j ),
                    testPattern<Value>( targetLonlatView( i, 0 ),
                                        targetLonlatView( i, 1 ), j ), tolerance<Value>() );
            }
        }

    }
};

// Test rank 3 fields with lonlat method.
template <typename Value>
struct TestRedistributionPoints3 : public TestRedistribution<Value, 3> {
    using TestRedistribution<Value, 3>::TestRedistribution;
    void execute() {

        auto sourceLonlatView = array::make_view<double, 2>( this->sourceFunctionSpace_.lonlat() );
        auto targetLonlatView = array::make_view<double, 2>( this->targetFunctionSpace_.lonlat() );

        // Set source field.
        for ( idx_t i = 0; i < this->sourceView_.shape( 0 ); ++i ) {
            for ( idx_t j = 0; j < this->sourceView_.shape( 1 ); ++j ) {
                this->sourceView_( i, j, 0 ) =
                    testPattern<Value>( sourceLonlatView( i, 0 ),
                                        sourceLonlatView( i, 1 ), j );
                this->sourceView_( i, j, 1 ) =
                    -testPattern<Value>( sourceLonlatView( i, 0 ),
                                         sourceLonlatView( i, 1 ), j );
            }
        }

        // Perform redistribution.
        this->redist_.execute( this->sourceField_, this->targetField_ );

        // Perform halo excahnge;
        this->targetFunctionSpace_.haloExchange( this->targetField_);

        // Check target field.
        for ( idx_t i = 0; i < this->targetView_.shape( 0 ); ++i ) {
            for ( idx_t j = 0; j < this->targetView_.shape( 1 ); ++j ) {
                EXPECT_APPROX_EQ( this->targetView_( i, j, 0 ),
                    testPattern<Value>( targetLonlatView( i, 0 ),
                                        targetLonlatView( i, 1 ), j ), tolerance<Value>() );
                EXPECT_APPROX_EQ( this->targetView_( i, j, 1 ),
                    -testPattern<Value>( targetLonlatView( i, 0 ),
                                         targetLonlatView( i, 1 ), j ), tolerance<Value>() );
            }
        }

    }
};

// Test CellColumns or EdgeColumns fields.
struct TestRedistributionElems : public TestRedistribution<uidx_t, 1> {
    using TestRedistribution<uidx_t, 1>::TestRedistribution;
    void execute() {

        // Non lonlat method defined. Test redistrubtion of fields of unique
        // IDs instead of test_pattern.

        // Get source mesh and elements from functionspace.
        const Mesh* sourceMeshPtr = getMesh( this->sourceFunctionSpace_ ) ;
        const mesh::HybridElements* sourceElemsPtr = getHybridElems( this->sourceFunctionSpace_ );

        // Get target mesh and elements from functionspace.
        const Mesh* targetMeshPtr = getMesh( this->targetFunctionSpace_ ) ;
        const mesh::HybridElements* targetElemsPtr = getHybridElems( this->targetFunctionSpace_ );

        // Create unique ID generators and connectivity graphs for source functionspace.
        const auto sourceUidGenerator = util::UniqueLonLat( *sourceMeshPtr );
        const auto& sourceConnectivity = sourceElemsPtr->node_connectivity();

        // Create unique ID generators and connectivity graphs for target functionspace.
        const auto targetUidGenerator = util::UniqueLonLat( *targetMeshPtr );
        const auto& targetConnectivity = targetElemsPtr->node_connectivity();

        // Set source field.
        for ( idx_t i = 0; i < this->sourceView_.shape( 0 ); ++i ) {
            this->sourceView_( i ) = sourceUidGenerator( sourceConnectivity.row( i ) );
        }

        // Perform redistribution.
        this->redist_.execute( this->sourceField_, this->targetField_ );

        // Perform halo excahnge;
        this->targetFunctionSpace_.haloExchange(this->targetField_);

        // Check target field.
        for ( idx_t i = 0; i < this->targetView_.shape( 0 ); ++i ) {
            EXPECT( this->targetView_( i ) == targetUidGenerator( targetConnectivity.row( i ) ));
        }

    }
};


CASE( "Generic" ) {

    auto grid = atlas::Grid( "L48x37" );

    // Set mesh config.
    const auto sourceMeshConfig = util::Config( "partitioner", "checkerboard" );
    const auto targetMeshConfig = util::Config( "partitioner", "equal_regions" );

    auto sourceMesh = MeshGenerator( "structured", sourceMeshConfig ).generate(grid);
    auto targetMesh = MeshGenerator( "structured", targetMeshConfig ).generate(grid);

    const auto sourceFunctionSpace = functionspace::NodeColumns(sourceMesh);
    const auto targetFunctionSpace = functionspace::NodeColumns(targetMesh);

    const auto sourceCells = functionspace::CellColumns(sourceMesh);
    const auto targetCells = functionspace::CellColumns(targetMesh);

    TestRedistributionPoints1<double>(sourceFunctionSpace, targetFunctionSpace).execute();
    TestRedistributionPoints2<double>(sourceFunctionSpace, targetFunctionSpace).execute();
    TestRedistributionPoints3<double>(sourceFunctionSpace, targetFunctionSpace).execute();

    TestRedistributionPoints1<float>(sourceFunctionSpace, targetFunctionSpace).execute();
    TestRedistributionPoints2<float>(sourceFunctionSpace, targetFunctionSpace).execute();
    TestRedistributionPoints3<float>(sourceFunctionSpace, targetFunctionSpace).execute();

    TestRedistributionElems(sourceCells, targetCells).execute();

}

}  // namespace test
}  // namespace atlas


int main( int argc, char** argv ) {
    return atlas::test::run( argc, argv );
}
