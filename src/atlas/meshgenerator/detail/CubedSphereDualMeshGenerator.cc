/*
 * (C) Crown Copyright 2021 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <string>
#include <vector>

#include "eckit/utils/Hash.h"


#include "atlas/array/helpers/ArrayCopier.h"
#include "atlas/functionspace/CubedSphereColumns.h"
#include "atlas/grid/CubedSphereGrid.h"
#include "atlas/grid/Distribution.h"
#include "atlas/grid/Partitioner.h"
#include "atlas/library/config.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/mesh/Nodes.h"
#include "atlas/meshgenerator/MeshGenerator.h"
#include "atlas/meshgenerator/detail/CubedSphereDualMeshGenerator.h"
#include "atlas/meshgenerator/detail/MeshGeneratorFactory.h"
#include "atlas/meshgenerator/detail/cubedsphere/CubedSphereUtility.h"
#include "atlas/parallel/mpi/mpi.h"

#define DEBUG_OUTPUT 0
#define DEBUG_OUTPUT_DETAIL 0

namespace atlas {
namespace meshgenerator {

// -----------------------------------------------------------------------------

CubedSphereDualMeshGenerator::CubedSphereDualMeshGenerator( const eckit::Parametrisation& p ) {
    configure_defaults();

    // Get number of partitions.
    size_t nb_parts;
    if ( p.get( "nb_parts", nb_parts ) ) {
        options.set( "nb_parts", nb_parts );
    }

    // Get this partition.
    int part;
    if ( p.get( "part", part ) ) {
        options.set( "part", part );
    }

    // Get halo size.
    int halo;
    if ( p.get( "halo", halo ) ) {
        options.set( "halo", halo );
    }

    // Get partitioner.
    std::string partitioner;
    if ( p.get( "partitioner", partitioner ) ) {
        options.set( "partitioner", partitioner );
    }
}

// -----------------------------------------------------------------------------


void CubedSphereDualMeshGenerator::configure_defaults() {
    // This option sets number of partitions.
    options.set( "nb_parts", mpi::size() );

    // This option sets the part that will be generated.
    options.set( "part", mpi::rank() );

    // This options sets the number of halo elements around each partition.
    options.set( "halo", 1 );

    // This options sets the default partitioner.
    options.set<std::string>( "partitioner", "cubedsphere" );
}

// -----------------------------------------------------------------------------

void CubedSphereDualMeshGenerator::generate( const Grid& grid, Mesh& mesh ) const {
    // Get partitioner type and number of partitions from config.
    const idx_t nParts         = static_cast<idx_t>( options.get<size_t>( "nb_parts" ) );
    const std::string partType = options.get<std::string>( "partitioner" );

    auto partConfig = util::Config{};
    partConfig.set( "type", partType );
    partConfig.set( "partitions", nParts );

    // Use lonlat instead of xy for non cubedsphere partitioner.
    if ( partType != "cubedsphere" ) {
        partConfig.set( "coordinates", "lonlat" );
    }

    // Set distribution.
    const auto partitioner  = grid::Partitioner( partConfig );
    const auto distribution = grid::Distribution( grid, partitioner );

    generate( grid, distribution, mesh );
}

// -----------------------------------------------------------------------------

void CubedSphereDualMeshGenerator::generate( const Grid& grid, const grid::Distribution& distribution, Mesh& mesh ) const {
    // Check for correct grid and need for mesh
    ATLAS_ASSERT( !mesh.generated() );

    // Cast grid to cubed sphere grid.
    const auto csGrid = CubedSphereGrid( grid );

    // Check for successful cast.
    if ( !csGrid ) {
        throw_Exception(
            "CubedSphereDualMeshGenerator can only work with a cubedsphere grid.", Here() );
    }

    // Check for correct grid stagger.
    if ( csGrid.stagger() != "C" ) {
        throw_Exception(
            "CubedSphereDualMeshGenerator will only work with a cell-centroid grid.", Here() );
    }

    // Enforce compatible halo size.
    if ( options.get<idx_t>( "halo" ) != 1 ) {
        throw_Exception( "Halo size CubedSphereDualMeshGenerator is currently "
                         "limited to 1.", Here() );
    }

    // Clone some grid properties.
    setGrid( mesh, csGrid, distribution );

    generate_mesh( csGrid, distribution, mesh );
}

// -----------------------------------------------------------------------------

namespace  {

// Helper function to copy fields.
template<typename Value, int Rank>
void copyField(const Field& sourceField, Field& targetField) {

    // Make views to field.
    const auto sourceView = array::make_view<Value, Rank>( sourceField );
    auto targetView = array::make_view<Value, Rank>( targetField );

    // Assign source field values to target field.
    array::helpers::array_copier<Value, Rank>::apply( sourceView, targetView );
}

}

void CubedSphereDualMeshGenerator::generate_mesh( const CubedSphereGrid& csGrid, const grid::Distribution& distribution,
                                              Mesh& mesh ) const {
    ATLAS_TRACE( "CubedSphereDualMeshGenerator::generate" );

    using namespace detail::cubedsphere;

    const idx_t N         = csGrid.N();
    const idx_t nHalo     = options.get<int>( "halo" );

    //--------------------------------------------------------------------------
    // Create a cubed-sphere mesh.
    //--------------------------------------------------------------------------

    // Generate cubed sphere mesh.
    const auto csMesh =
        MeshGenerator("cubedsphere", options).generate( csGrid, distribution );

    // Generate fucntionspaces cubed sphere mesh.
    const auto csCellsFunctionSpace = functionspace::CubedSphereCellColumns( csMesh );
    const auto csNodesFunctionSpace = functionspace::CubedSphereNodeColumns( csMesh );
    const auto& csCells = csCellsFunctionSpace.cells();
    const auto& csNodes = csNodesFunctionSpace.nodes();

    //--------------------------------------------------------------------------
    // Set dual mesh nodes (easy part).
    //--------------------------------------------------------------------------

    auto& nodes = mesh.nodes();
    nodes.resize( csCellsFunctionSpace.size() );

    nodes.add( Field("tij", array::make_datatype<idx_t>(),
                            array::make_shape( nodes.size(), 3) ) );

    // Copy mesh fields to dual mesh.
    copyField<gidx_t, 1>( csCells.global_index(), nodes.global_index() );
    copyField<idx_t, 1>( csCells.remote_index(), nodes.remote_index() );
    copyField<int, 1>( csCells.partition(), nodes.partition() );
    copyField<int, 1>( csCells.halo(), nodes.halo() );
    copyField<int, 1>( csCells.flags(), nodes.flags() );
    copyField<double, 2>( csCells.field( "xy" ), nodes.xy() );
    copyField<double, 2>( csCells.field( "lonlat" ), nodes.lonlat() );
    copyField<idx_t, 2>( csCells.field( "tij"), nodes.field("tij") );

    // Need to set a ghost field and decrement halo by one.
    auto nodesGhost = array::make_view<int, 1>( nodes.ghost() );
    auto nodesHalo = array::make_view<int, 1>( nodes.halo() );

    for ( idx_t idx = 0; idx < nodes.size(); ++idx ) {
        nodesGhost( idx ) = nodesHalo( idx ) > 0;
    }

    //--------------------------------------------------------------------------
    // Set dual mesh cells (not so easy part).
    //--------------------------------------------------------------------------

    // Make views to cubed sphere nodes.
    const auto csNodesGlobalIdx = array::make_view<gidx_t, 1>( csNodes.global_index() );
    const auto csNodesRemoteIdx = array::make_indexview<idx_t, 1>( csNodes.remote_index() );
    const auto csNodesXy        = array::make_view<double, 2>( csNodes.xy() );
    const auto csNodesLonLat    = array::make_view<double, 2>( csNodes.lonlat() );
    const auto csNodesPart      = array::make_view<int, 1>( csNodes.partition() );
    const auto csNodesHalo      = array::make_view<int, 1>( csNodes.halo() );
    const auto csNodesFlags     = array::make_view<int, 1>( csNodes.flags() );
    const auto csNodesTij       = array::make_view<idx_t, 2>( csNodes.field( "tij" ) );

    // Check if cell is in tile corner.
    const auto cornerCell = [&]( idx_t idxCell ) -> bool {

        // Nodes of csMesh are cells of dual mesh.

        const idx_t i = csNodesTij( idxCell, Coordinates::I );
        const idx_t j = csNodesTij( idxCell, Coordinates::J );
        return ( i == 0 && j == 0 ) || ( i == N && j == 0 ) ||
               ( i == N && j == N ) || ( i == 0 && j == N );
    };

    // Check if cell is out of bounds

    // Figure out the number of quads and triangles we need to add.
    idx_t nQuadCell = 0;
    idx_t nTriCell = 0;
    for ( idx_t idxCsNode = 0; idxCsNode < csNodes.size(); ++idxCsNode ) {

        // Exclude outer halo ring of cubed sphere mesh nodes.
        if ( csNodesHalo( idxCsNode ) < nHalo ) {

            // Tile corner elements are triangles. Rest are quads.
            if ( cornerCell( idxCsNode ) ) {
                ++nTriCell;
            }
            else {
                ++nQuadCell;
            }
        }
    }

    // Add cells to mesh.
    auto& cells = mesh.cells();
    cells.add( new mesh::temporary::Quadrilateral(), nQuadCell );
    cells.add( new mesh::temporary::Triangle(), nTriCell );
    auto& nodeConnectivity = cells.node_connectivity();

    // Add extra fields.
    auto tijField = cells.add( Field( "tij", array::make_datatype<idx_t>(),
                                      array::make_shape( cells.size(), 3 ) ) );
    tijField.set_variables( 3 );

    Field xyField = cells.add( Field( "xy", array::make_datatype<double>(),
                                      array::make_shape( cells.size(), 2 ) ) );
    xyField.set_variables( 2 );

    Field lonLatField = cells.add( Field( "lonlat", array::make_datatype<double>(),
                                          array::make_shape( cells.size(), 2 ) ) );
    lonLatField.set_variables( 2 );

    // Make views to dual mesh fields.
    auto cellsGlobalIdx = array::make_view<gidx_t, 1>( cells.global_index() );
    auto cellsRemoteIdx = array::make_indexview<idx_t, 1>( cells.remote_index() );
    auto cellsPart      = array::make_view<int, 1>( cells.partition() );
    auto cellsHalo      = array::make_view<int, 1>( cells.halo() );
    auto cellsFlags     = array::make_view<int, 1>( cells.flags() );
    auto cellsTij       = array::make_view<idx_t, 2>( tijField );
    auto cellsXy        = array::make_view<double, 2>( xyField );
    auto cellsLonLat    = array::make_view<double, 2>( lonLatField );

    // Set separate counters for quads and triangles.
    idx_t idxQuadCell = cells.elements( 0 ).begin();
    idx_t idxTriCell = cells.elements( 1 ).begin();

    // Loop over mesh nodes and set dual mesh cells.
    for ( idx_t idxCsNode = 0; idxCsNode < nQuadCell + nTriCell; ++idxCsNode ) {

        const idx_t tCell = csNodesTij( idxCsNode, Coordinates::T );
        const idx_t iCell = csNodesTij( idxCsNode, Coordinates::I );
        const idx_t jCell = csNodesTij( idxCsNode, Coordinates::J );

        idx_t idxCell;

        // Set corner triangle.
        if ( cornerCell( idxCsNode ) ) {
            idxCell = idxTriCell++;

            auto triCellNodes = std::array<idx_t, 3>{};

            // Bottom-left corner.
            if ( iCell == 0 && jCell == 0 ) {
                triCellNodes = {
                    csCellsFunctionSpace.index( tCell,  0, -1 ),
                    csCellsFunctionSpace.index( tCell,  0,  0 ),
                    csCellsFunctionSpace.index( tCell, -1,  0 )};
            }

            // Bottom-right corner.
            else if ( iCell == N && jCell == 0 ) {
                triCellNodes = {
                    csCellsFunctionSpace.index( tCell, N - 1, -1 ),
                    csCellsFunctionSpace.index( tCell, N    ,  0 ),
                    csCellsFunctionSpace.index( tCell, N - 1,  0 )};
            }

            // Top-right corner.
            else if ( iCell == N && jCell == N ) {
                triCellNodes = {
                    csCellsFunctionSpace.index( tCell, N - 1, N - 1 ),
                    csCellsFunctionSpace.index( tCell, N    , N - 1 ),
                    csCellsFunctionSpace.index( tCell, N - 1, N     )};
            }

            // Top-left corner.
            else {
                triCellNodes = {
                    csCellsFunctionSpace.index( tCell, -1, N - 1 ),
                    csCellsFunctionSpace.index( tCell,  0, N - 1 ),
                    csCellsFunctionSpace.index( tCell,  0, N     )};
            }

            nodeConnectivity.set( idxCell, triCellNodes.data());

        }

        // Set quad.
        else {
            idxCell = idxQuadCell++;

            auto quadCellNodes = std::array<idx_t, 4>{
                csCellsFunctionSpace.index( tCell, iCell - 1, jCell - 1 ),
                csCellsFunctionSpace.index( tCell, iCell    , jCell - 1 ),
                csCellsFunctionSpace.index( tCell, iCell    , jCell     ),
                csCellsFunctionSpace.index( tCell, iCell - 1, jCell     ),
            };

            nodeConnectivity.set( idxCell, quadCellNodes.data() );

        }

        // Copy field data.
        cellsGlobalIdx(idxCell) = csNodesGlobalIdx(idxCsNode);
        cellsRemoteIdx(idxCell) = csNodesRemoteIdx(idxCsNode);
        cellsPart(idxCell) = csNodesPart(idxCsNode);
        cellsHalo(idxCell) = csNodesHalo(idxCsNode);
        cellsFlags(idxCell) = csNodesFlags(idxCsNode);
        cellsTij(idxCell, Coordinates::T) = csNodesTij(idxCsNode, Coordinates::T);
        cellsTij(idxCell, Coordinates::I) = csNodesTij(idxCsNode, Coordinates::I);
        cellsTij(idxCell, Coordinates::J) = csNodesTij(idxCsNode, Coordinates::J);
        cellsXy(idxCell, XX) = csNodesXy(idxCsNode, XX);
        cellsXy(idxCell, YY) = csNodesXy(idxCsNode, YY);
        cellsLonLat(idxCell, LON) = csNodesLonLat(idxCsNode, LON);
        cellsLonLat(idxCell, LAT) = csNodesLonLat(idxCsNode, LAT);

    }

    // Set metadata.
    mesh.metadata().set( "halo", nHalo );
    mesh.metadata().set( "halo_locked", true );
    mesh.nodes().metadata().set( "parallel", true );
    mesh.cells().metadata().set( "parallel", true );

    return;
}


// -----------------------------------------------------------------------------

void CubedSphereDualMeshGenerator::hash( eckit::Hash& h ) const {
    h.add( "CubedSphereDualMeshGenerator" );
    options.hash( h );
}

// -----------------------------------------------------------------------------

namespace {
static MeshGeneratorBuilder<CubedSphereDualMeshGenerator> CubedSphereDualMeshGenerator(
    CubedSphereDualMeshGenerator::static_type() );
}

// -----------------------------------------------------------------------------

}  // namespace meshgenerator
}  // namespace atlas
