/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <vector>

#include "eckit/utils/Hash.h"

#include "atlas/array/Array.h"
#include "atlas/array/ArrayView.h"
#include "atlas/array/IndexView.h"
#include "atlas/field/Field.h"
#include "atlas/grid/Distribution.h"
#include "atlas/grid/Partitioner.h"
#include "atlas/grid/StructuredGrid.h"
#include "atlas/library/config.h"
#include "atlas/mesh/ElementType.h"
#include "atlas/mesh/Elements.h"
#include "atlas/mesh/HybridElements.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/mesh/Nodes.h"
#include "atlas/meshgenerator/detail/MeshGeneratorFactory.h"
#include "atlas/meshgenerator/detail/RegularMeshGenerator.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/runtime/Exception.h"
#include "atlas/runtime/Log.h"
#include "atlas/util/CoordinateEnums.h"

#define DEBUG_OUTPUT 0
#define DEBUG_OUTPUT_DETAIL 0

using Topology = atlas::mesh::Nodes::Topology;

namespace atlas {
namespace meshgenerator {

RegularMeshGenerator::RegularMeshGenerator( const eckit::Parametrisation& p ) {
    configure_defaults();

    // options copied from Structured MeshGenerator
    size_t nb_parts;
    if ( p.get( "nb_parts", nb_parts ) ) {
        options.set( "nb_parts", nb_parts );
    }

    size_t part;
    if ( p.get( "part", part ) ) {
        options.set( "part", part );
    }

    std::string partitioner;
    if ( p.get( "partitioner", partitioner ) ) {
        if ( not grid::Partitioner::exists( partitioner ) ) {
            Log::warning() << "Atlas does not have support for partitioner " << partitioner << ". "
                           << "Defaulting to use partitioner EqualRegions" << std::endl;
            partitioner = "equal_regions";
        }
        options.set( "partitioner", partitioner );
    }

    // options specifically for this MeshGenerator
    bool periodic_x;
    if ( p.get( "periodic_x", periodic_x ) ) {
        options.set( "periodic_x", periodic_x );
    }

    bool periodic_y;
    if ( p.get( "periodic_y", periodic_y ) ) {
        options.set( "periodic_y", periodic_y );
    }

    bool biperiodic;
    if ( p.get( "biperiodic", biperiodic ) ) {
        options.set( "periodic_x", biperiodic );
        options.set( "periodic_y", biperiodic );
    }
}

void RegularMeshGenerator::configure_defaults() {
    // This option sets number of parts the mesh will be split in
    options.set( "nb_parts", mpi::size() );

    // This option sets the part that will be generated
    options.set( "part", mpi::rank() );

    // This options sets the default partitioner
    std::string partitioner;
    if ( grid::Partitioner::exists( "trans" ) && mpi::size() > 1 ) {
        partitioner = "trans";
    }
    else {
        partitioner = "checkerboard";
    }
    options.set<std::string>( "partitioner", partitioner );

    // Options for for periodic grids
    options.set<bool>( "periodic_x", false );
    options.set<bool>( "periodic_y", false );
}

void RegularMeshGenerator::generate( const Grid& grid, Mesh& mesh ) const {
    ATLAS_ASSERT( !mesh.generated() );

    const RegularGrid rg = RegularGrid( grid );
    if ( !rg ) {
        throw_Exception( "RegularMeshGenerator can only work with a Regular grid", Here() );
    }

    size_t nb_parts = options.get<size_t>( "nb_parts" );

    std::string partitioner_type = "checkerboard";
    options.get( "checkerboard", partitioner_type );

    // if ( rg->nlat()%2 == 1 ) partitioner_factory = "equal_regions"; // Odd
    // number of latitudes
    // if ( nb_parts == 1 || eckit::mpi::size() == 1 ) partitioner_factory =
    // "equal_regions"; // Only one part --> Trans is slower

    grid::Partitioner partitioner( partitioner_type, nb_parts );
    grid::Distribution distribution( partitioner.partition( grid ) );
    generate( grid, distribution, mesh );
}

void RegularMeshGenerator::hash( eckit::Hash& h ) const {
    h.add( "RegularMeshGenerator" );
    options.hash( h );
}

void RegularMeshGenerator::generate( const Grid& grid, const grid::Distribution& distribution, Mesh& mesh ) const {
    const auto rg = RegularGrid( grid );
    if ( !rg ) {
        throw_Exception( "Grid could not be cast to a Regular", Here() );
    }

    ATLAS_ASSERT( !mesh.generated() );

    if ( grid.size() != static_cast<idx_t>( distribution.size() ) ) {
        std::stringstream msg;
        msg << "Number of points in grid (" << grid.size()
            << ") different from "
               "number of points in grid distribution ("
            << distribution.size() << ")";
        throw_AssertionFailed( msg.str(), Here() );
    }

    // clone some grid properties
    setGrid( mesh, rg, distribution );

    generate_mesh( rg, distribution, mesh );
}

void RegularMeshGenerator::generate_mesh( const RegularGrid& rg, const grid::Distribution& distribution,
                                          // const Region& region,
                                          Mesh& mesh ) const {
    int mypart = options.get<size_t>( "part" );
    int nparts = options.get<size_t>( "nb_parts" );
    int nx     = rg.nx();
    int ny     = rg.ny();

    bool periodic_x = options.get<bool>( "periodic_x" ) or rg.periodic();
    bool periodic_y = options.get<bool>( "periodic_y" );

    Log::debug() << Here() << " periodic_x = " << periodic_x << std::endl;
    Log::debug() << Here() << " periodic_y = " << periodic_y << std::endl;

// for asynchronous output
#if DEBUG_OUTPUT
    sleep( mypart );
#endif

    // this function should do the following:
    // - define nodes with
    //      mesh.nodes().resize(nnodes);
    //      mesh::Nodes& nodes = mesh.nodes();
    //    following properties should be defined:
    //      array::ArrayView<double,2> xy            ( nodes.xy() );
    //      array::ArrayView<gidx_t,1> glb_idx       ( nodes.global_index() );
    //      array::ArrayView<int,   1> part          ( nodes.partition() );
    //      array::ArrayView<int,   1> ghost         ( nodes.ghost() );
    //      array::ArrayView<int,   1> flags         ( nodes.flags() );
    // - define cells (only quadrilaterals for now) with
    //      mesh.cells().add( new mesh::temporary::Quadrilateral(), nquads  );
    //    further define cells with
    //      array::ArrayView<gidx_t,1> cells_glb_idx( mesh.cells().global_index()
    //      );
    //      array::ArrayView<int,1>    cells_part(    mesh.cells().partition() );
    // - define connectivity with
    //      mesh::HybridElements::Connectivity& node_connectivity =
    //      mesh.cells().node_connectivity();
    //      node_connectivity.set( jcell, quad_nodes );
    //    where quad_nodes is a 4-element integer array containing the LOCAL
    //    indices of the nodes

    // Start with calculating number of quadrilaterals
    // The rule do determine if a cell belongs to a proc is the following: if the
    // lowerleft corner of the cell belongs to that proc.
    // so we loop over all gridpoints, select those that belong to the proc, and
    // determine the number of cells
    int ii_glb;  // global index
    int ncells;

    // vector of local indices: necessary for remote indices of ghost nodes
    std::vector<int> local_idx( rg.size(), -1 );
    std::vector<int> current_idx( nparts, 0 );  // index counter for each proc

    // determine rectangle (ix_min:ix_max) x (iy_min:iy_max) surrounding the nodes
    // on this processor
    int ix_min, ix_max, iy_min, iy_max, ix_glb, iy_glb, ix, iy;
    int nnodes_nonghost, nnodes;  // number of nodes: non-ghost; total;  inside
                                  // surrounding rectangle
    int nnodes_SR, ii;

    // loop over all points to determine local indices and surroundig rectangle
    ix_min          = nx + 1;
    ix_max          = 0;
    iy_min          = ny + 1;
    iy_max          = 0;
    nnodes_nonghost = 0;

    ii_glb = 0;
    for ( iy = 0; iy < ny; iy++ ) {
        for ( ix = 0; ix < nx; ix++ ) {
            local_idx[ii_glb] = current_idx[distribution.partition( ii_glb )]++;  // store local index on
                                                                                  // the local proc of
                                                                                  // this point
            if ( distribution.partition( ii_glb ) == mypart ) {
                ++nnodes_nonghost;  // non-ghost node: belongs to this part
                ix_min = std::min( ix_min, ix );
                ix_max = std::max( ix_max, ix );
                iy_min = std::min( iy_min, iy );
                iy_max = std::max( iy_max, iy );
            }
            ++ii_glb;  // global index
        }
    }

    // add one row/column for ghost nodes (which include periodicity points)
    ix_max = ix_max + 1;
    iy_max = iy_max + 1;

#if DEBUG_OUTPUT_DETAIL
    std::cout << "[" << mypart << "] : "
              << "SR = " << ix_min << ":" << ix_max << " x " << iy_min << ":" << iy_max << std::endl;
#endif

    // dimensions of surrounding rectangle (SR)
    int nxl = ix_max - ix_min + 1;
    int nyl = iy_max - iy_min + 1;

    // upper estimate for number of nodes
    nnodes_SR = nxl * nyl;

    // partitions and local indices in SR
    std::vector<int> parts_SR( nnodes_SR, -1 );
    std::vector<int> local_idx_SR( nnodes_SR, -1 );
    std::vector<bool> is_ghost_SR( nnodes_SR, true );
    ii = 0;  // index inside SR
    for ( iy = 0; iy < nyl; iy++ ) {
        iy_glb = ( iy_min + iy );  // global y-index
        for ( ix = 0; ix < nxl; ix++ ) {
            ix_glb          = ( ix_min + ix );  // global x-index
            is_ghost_SR[ii] = !( ( parts_SR[ii] == mypart ) && ix < nxl - 1 && iy < nyl - 1 );
            if ( ix_glb < nx && iy_glb < ny ) {
                ii_glb           = (iy_glb)*nx + ix_glb;  // global index
                parts_SR[ii]     = distribution.partition( ii_glb );
                local_idx_SR[ii] = local_idx[ii_glb];
                is_ghost_SR[ii]  = !( ( parts_SR[ii] == mypart ) && ix < nxl - 1 && iy < nyl - 1 );
            }
            else if ( ix_glb == nx && iy_glb < ny ) {
                // take properties from the point to the left
                parts_SR[ii]     = distribution.partition( iy_glb * nx + ix_glb - 1 );
                local_idx_SR[ii] = -1;
                is_ghost_SR[ii]  = true;
            }
            else if ( iy_glb == ny && ix_glb < nx ) {
                // take properties from the point below
                parts_SR[ii]     = distribution.partition( ( iy_glb - 1 ) * nx + ix_glb );
                local_idx_SR[ii] = -1;
                is_ghost_SR[ii]  = true;
            }
            else {
                // take properties from the point belowleft
                parts_SR[ii]     = distribution.partition( ( iy_glb - 1 ) * nx + ix_glb - 1 );
                local_idx_SR[ii] = -1;
                is_ghost_SR[ii]  = true;
            }
            ++ii;
        }
    }

#if DEBUG_OUTPUT_DETAIL
    std::cout << "[" << mypart << "] : "
              << "parts_SR = ";
    for ( ii = 0; ii < nnodes_SR; ii++ )
        std::cout << parts_SR[ii] << ",";
    std::cout << std::endl;
    std::cout << "[" << mypart << "] : "
              << "local_idx_SR = ";
    for ( ii = 0; ii < nnodes_SR; ii++ )
        std::cout << local_idx_SR[ii] << ",";
    std::cout << std::endl;
    std::cout << "[" << mypart << "] : "
              << "is_ghost_SR = ";
    for ( ii = 0; ii < nnodes_SR; ii++ )
        std::cout << is_ghost_SR[ii] << ",";
    std::cout << std::endl;
#endif

    // vectors marking nodes that are necessary for this proc's cells
    std::vector<bool> is_node_SR( nnodes_SR, false );

    // determine number of cells and number of nodes
    nnodes = 0;
    ncells = 0;
    for ( iy = 0; iy < nyl - 1; iy++ ) {      // don't loop into ghost/periodicity row
        for ( ix = 0; ix < nxl - 1; ix++ ) {  // don't loop into ghost/periodicity column
            ii = iy * nxl + ix;
            if ( !is_ghost_SR[ii] ) {
                // mark this node as being used
                if ( !is_node_SR[ii] ) {
                    ++nnodes;
                    is_node_SR[ii] = true;
                }
                // check if this node is the lowerleft corner of a new cell
                if ( ( ix_min + ix < nx - 1 || periodic_x ) && ( iy_min + iy < ny - 1 || periodic_y ) ) {
                    ++ncells;
                    // mark lowerright corner
                    ii = iy * nxl + ix + 1;
                    if ( !is_node_SR[ii] ) {
                        ++nnodes;
                        is_node_SR[ii] = true;
                    }
                    // mark upperleft corner
                    ii = ( iy + 1 ) * nxl + ix;
                    if ( !is_node_SR[ii] ) {
                        ++nnodes;
                        is_node_SR[ii] = true;
                    }
                    // mark upperright corner
                    ii = ( iy + 1 ) * nxl + ix + 1;
                    if ( !is_node_SR[ii] ) {
                        ++nnodes;
                        is_node_SR[ii] = true;
                    }
                }
                // periodic points are always needed, even if they don't belong to a
                // cell
                ii = iy * nxl + ix + 1;
                if ( periodic_x && ix_min + ix == nx - 1 && !is_node_SR[ii] ) {
                    ++nnodes;
                    is_node_SR[ii] = true;
                }
                ii = ( iy + 1 ) * nxl + ix;
                if ( periodic_y && iy_min + iy == ny - 1 && !is_node_SR[ii] ) {
                    ++nnodes;
                    is_node_SR[ii] = true;
                }
                ii = ( iy + 1 ) * nxl + ix + 1;
                if ( periodic_x && periodic_y && ix_min + ix == nx - 1 && iy_min + iy == ny - 1 && !is_node_SR[ii] ) {
                    ++nnodes;
                    is_node_SR[ii] = true;
                }
            }
        }
    }

#if DEBUG_OUTPUT_DETAIL
    std::cout << "[" << mypart << "] : "
              << "nnodes = " << nnodes << std::endl;
    std::cout << "[" << mypart << "] : "
              << "is_node_SR = ";
    for ( ii = 0; ii < nnodes_SR; ii++ )
        std::cout << is_node_SR[ii] << ",";
    std::cout << std::endl;
#endif

    // define nodes and associated properties
    mesh.nodes().resize( nnodes );
    mesh::Nodes& nodes = mesh.nodes();
    auto xy            = array::make_view<double, 2>( nodes.xy() );
    auto lonlat        = array::make_view<double, 2>( nodes.lonlat() );
    auto glb_idx       = array::make_view<gidx_t, 1>( nodes.global_index() );
    auto remote_idx    = array::make_indexview<idx_t, 1>( nodes.remote_index() );
    auto part          = array::make_view<int, 1>( nodes.partition() );
    auto ghost         = array::make_view<int, 1>( nodes.ghost() );
    auto flags         = array::make_view<int, 1>( nodes.flags() );

    // define cells and associated properties
    mesh.cells().add( new mesh::temporary::Quadrilateral(), ncells );
    int quad_begin                                        = mesh.cells().elements( 0 ).begin();
    auto cells_part                                       = array::make_view<int, 1>( mesh.cells().partition() );
    mesh::HybridElements::Connectivity& node_connectivity = mesh.cells().node_connectivity();

    idx_t quad_nodes[4];
    int jcell = quad_begin;
    int inode, inode_nonghost, inode_ghost;

    // global indices for periodicity points
    inode = nx * ny;
    std::vector<int> glb_idx_px( ny + 1, -1 );
    std::vector<int> glb_idx_py( nx + 1, -1 );
    if ( periodic_x ) {
        for ( iy = 0; iy < ny + ( periodic_y ? 1 : 0 ); iy++ ) {
            glb_idx_px[iy] = inode++;
        }
    }
    if ( periodic_y ) {
        for ( ix = 0; ix < nx; ix++ ) {
            glb_idx_py[ix] = inode++;
        }
    }

    // loop over nodes and set properties
    ii             = 0;
    inode_nonghost = 0;
    inode_ghost    = nnodes_nonghost;  // ghost nodes start counting after nonghost nodes
    for ( iy = 0; iy < nyl; iy++ ) {
        for ( ix = 0; ix < nxl; ix++ ) {
            // node properties
            if ( is_node_SR[ii] ) {
                // set node counter
                if ( is_ghost_SR[ii] ) {
                    inode = inode_ghost++;
                }
                else {
                    inode = inode_nonghost++;
                }
                // global index
                ix_glb = ( ix_min + ix );  // don't take modulus here: periodicity points
                                           // have their own global index.
                iy_glb = ( iy_min + iy );
                if ( ix_glb < nx && iy_glb < ny ) {
                    ii_glb = iy_glb * nx + ix_glb;  // no periodic point
                }
                else {
                    if ( ix_glb == nx ) {
                        // periodicity point in x-direction
                        ii_glb = glb_idx_px[iy_glb];
                    }
                    else {
                        // periodicity point in x-direction
                        ii_glb = glb_idx_py[ix_glb];
                    }
                }
                glb_idx( inode ) = ii_glb + 1;  // starting from 1
                // grid coordinates
                double _xy[2];
                if ( iy_glb < ny ) {
                    // normal calculation
                    rg.xy( ix_glb, iy_glb, _xy );
                }
                else {
                    // for periodic_y grids, iy_glb==ny lies outside the range of
                    // latitudes in the Structured grid...
                    // so we extrapolate from two other points -- this is okay for regular
                    // grids with uniform spacing.
                    double xy1[2], xy2[2];
                    rg.xy( ix_glb, iy_glb - 1, xy1 );
                    rg.xy( ix_glb, iy_glb - 2, xy2 );
                    _xy[0] = 2 * xy1[0] - xy2[0];
                    _xy[1] = 2 * xy1[1] - xy2[1];
                }
                xy( inode, LON ) = _xy[LON];
                xy( inode, LAT ) = _xy[LAT];

                // geographic coordinates by using projection
                rg.projection().xy2lonlat( _xy );
                lonlat( inode, LON ) = _xy[LON];
                lonlat( inode, LAT ) = _xy[LAT];

                // part
                part( inode ) = parts_SR[ii];
                // ghost nodes
                ghost( inode ) = is_ghost_SR[ii];
                // flags
                Topology::reset( flags( inode ) );
                if ( ghost( inode ) ) {
                    Topology::set( flags( inode ), Topology::GHOST );
                    remote_idx( inode ) = local_idx_SR[ii];
                    // change local index -- required for cells
                    local_idx_SR[ii] = inode;
                }
                else {
                    remote_idx( inode ) = -1;
                }

#if DEBUG_OUTPUT_DETAIL
                std::cout << "[" << mypart << "] : "
                          << "New node "
                          << "\n\t";
                std::cout << "[" << mypart << "] : "
                          << "\tinode=" << inode << "; ix_glb=" << ix_glb << "; iy_glb=" << iy_glb
                          << "; glb_idx=" << ii_glb << std::endl;
                std::cout << "[" << mypart << "] : "
                          << "\tglon=" << lonlat( inode, 0 ) << "; glat=" << lonlat( inode, 1 )
                          << "; glb_idx=" << glb_idx( inode ) << std::endl;
#endif
            }
            ++ii;
        }
    }

    // loop over nodes and define cells
    for ( iy = 0; iy < nyl - 1; iy++ ) {      // don't loop into ghost/periodicity row
        for ( ix = 0; ix < nxl - 1; ix++ ) {  // don't loop into ghost/periodicity column
            ii = iy * nxl + ix;
            if ( !is_ghost_SR[ii] ) {
                if ( ( ix_min + ix < nx - 1 || periodic_x ) && ( iy_min + iy < ny - 1 || periodic_y ) ) {
                    // define cell corners (local indices)
                    quad_nodes[0] = local_idx_SR[ii];
                    quad_nodes[3] = local_idx_SR[iy * nxl + ix + 1];          // point to the right
                    quad_nodes[2] = local_idx_SR[( iy + 1 ) * nxl + ix + 1];  // point above right
                    quad_nodes[1] = local_idx_SR[( iy + 1 ) * nxl + ix];      // point above
                    node_connectivity.set( jcell, quad_nodes );
                    cells_part( jcell ) = mypart;
#if DEBUG_OUTPUT_DETAIL
                    std::cout << "[" << mypart << "] : "
                              << "New quad " << jcell << "\n\t";
                    std::cout << "[" << mypart << "] : " << quad_nodes[0] << "," << quad_nodes[1] << ","
                              << quad_nodes[2] << "," << quad_nodes[3] << std::endl;
#endif
                    ++jcell;
                }
            }
        }
    }

#if DEBUG_OUTPUT
    // list nodes
    for ( inode = 0; inode < nnodes; inode++ ) {
        std::cout << "[" << mypart << "] : "
                  << " node " << inode << ": ghost = " << ghost( inode ) << ", glb_idx = " << glb_idx( inode ) - 1
                  << ", part = " << part( inode ) << ", lon = " << lonlat( inode, 0 )
                  << ", lat = " << lonlat( inode, 1 ) << ", remote_idx = " << remote_idx( inode ) << std::endl;
    }

    int* cell_nodes;
    for ( jcell = 0; jcell < ncells; jcell++ ) {
        std::cout << "[" << mypart << "] : "
                  << " cell " << jcell << ": " << node_connectivity( jcell, 0 ) << "," << node_connectivity( jcell, 1 )
                  << "," << node_connectivity( jcell, 2 ) << "," << node_connectivity( jcell, 3 ) << std::endl;
    }
#endif

    generateGlobalElementNumbering( mesh );

    nodes.metadata().set( "parallel", true );
}

namespace {
static MeshGeneratorBuilder<RegularMeshGenerator> __RegularMeshGenerator( "regular" );
}

}  // namespace meshgenerator
}  // namespace atlas
