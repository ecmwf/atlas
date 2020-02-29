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
#include <sstream>

#include "atlas/parallel/mpi/mpi.h"

#include "atlas/array.h"
#include "atlas/field/Field.h"
#include "atlas/grid.h"
#include "atlas/grid/detail/partitioner/EqualRegionsPartitioner.h"
#include "atlas/library/config.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/mesh/Nodes.h"
#include "atlas/mesh/actions/BuildParallelFields.h"
#include "atlas/mesh/actions/BuildPeriodicBoundaries.h"
#include "atlas/meshgenerator.h"
#include "atlas/output/Gmsh.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/util/CoordinateEnums.h"
#include "atlas/util/Metadata.h"

#include "tests/AtlasTestEnvironment.h"

using Topology = atlas::mesh::Nodes::Topology;

using namespace atlas::array;
using namespace atlas::output;
using namespace atlas::meshgenerator;

namespace atlas {
namespace test {

//-----------------------------------------------------------------------------

class IsGhost {
public:
    IsGhost( const mesh::Nodes& nodes ) :
        part_( make_view<int, 1>( nodes.partition() ) ),
        ridx_( make_indexview<idx_t, 1>( nodes.remote_index() ) ),
        mypart_( mpi::comm().rank() ) {}

    bool operator()( idx_t idx ) const {
        if ( part_( idx ) != mypart_ ) {
            return true;
        }
        if ( ridx_( idx ) != idx ) {
            return true;
        }
        return false;
    }

private:
    array::ArrayView<const int, 1> part_;
    IndexView<const idx_t, 1> ridx_;
    int mypart_;
};

#define DISABLE if ( 0 )
#define ENABLE if ( 1 )

//-----------------------------------------------------------------------------

CASE( "test1" ) {
    Mesh m;

    mesh::Nodes& nodes = m.nodes();
    nodes.resize( 10 );
    auto xy      = make_view<double, 2>( nodes.xy() );
    auto lonlat  = make_view<double, 2>( nodes.lonlat() );
    auto glb_idx = make_view<gidx_t, 1>( nodes.global_index() );
    auto part    = make_view<int, 1>( nodes.partition() );
    auto flags   = make_view<int, 1>( nodes.flags() );
    flags.assign( Topology::NONE );

    // This is typically available
    glb_idx( 0 ) = 1;
    part( 0 )    = 0;
    glb_idx( 1 ) = 2;
    part( 1 )    = 0;
    glb_idx( 2 ) = 3;
    part( 2 )    = 0;
    glb_idx( 3 ) = 4;
    part( 3 )    = 0;
    glb_idx( 4 ) = 5;
    part( 4 )    = 0;
    glb_idx( 5 ) = 6;
    part( 5 )    = std::min( 1, (int)mpi::comm().size() - 1 );
    glb_idx( 6 ) = 7;
    part( 6 )    = std::min( 1, (int)mpi::comm().size() - 1 );
    glb_idx( 7 ) = 8;
    part( 7 )    = std::min( 1, (int)mpi::comm().size() - 1 );
    glb_idx( 8 ) = 9;
    part( 8 )    = std::min( 1, (int)mpi::comm().size() - 1 );
    glb_idx( 9 ) = 10;
    part( 9 )    = std::min( 1, (int)mpi::comm().size() - 1 );

    xy( 0, XX ) = 0.;
    xy( 0, YY ) = 80.;
    Topology::set( flags( 0 ), Topology::BC | Topology::WEST );
    xy( 1, XX ) = 0.;
    xy( 1, YY ) = -80.;
    Topology::set( flags( 1 ), Topology::BC | Topology::WEST );
    xy( 2, XX ) = 90.;
    xy( 2, YY ) = 80.;
    xy( 3, XX ) = 90.;
    xy( 3, YY ) = -80.;
    xy( 4, XX ) = 180.;
    xy( 4, YY ) = 80.;
    xy( 5, XX ) = 180.;
    xy( 5, YY ) = -80.;
    xy( 6, XX ) = 270.;
    xy( 6, YY ) = 80.;
    xy( 7, XX ) = 270.;
    xy( 7, YY ) = -80.;
    xy( 8, XX ) = 360.;
    xy( 8, YY ) = 80.;
    Topology::set( flags( 8 ), Topology::BC | Topology::EAST );
    xy( 9, XX ) = 360.;
    xy( 9, YY ) = -80.;
    Topology::set( flags( 9 ), Topology::BC | Topology::EAST );

    for ( idx_t n = 0; n < xy.shape( 0 ); ++n ) {
        lonlat( n, LON ) = xy( n, XX );
        lonlat( n, LAT ) = xy( n, YY );
    }

    mesh::actions::build_parallel_fields( m );

    EXPECT( nodes.has_field( "remote_idx" ) );

    auto loc = make_indexview<idx_t, 1>( nodes.remote_index() );
    EXPECT( loc( 0 ) == 0 );
    EXPECT( loc( 1 ) == 1 );
    EXPECT( loc( 2 ) == 2 );
    EXPECT( loc( 3 ) == 3 );
    EXPECT( loc( 4 ) == 4 );
    EXPECT( loc( 5 ) == 5 );
    EXPECT( loc( 6 ) == 6 );
    EXPECT( loc( 7 ) == 7 );
    EXPECT( loc( 8 ) == 8 );
    EXPECT( loc( 9 ) == 9 );

    test::IsGhost is_ghost( m.nodes() );

    switch ( mpi::comm().rank() ) {
        case 0:
            EXPECT( is_ghost( 0 ) == false );
            EXPECT( is_ghost( 1 ) == false );
            EXPECT( is_ghost( 2 ) == false );
            EXPECT( is_ghost( 3 ) == false );
            EXPECT( is_ghost( 4 ) == false );
            EXPECT( is_ghost( 5 ) == true );
            EXPECT( is_ghost( 6 ) == true );
            EXPECT( is_ghost( 7 ) == true );
            EXPECT( is_ghost( 8 ) == true );
            EXPECT( is_ghost( 9 ) == true );
            break;
        case 1:
            EXPECT( is_ghost( 0 ) == true );
            EXPECT( is_ghost( 1 ) == true );
            EXPECT( is_ghost( 2 ) == true );
            EXPECT( is_ghost( 3 ) == true );
            EXPECT( is_ghost( 4 ) == true );
            EXPECT( is_ghost( 5 ) == false );
            EXPECT( is_ghost( 6 ) == false );
            EXPECT( is_ghost( 7 ) == false );
            EXPECT( is_ghost( 8 ) == false );
            EXPECT( is_ghost( 9 ) == false );
            break;
    }

    mesh::actions::build_periodic_boundaries( m );

    // points 8 and 9 are periodic slave points of points 0 and 1
    EXPECT( part( 8 ) == 0 );
    EXPECT( part( 9 ) == 0 );
    EXPECT( loc( 8 ) == 0 );
    EXPECT( loc( 9 ) == 1 );
    if ( mpi::comm().rank() == 1 ) {
        EXPECT( is_ghost( 8 ) == true );
        EXPECT( is_ghost( 9 ) == true );
    }
}

CASE( "test2" ) {
    util::Config meshgen_options;
    meshgen_options.set( "angle", 27.5 );
    meshgen_options.set( "triangulate", false );
    meshgen_options.set( "partitioner", "equal_regions" );
    StructuredMeshGenerator generate( meshgen_options );
    Mesh m = generate( Grid( "N32" ) );
    mesh::actions::build_parallel_fields( m );

    mesh::Nodes& nodes = m.nodes();

    test::IsGhost is_ghost( nodes );

    idx_t nb_ghost = 0;
    for ( idx_t jnode = 0; jnode < nodes.size(); ++jnode ) {
        if ( is_ghost( jnode ) ) {
            ++nb_ghost;
        }
    }

    ATLAS_DEBUG_VAR( nb_ghost );
    if ( mpi::comm().rank() == 0 ) {
        EXPECT( nb_ghost == 128 );  // South boundary of Northern hemisphere
    }
    if ( mpi::comm().rank() == 1 ) {
        EXPECT( nb_ghost == 0 );  // Southern hemisphere has no ghosts
    }

    mesh::actions::build_periodic_boundaries( m );

    int nb_periodic = -nb_ghost;
    for ( idx_t jnode = 0; jnode < nodes.size(); ++jnode ) {
        if ( is_ghost( jnode ) ) {
            ++nb_periodic;
        }
    }

    ATLAS_DEBUG_VAR( nb_periodic );

    if ( mpi::comm().rank() == 0 ) {
        EXPECT( nb_periodic == 33 );  // Periodic East boundary of Northern hemisphere
    }
    // (plus one point south)
    if ( mpi::comm().rank() == 1 ) {
        EXPECT( nb_periodic == 32 );  // Periodic East boundary of Southern hemisphere
    }

    Gmsh( "periodic.msh", util::Config( "info", true ) ).write( m );
}

//-----------------------------------------------------------------------------

}  // namespace test
}  // namespace atlas

int main( int argc, char** argv ) {
    return atlas::test::run( argc, argv );
}
