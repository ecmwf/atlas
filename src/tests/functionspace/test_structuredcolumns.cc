/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "eckit/log/Bytes.h"
#include "eckit/types/Types.h"

#include "atlas/array/ArrayView.h"
#include "atlas/array/MakeView.h"
#include "atlas/field/Field.h"
#include "atlas/functionspace/NodeColumns.h"
#include "atlas/functionspace/StructuredColumns.h"
#include "atlas/grid/Partitioner.h"
#include "atlas/grid/StructuredGrid.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/meshgenerator.h"
#include "atlas/output/Gmsh.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/util/CoordinateEnums.h"
#include "atlas/util/MicroDeg.h"

#include "tests/AtlasTestEnvironment.h"

using namespace eckit;
using namespace atlas::functionspace;
using namespace atlas::util;

namespace atlas {
namespace test {

//-----------------------------------------------------------------------------

CASE( "test_functionspace_StructuredColumns_no_halo" ) {
    size_t root          = 0;
    std::string gridname = eckit::Resource<std::string>( "--grid", "O8" );
    Grid grid( gridname );
    util::Config config;
    config.set( "halo", 0 );
    config.set( "periodic_points", true );
    functionspace::StructuredColumns fs( grid, grid::Partitioner( "equal_regions" ), config );
    ATLAS_DEBUG_VAR( fs.size() );
    ATLAS_DEBUG_VAR( eckit::Bytes( fs.footprint() ) );

    Field field     = fs.createField<double>( option::name( "field" ) );
    Field field_glb = fs.createField<double>( option::name( "field_global" ) | option::global( root ) );

    auto value     = array::make_view<double, 1>( field );
    auto value_glb = array::make_view<double, 1>( field_glb );

    value.assign( mpi::comm().rank() );

    fs.gather( field, field_glb );

    Log::info() << "field checksum = " << fs.checksum( field ) << std::endl;

    //  for( size_t j=0; j<value_glb.size(); ++j )
    //    Log::info() << value_glb(j) << " ";
    //  Log::info() << std::endl;

    if ( mpi::comm().rank() == root && mpi::comm().size() == 5 ) {
        std::vector<double> check{
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3,
            3, 3, 3, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3,
            3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2,
            2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
            1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
            1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3,
            3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2,
            2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
            1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 1, 1, 1, 1,
            1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4,
            4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
            4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
            4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4};

        EXPECT( value_glb.size() == idx_t( check.size() ) );
        for ( idx_t j = 0; j < value_glb.size(); ++j ) {
            EXPECT( value_glb( j ) == check[j] );
        }
    }
    ATLAS_TRACE_SCOPE( "output gmsh" ) {
        output::Gmsh gmsh( "structured.msh" );

        gmsh.write( MeshGenerator( "structured" ).generate( grid ) );
        gmsh.write( field );
    }
}

CASE( "test_functionspace_StructuredColumns_halo with output" ) {
    ATLAS_DEBUG_VAR( mpi::comm().size() );
    //  grid::StructuredGrid grid(
    //      grid::StructuredGrid::XSpace( {0.,360.} , {2,4,6,6,4,2} , false ),
    //      grid::StructuredGrid::YSpace( grid::LinearSpacing( {80.,-80.}, 6 ) ),
    //      Projection(),
    //      Domain() );

    std::string gridname = eckit::Resource<std::string>( "--grid", "O8" );

    StructuredGrid grid( gridname );

    int halo = eckit::Resource<int>( "--halo", 2 );
    util::Config config;
    config.set( "halo", halo );
    config.set( "periodic_points", true );
    functionspace::StructuredColumns fs( grid, grid::Partitioner( "equal_regions" ), config );

    Field field = fs.createField<long>( option::name( "field" ) );

    auto value = array::make_view<long, 1>( field );
    auto xy    = array::make_view<double, 2>( fs.xy() );
    auto g     = array::make_view<gidx_t, 1>( fs.global_index() );
    auto r     = array::make_view<idx_t, 1>( fs.remote_index() );
    auto p     = array::make_view<int, 1>( fs.partition() );

    for ( idx_t j = fs.j_begin(); j < fs.j_end(); ++j ) {
        for ( idx_t i = fs.i_begin( j ); i < fs.i_end( j ); ++i ) {
            idx_t n    = fs.index( i, j );
            value( n ) = util::microdeg( xy( n, XX ) );
        }
    }

    // EXPECT( fs.checksum(field) == "cef2694016492d408fa157b7c59ce741" );

    fs.haloExchange( field );

    // EXPECT( fs.checksum(field) == "cef2694016492d408fa157b7c59ce741" );

    ATLAS_TRACE_SCOPE( "Output python" ) {
        eckit::PathName filepath( "test_functionspace_StructuredColumns_halo_p" + std::to_string( mpi::comm().rank() ) +
                                  ".py" );

        std::ofstream f( filepath.asString().c_str(), std::ios::trunc );

        f << "\n"
             "import matplotlib.pyplot as plt"
             "\n"
             "from matplotlib.path import Path"
             "\n"
             "import matplotlib.patches as patches"
             "\n"
             ""
             "\n"
             "from itertools import cycle"
             "\n"
             "import matplotlib.cm as cm"
             "\n"
             "import numpy as np"
             "\n"
             ""
             "\n"
             "fig = plt.figure(figsize=(20,10))"
             "\n"
             "ax = fig.add_subplot(111,aspect='equal')"
             "\n"
             "";

        double xmin = std::numeric_limits<double>::max();
        double xmax = -std::numeric_limits<double>::max();
        double ymin = std::numeric_limits<double>::max();
        double ymax = -std::numeric_limits<double>::max();
        f << "\n"
             "x = [";
        for ( idx_t j = fs.j_begin_halo(); j < fs.j_end_halo(); ++j ) {
            for ( idx_t i = fs.i_begin_halo( j ); i < fs.i_end_halo( j ); ++i ) {
                idx_t n = fs.index( i, j );
                f << xy( n, XX ) << ", ";
                xmin = std::min( xmin, xy( n, XX ) );
                xmax = std::max( xmax, xy( n, XX ) );
            }
        }
        f << "]";

        f << "\n"
             "y = [";
        for ( idx_t j = fs.j_begin_halo(); j < fs.j_end_halo(); ++j ) {
            for ( idx_t i = fs.i_begin_halo( j ); i < fs.i_end_halo( j ); ++i ) {
                idx_t n = fs.index( i, j );
                f << xy( n, YY ) << ", ";
                ymin = std::min( ymin, xy( n, YY ) );
                ymax = std::max( ymax, xy( n, YY ) );
            }
        }
        f << "]";

        f << "\n"
             "g = [";
        for ( idx_t j = fs.j_begin_halo(); j < fs.j_end_halo(); ++j ) {
            for ( idx_t i = fs.i_begin_halo( j ); i < fs.i_end_halo( j ); ++i ) {
                idx_t n = fs.index( i, j );
                f << g( n ) << ", ";
            }
        }
        f << "]";

        f << "\n"
             "p = [";
        for ( idx_t j = fs.j_begin_halo(); j < fs.j_end_halo(); ++j ) {
            for ( idx_t i = fs.i_begin_halo( j ); i < fs.i_end_halo( j ); ++i ) {
                idx_t n = fs.index( i, j );
                f << p( n ) << ", ";
            }
        }
        f << "]";

        f << "\n"
             "r = [";
        for ( idx_t j = fs.j_begin_halo(); j < fs.j_end_halo(); ++j ) {
            for ( idx_t i = fs.i_begin_halo( j ); i < fs.i_end_halo( j ); ++i ) {
                idx_t n = fs.index( i, j );
                f << r( n ) << ", ";
            }
        }
        f << "]";

        f << "\n"
             ""
             "\n"
             "c = [ cm.Paired( float(pp%13)/12. ) for pp in p ]"
             "\n"
             "ax.scatter(x, y, color=c, marker='o')"
             "\n"
             "for i in range("
          << fs.size()
          << "):"
             "\n"
             "  ax.annotate(g[i], (x[i],y[i]), fontsize=8)"
             "\n"
             "";
        f << "\n"
             "ax.set_xlim( "
          << std::min( 0., xmin ) << "-5, " << std::max( 360., xmax )
          << "+5)"
             "\n"
             "ax.set_ylim( "
          << std::min( -90., ymin ) << "-5, " << std::max( 90., ymax )
          << "+5)"
             "\n"
             "ax.set_xticks([0,45,90,135,180,225,270,315,360])"
             "\n"
             "ax.set_yticks([-90,-45,0,45,90])"
             "\n"
             "plt.grid()"
             "\n"
             "plt.show()"
             "\n";
    }
}

//-----------------------------------------------------------------------------

CASE( "test_functionspace_StructuredColumns_halo checks without output" ) {
    std::string gridname = eckit::Resource<std::string>( "--grid", "O8" );

    StructuredGrid grid( gridname );

    int halo = eckit::Resource<int>( "--halo", 2 );
    util::Config config;
    config.set( "halo", halo );
    config.set( "levels", 10 );
    config.set( "periodic_points", true );
    functionspace::StructuredColumns fs( grid, grid::Partitioner( "equal_regions" ), config );

    Field field = fs.createField<long>( option::name( "field" ) );

    auto value = array::make_view<long, 2>( field );
    auto xy    = array::make_view<double, 2>( fs.xy() );

    for ( idx_t j = fs.j_begin(); j < fs.j_end(); ++j ) {
        for ( idx_t i = fs.i_begin( j ); i < fs.i_end( j ); ++i ) {
            idx_t n = fs.index( i, j );
            for ( idx_t k = 0; k < fs.levels(); ++k ) {
                value( n, k ) = util::microdeg( xy( n, XX ) );
            }
        }
    }

    ATLAS_TRACE_SCOPE( "control each value " )
    fs.parallel_for( [&]( idx_t n, idx_t k ) { EXPECT( value( n, k ) == util::microdeg( xy( n, XX ) ) ); } );
    fs.parallel_for(
        [&]( idx_t n, idx_t i, idx_t j, idx_t k ) { EXPECT( value( n, k ) == util::microdeg( grid.x( i, j ) ) ); } );

    Field fieldg = fs.createField( field, option::global() );
    fs.gather( field, fieldg );

    ATLAS_TRACE_SCOPE( "control_global" ) {
        auto valueg = array::make_view<long, 2>( fieldg );
        fs.parallel_for( option::global(), [&]( idx_t n, idx_t i, idx_t j, idx_t k ) {
            EXPECT( valueg( n, k ) == util::microdeg( grid.x( i, j ) ) );
        } );
    }
}

//-----------------------------------------------------------------------------

CASE( "test_functionspace_StructuredColumns halo exchange registration" ) {
    // Test by observing log with ATLAS_DEBUG=1
    // The HaloExchange Cache should be created twice, already found 4 times,
    // erased twice.

    std::string gridname = eckit::Resource<std::string>( "--grid", "O8" );

    StructuredGrid grid( gridname );

    util::Config config;
    config.set( "levels", 10 );
    config.set( "periodic_points", true );
    for ( idx_t i = 0; i < 3; ++i ) {
        config.set( "halo", 2 );
        functionspace::StructuredColumns fs1( grid, grid::Partitioner( "equal_regions" ), config );
        config.set( "halo", 4 );
        functionspace::StructuredColumns fs2( grid, grid::Partitioner( "equal_regions" ), config );

        Field field1 = fs1.createField<long>( option::name( "field" ) );
        Field field2 = fs2.createField<long>( option::name( "field" ) );

        field1.haloExchange();
        field2.haloExchange();
    }
}

//-----------------------------------------------------------------------------

}  // namespace test
}  // namespace atlas

int main( int argc, char** argv ) {
    return atlas::test::run( argc, argv );
}
