/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "eckit/memory/ScopedPtr.h"
#include "eckit/types/Types.h"

#include "atlas/array/ArrayView.h"
#include "atlas/array/MakeView.h"
#include "atlas/field/Field.h"
#include "atlas/functionspace/NodeColumns.h"
#include "atlas/functionspace/StructuredColumns.h"
#include "atlas/grid/Grid.h"
#include "atlas/library/Library.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/meshgenerator/MeshGenerator.h"
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
    int root = 0;
    Grid grid( "O8" );
    util::Config config;
    config.set( "halo", 0 );
    config.set( "periodic_points", true );
    functionspace::StructuredColumns fs( grid, grid::Partitioner( "equal_regions" ), config );

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

        EXPECT( value_glb.size() == check.size() );
        for ( size_t j = 0; j < value_glb.size(); ++j ) {
            EXPECT( value_glb( j ) == check[j] );
        }
    }

    output::Gmsh gmsh( "structured.msh" );

    gmsh.write( MeshGenerator( "structured" ).generate( grid ) );
    gmsh.write( field );
}

CASE( "test_functionspace_StructuredColumns_halo" ) {
    ATLAS_DEBUG_VAR( mpi::comm().size() );
    //  grid::StructuredGrid grid(
    //      grid::StructuredGrid::XSpace( {0.,360.} , {2,4,6,6,4,2} , false ),
    //      grid::StructuredGrid::YSpace( grid::LinearSpacing( {80.,-80.}, 6 ) ),
    //      Projection(),
    //      Domain() );

    std::string gridname = eckit::Resource<std::string>( "--grid", "O8" );

    grid::StructuredGrid grid( gridname );

    int halo = eckit::Resource<int>( "--halo", 2 );
    util::Config config;
    config.set( "halo", halo );
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

//-----------------------------------------------------------------------------

CASE( "test_functionspace_StructuredColumns_halo" ) {
    std::string gridname = eckit::Resource<std::string>( "--grid", "O8" );

    grid::StructuredGrid grid( gridname );

    int halo = eckit::Resource<int>( "--halo", 2 );
    util::Config config;
    config.set( "halo", halo );
    config.set( "levels", 10 );
    config.set( "periodic_points", true );
    functionspace::StructuredColumns fs( grid, grid::Partitioner( "equal_regions" ), config );
    auto for_ij = fs.for_ij();

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
    for_ij( [=]( idx_t i, idx_t j ) {
        idx_t n = fs.index( i, j );
        for ( idx_t k = 0; k < fs.levels(); ++k ) {
            EXPECT( value( n, k ) == util::microdeg( xy( n, XX ) ) );
        }
    } );

    PointXY dp{180., 45.};
}

//-----------------------------------------------------------------------------


class ComputeVertical {
    std::vector<idx_t> nvaux_;
    idx_t nlev_;
    idx_t nlevaux_;

public:
    ComputeVertical( const std::vector<double>& zcoord ) {
        nlev_          = zcoord.size() - 2;
        double dzcoord = zcoord.back() - zcoord.front();
        ASSERT( dzcoord > 0 );
        for ( idx_t jlev = 0; jlev <= nlev_; ++jlev ) {
            dzcoord = std::min( dzcoord, zcoord[jlev + 1] - zcoord[jlev] );
        }
        nlevaux_ = std::round( 2. / dzcoord + 0.5 ) + 1;
        nvaux_.resize( nlevaux_ + 1 );
        double dzaux = 1. / double( nlevaux_ );

        double zaux = 0.;
        for ( idx_t jlevaux = 0; jlevaux <= nlevaux_; ++jlevaux ) {
            for ( idx_t jlev = 0; jlev <= nlev_; ++jlev ) {
                if ( zaux <= zcoord[jlev + 1] ) {
                    nvaux_[jlevaux] = jlev;
                    break;
                }
            }
            zaux += dzaux;
        }
    }
    idx_t operator()( double z ) { return nvaux_[std::floor( z * nlevaux_ )]; }
};


class ComputeNorth {
    std::vector<double> y_;
    double dy_;
    const double tol{0.5e-6};

public:
    ComputeNorth( const grid::StructuredGrid& grid ) {
        y_.resize( grid.ny() );
        for ( idx_t j = 0; j < grid.ny(); ++j ) {
            y_[j] = grid.y( j );
        }
        dy_ = std::abs( y_[1] - y_[0] );
    }
    idx_t operator()( double y ) {
        idx_t j = std::floor( ( y_[0] - ( y - tol ) ) / dy_ );
#ifndef NDEBUG
        ASSERT( j >= -1 );
#endif
        while ( y_[std::max( j, 0 )] > y - tol ) {
            ++j;
        }
        if ( j >= 0 ) {
            do {
                --j;
            } while ( j >= y_.size() || y_[j] < y - tol );
        }
        return j;
    }
};

class ComputeWest {
    std::vector<double> dx;
    std::vector<double> xref;
    const double tol{0.5e-6};
    idx_t ny;

public:
    ComputeWest( const grid::StructuredGrid& grid ) {
        if ( not grid::RegularGrid( grid ) &&
             std::abs( std::max( std::abs( grid.y().front() ), std::abs( grid.y().back() ) ) - 90. ) < tol ) {
            throw eckit::NotImplemented( "ComputeWest not yet implemented for irregular grids with latitudes at pole",
                                         Here() );
        }
        dx.resize( grid.ny() );
        xref.resize( grid.ny() );
        for ( idx_t j = 0; j < grid.ny(); ++j ) {
            dx[j]   = std::abs( grid.x( 1, j ) - grid.x( 0, j ) );
            xref[j] = grid.x( 0, j );
        }
        ny = grid.ny();
    }
    idx_t operator()( const double& x, idx_t j ) {
        idx_t jj{j};
        if ( jj < 0 ) { jj = -j - 1; }
        if ( jj >= ny ) { jj = ny - 1 - ( jj - ny ); }
        idx_t i = std::floor( ( x + tol - xref[jj] ) / dx[jj] );
        return i;
    }
};


// @class ComputeHorizontalStencil
// @brief Compute stencil in horizontal direction (i,j)
//
// Given a stencil width, the stencil for a given P{x,y} is:
//
//        i[0]     i[1]     i[2]    i[3]
//         x        x        x         x       j + 0
//          x       x       x        x         j + 1
//                     P
//          x       x       x        x         j + 2
//         x        x        x         x       j + 3
//
//   In case the x-component of P is aligned with any
//   stencil, gridpoint, the stencil will assume the grid point
//   is on the point P's left side:
//
//        i[0]     i[1]     i[2]    i[3]
//         x        x        x         x       j + 0
//          x       x       x        x         j + 1
//                  P
//          x       x       x        x         j + 2
//         x        x        x         x       j + 3


class ComputeHorizontalStencil {
    ComputeNorth compute_north_;
    ComputeWest compute_west_;
    idx_t stencil_width_;
    idx_t stencil_begin_;

public:
    ComputeHorizontalStencil( const grid::StructuredGrid& grid, idx_t stencil_width ) :
        compute_north_( grid ),
        compute_west_( grid ),
        stencil_width_( stencil_width ) {
        stencil_begin_ = stencil_width_ - idx_t( double( stencil_width_ ) / 2. + 1. );
    }
    template <typename Vector>
    void operator()( const double& x, const double& y, Vector& i, idx_t& j ) {
        j = compute_north_( y ) - stencil_begin_;
        for ( idx_t jj = 0; jj < stencil_width_; ++jj ) {
            i[jj] = compute_west_( x, j + jj ) - stencil_begin_;
        }
    }
};

CASE( "test_departurepoint" ) {
    if ( mpi::comm().size() == 1 ) {
        std::string gridname = eckit::Resource<std::string>( "--grid", "O8" );

        grid::StructuredGrid grid( gridname );

        int halo = eckit::Resource<int>( "--halo", 2 );
        util::Config config;
        config.set( "halo", halo );
        config.set( "levels", 9 );
        config.set( "periodic_points", true );
        functionspace::StructuredColumns fs( grid, grid::Partitioner( "equal_regions" ), config );

        double tol = 0.5e-6;

        ComputeNorth compute_j_north( grid );
        ComputeWest compute_i_west( grid );

        auto dp = PointXY{180., grid.y( 0 ) + 0.5 * tol};
        dp      = PointXY{0., grid.y( 0 )};
        dp      = PointXY{0., grid.y().back()};
        struct IJ {
            idx_t i, j;
        };
        std::vector<IJ> ij;
        ij.reserve( 4 );
        auto j0 = compute_j_north( dp.y() );
        Log::info() << "j0 = " << j0 << std::endl;
        for ( idx_t j = j0 - 1; j < j0 + 3; ++j ) {
            auto i0 = compute_i_west( dp.x(), j );
            idx_t i = i0 - 1;
            ij.push_back( {i, j} );
        }

        auto glb_idx = array::make_view<gidx_t, 1>( fs.global_index() );
        for ( auto p : ij ) {
            idx_t n = fs.index( p.i, p.j );
            EXPECT( n < glb_idx.shape( 0 ) );
            Log::info() << "dp i,j = " << p.i << ", " << p.j << " glb_idx = " << glb_idx( n ) << std::endl;
        }

        constexpr int stencil_width = 4;

        ComputeHorizontalStencil compute_stencil( grid, stencil_width );
        std::array<idx_t, stencil_width> stencil_i;
        idx_t stencil_j;

        auto departure_points = {
            PointXY( 0., 90. ),
            PointXY( 0., -90. ),
            PointXY( 0., 0. ),
            PointXY( 360., 0. ),
        };
        for ( auto p : departure_points ) {
            Log::info() << p << std::endl;
            compute_stencil( p.x(), p.y(), stencil_i, stencil_j );
            for ( idx_t j = 0; j < stencil_width; ++j ) {
                Log::info() << stencil_i[j] << " " << stencil_j + j << "   --   "
                            << glb_idx( fs.index( stencil_i[j], stencil_j + j ) ) << std::endl;
            }
            Log::info() << std::endl;
        }


        // vertical
        idx_t nlev = fs.levels();
        std::vector<double> zcoord( fs.levels() + 2 );
        double dzcoord = 1. / double( fs.levels() + 1 );
        for ( idx_t jlev = 0; jlev <= nlev + 1; ++jlev ) {
            zcoord[jlev] = jlev * dzcoord;
        }
    }
}

//-----------------------------------------------------------------------------

}  // namespace test
}  // namespace atlas

int main( int argc, char** argv ) {
    return atlas::test::run( argc, argv );
}
