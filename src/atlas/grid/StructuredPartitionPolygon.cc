/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "atlas/grid/StructuredPartitionPolygon.h"

#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>

#include "atlas/array/MakeView.h"
#include "atlas/functionspace/StructuredColumns.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/runtime/Trace.h"
#include "atlas/util/CoordinateEnums.h"

namespace atlas {
namespace grid {

util::Polygon::edge_set_t compute_edges( const functionspace::FunctionSpaceImpl& _fs, idx_t _halo,
                                         std::vector<Point2>& points, std::vector<Point2>& bb ) {
    if ( not dynamic_cast<const functionspace::detail::StructuredColumns*>( &_fs ) ) {
        throw_Exception( "Could not cast functionspace to StructuredColumns", Here() );
    }
    const auto& fs  = dynamic_cast<const functionspace::detail::StructuredColumns&>( _fs );
    const auto grid = fs.grid();
    const auto dom  = RectangularDomain( grid.domain() );

    if ( _halo > 0 ) {
        throw_Exception( "halo must be smaller than that of StructuredColumns", Here() );
    }

    auto equal = []( const double& a, const double& b ) { return std::abs( a - b ) < 1.e-12; };

    util::Polygon::edge_set_t edges;

    PointXY p;
    idx_t c{0};
    idx_t i, j;

    auto add_edge = [&]( idx_t p1, idx_t p2 ) {
        util::Polygon::edge_t edge = {p1, p2};
        edges.insert( edge );
        //Log::info() << edge.first << "  " << edge.second << std::endl;
    };

    auto make_corner = [&]() {
        PointXY ptmp;
        ptmp  = p;
        p[XX] = points.back()[XX];
        p[YY] = 0.5 * ( points.back()[YY] + grid.y( j ) );
        points.emplace_back( p );
        add_edge( c, c + 1 );
        c++;
        p = ptmp;


        ptmp  = p;
        p[YY] = points.back()[YY];
        points.emplace_back( p );
        add_edge( c, c + 1 );
        c++;
        p = ptmp;
    };

    bool debug = true;

    double ymax, ymin, xmax, xmin;

    // Top
    // Top left point
    {
        j = fs.j_begin();
        i = fs.i_begin( j );
        if ( j == 0 ) {
            p[YY] = dom.ymax();
        }
        else {
            p[YY] = 0.5 * ( grid.y( j - 1 ) + grid.y( j ) );
        }
        if ( i == 0 ) {
            p[XX] = dom.xmin();
        }
        else {
            p[XX] = 0.5 * ( grid.x( i - 1, j ) + grid.x( i, j ) );
        }
        points.emplace_back( p );
    }

    // Top right point
    {
        j = fs.j_begin();
        i = fs.i_end( j ) - 1;
        if ( j == 0 ) {
            p[YY] = dom.ymax();
        }
        else {
            p[YY] = 0.5 * ( grid.y( j - 1 ) + grid.y( j ) );
        }
        if ( i == grid.nx( j ) - 1 ) {
            p[XX] = dom.xmax();
        }
        else {
            p[XX] = 0.5 * ( grid.x( i, j ) + grid.x( i + 1, j ) );
        }
        points.emplace_back( p );
        add_edge( c, c + 1 );
        c++;

        ymax = p[YY];
        xmax = p[XX];
    }
    // Right side
    {
        for ( j = fs.j_begin(); j < fs.j_end() - 1; ++j ) {
            p[YY] = grid.y( j );
            i     = fs.i_end( j ) - 1;
            if ( i == grid.nx( j ) - 1 ) {
                p[XX] = dom.xmax();
            }
            else {
                p[XX] = 0.5 * ( grid.x( i, j ) + grid.x( i + 1, j ) );
            }
            if ( j > fs.j_begin() && not equal( p[XX], points.back()[XX] ) ) {
                if ( not debug ) {
                    make_corner();
                }
                else {
                    PointXY ptmp = p;
                    p[XX]        = points.back()[XX];
                    p[YY]        = 0.5 * ( points.back()[YY] + p[YY] );
                    points.emplace_back( p );
                    add_edge( c, c + 1 );
                    c++;
                    p = ptmp;

                    ptmp  = p;
                    p[YY] = points.back()[YY];
                    points.emplace_back( p );
                    add_edge( c, c + 1 );
                    c++;
                    p = ptmp;
                }
            }
            points.emplace_back( p );
            add_edge( c, c + 1 );
            c++;

            xmax = std::min( xmax, p[XX] );
        }
    }
    // Bottom
    // Bottom right point(s)
    {
        j = fs.j_end() - 1;
        i = fs.i_end( j ) - 1;
        if ( j == grid.ny() - 1 ) {
            p[YY] = dom.ymin();
        }
        else {
            p[YY] = 0.5 * ( grid.y( j ) + grid.y( j + 1 ) );
        }
        if ( i == grid.nx( j ) - 1 ) {
            p[XX] = dom.xmax();
        }
        else {
            p[XX] = 0.5 * ( grid.x( i, j ) + grid.x( i + 1, j ) );
        }

        ymin = p[YY];

        PointXY pmin = p;

        if ( not equal( p[XX], points.back()[XX] ) ) {
            if ( not debug ) {
                make_corner();
            }
            else {
                PointXY ptmp;
                ptmp  = p;
                p[XX] = points.back()[XX];
                p[YY] = 0.5 * ( points.back()[YY] + grid.y( j ) );
                points.emplace_back( p );
                add_edge( c, c + 1 );
                c++;
                pmin = p;
                xmax = std::min( xmax, p[XX] );

                p = ptmp;


                ptmp  = p;
                p[YY] = points.back()[YY];
                points.emplace_back( p );
                add_edge( c, c + 1 );
                c++;
                p = ptmp;
            }
        }
        if ( xmax - grid.xspace().dx()[j] < grid.x( i, j ) ) {
            xmax = std::min( xmax, 0.5 * ( grid.x( i + 1, j ) + grid.x( i, j ) ) );
        }
        else {
            ymin = pmin[YY];
        }

        points.emplace_back( p );
        add_edge( c, c + 1 );
        c++;
    }

    // Bottom left point
    {
        j = fs.j_end() - 1;
        i = fs.i_begin( j );
        if ( j == grid.ny() - 1 ) {
            p[YY] = dom.ymin();
        }
        else {
            p[YY] = 0.5 * ( grid.y( j ) + grid.y( j + 1 ) );
        }
        if ( i == 0 ) {
            p[XX] = dom.xmin();
        }
        else {
            p[XX] = 0.5 * ( grid.x( i - 1, j ) + grid.x( i, j ) );
        }
        points.emplace_back( p );
        add_edge( c, c + 1 );
        c++;
        xmin = p[XX];
    }
    // Left side
    {
        for ( j = fs.j_end() - 1; j >= fs.j_begin(); --j ) {
            p[YY] = grid.y( j );
            i     = fs.i_begin( j );
            if ( i == 0 ) {
                p[XX] = dom.xmin();
            }
            else {
                p[XX] = 0.5 * ( grid.x( i - 1, j ) + grid.x( i, j ) );
            }

            if ( j > fs.j_begin() ) {
                xmin = std::max( xmin, p[XX] );
            }
            if ( j == fs.j_begin() + 1 && xmin < points[0][XX] ) {
                idx_t jtop = fs.j_begin();
                idx_t itop = fs.i_begin( jtop );
                if ( xmin + grid.xspace().dx()[jtop] > grid.x( itop, jtop ) ) {
                    xmin = std::max( xmin, 0.5 * ( grid.x( itop - 1, jtop ) + grid.x( itop, jtop ) ) );
                }
                else {
                    ymax = 0.5 * ( p[YY] + grid.y( jtop ) );
                }
            }

            if ( j < fs.j_end() - 1 && not equal( p[XX], points.back()[XX] ) ) {
                if ( not debug ) {
                    make_corner();
                }
                else {
                    PointXY ptmp;
                    ptmp  = p;
                    p[XX] = points.back()[XX];
                    p[YY] = 0.5 * ( points.back()[YY] + grid.y( j ) );
                    points.emplace_back( p );
                    add_edge( c, c + 1 );
                    c++;
                    p = ptmp;


                    ptmp  = p;
                    p[YY] = points.back()[YY];
                    points.emplace_back( p );
                    add_edge( c, c + 1 );
                    c++;
                    p = ptmp;
                }
            }


            points.emplace_back( p );
            add_edge( c, c + 1 );
            c++;
        }
    }
    // Connect to top
    add_edge( c, 0 );


    bb = std::vector<Point2>{{xmin, ymax}, {xmin, ymin}, {xmax, ymin}, {xmax, ymax}};

    return edges;
}

StructuredPartitionPolygon::StructuredPartitionPolygon( const functionspace::FunctionSpaceImpl& fs, idx_t halo ) :
    fs_( fs ), halo_( halo ) {
    ATLAS_TRACE( "StructuredPartitionPolygon" );
    setup( compute_edges( fs, halo, points_, inner_bounding_box_ ) );
    points_.emplace_back( points_[0] );
    inner_bounding_box_.emplace_back( inner_bounding_box_[0] );
    auto min = Point2( std::numeric_limits<double>::max(), std::numeric_limits<double>::max() );
    auto max = Point2( std::numeric_limits<double>::lowest(), std::numeric_limits<double>::lowest() );
    for ( size_t i = 0; i < inner_bounding_box_.size(); ++i ) {
        min = Point2::componentsMin( min, inner_bounding_box_[i] );
        max = Point2::componentsMax( max, inner_bounding_box_[i] );
    }
    inscribed_domain_ = {{min[XX], max[XX]}, {min[YY], max[YY]}};
}

size_t StructuredPartitionPolygon::footprint() const {
    size_t size = sizeof( *this );
    size += capacity() * sizeof( idx_t );
    return size;
}

void StructuredPartitionPolygon::outputPythonScript( const eckit::PathName& filepath,
                                                     const eckit::Configuration& config ) const {
    ATLAS_TRACE( "Output PartitionPolygon" );
    const eckit::mpi::Comm& comm = atlas::mpi::comm();
    int mpi_rank                 = int( comm.rank() );
    int mpi_size                 = int( comm.size() );

    auto xy = array::make_view<double, 2>( dynamic_cast<const functionspace::detail::StructuredColumns&>( fs_ ).xy() );

    const std::vector<Point2>& points = config.getBool( "inner_bounding_box", false ) ? inner_bounding_box_ : points_;

    double xmin = std::numeric_limits<double>::max();
    double xmax = -std::numeric_limits<double>::max();
    for ( size_t i = 0; i < points.size(); ++i ) {
        xmin = std::min( xmin, points[i][XX] );
        xmax = std::max( xmax, points[i][XX] );
    }
    comm.allReduceInPlace( xmin, eckit::mpi::min() );
    comm.allReduceInPlace( xmax, eckit::mpi::max() );

    idx_t count     = dynamic_cast<const functionspace::detail::StructuredColumns&>( fs_ ).sizeOwned();
    idx_t count_all = count;
    comm.allReduceInPlace( count_all, eckit::mpi::sum() );

    bool plot_nodes = config.getBool( "nodes", false );
    for ( int r = 0; r < mpi_size; ++r ) {
        // clang-format off
        if ( mpi_rank == r ) {
            std::ofstream f( filepath.asString().c_str(), mpi_rank == 0 ? std::ios::trunc : std::ios::app );

            if ( mpi_rank == 0 ) {
                f << "\n" "import sys"
                     "\n"
                     "\n" "# Configuration option to plot nodes"
                     "\n" "plot_nodes = False"
                     "\n" "for argv in sys.argv[1:] :"
                     "\n" "  if argv == \"--nodes\" :"
                     "\n" "      plot_nodes = " + std::string( plot_nodes ? "True" : "False" ) +
                     "\n"
                     "\n" "import matplotlib.pyplot as plt"
                     "\n" "from matplotlib.path import Path"
                     "\n" "import matplotlib.patches as patches"
                     "\n" 
                     "\n" "from itertools import cycle"
                     "\n" "import matplotlib.cm as cm"
                     "\n" "import numpy as np"
                     "\n" "cycol = cycle([cm.Paired(i) for i in "
                          "np.linspace(0,1,12,endpoint=True)]).next"
                     "\n"
                     "\n" "fig = plt.figure()"
                     "\n" "ax = fig.add_subplot(111,aspect='equal')"
                     "\n";
            }
            f << "\n" "verts_" << r << " = [";
            for ( size_t i=0; i<points.size(); ++i ) {
                f << "\n  (" << points[i][XX] << ", " << points[i][YY] << "), ";
            }
            f << "\n]"
                 "\n"
                 "\n" "codes_" << r << " = [Path.MOVETO]"
                 "\n" "codes_" << r << ".extend([Path.LINETO] * " << ( points.size() - 2 ) << ")"
                 "\n" "codes_" << r << ".extend([Path.CLOSEPOLY])"
                 "\n"
                 "\n" "count_" << r << " = " << count <<
                 "\n" "count_all_" << r << " = " << count_all <<
                 "\n";
            if ( plot_nodes ) {
                f << "\n" "x_" << r << " = [";
                for ( idx_t i = 0; i < count; ++i ) {
                    f << xy( i, XX ) << ", ";
                }
                f << "]"
                     "\n" "y_" << r << " = [";
                for ( idx_t i = 0; i < count; ++i ) {
                    f << xy( i, YY ) << ", ";
                }
                f << "]";
            }
            f << "\n"
                 "\n" "c = cycol()"
                 "\n" "ax.add_patch(patches.PathPatch(Path(verts_" << r << ", codes_" << r << "), facecolor=c, color=c, alpha=0.3, lw=1))";
            if ( plot_nodes ) {
                f << "\n" "if plot_nodes:"
                     "\n" "    ax.scatter(x_" << r << ", y_" << r << ", color=c, marker='o')";
            }
            f << "\n";
            if ( mpi_rank == mpi_size - 1 ) {
                f << "\n" "ax.set_xlim( " << xmin << "-5, " << xmax << "+5)"
                     "\n" "ax.set_ylim(-90-5,  90+5)"
                     "\n" "ax.set_xticks([0,45,90,135,180,225,270,315,360])"
                     "\n" "ax.set_yticks([-90,-45,0,45,90])"
                     "\n" "plt.grid()"
                     "\n" "plt.show()";
            }
        }
        // clang-format on
        comm.barrier();
    }
}

// void StructuredPartitionPolygon::print( std::ostream& out ) const {
//     out << "polygon:{"
//         << "halo:" << halo_ << ",size:" << polygon_.size() << ",nodes:" << static_cast<const util::Polygon&>( polygon_ )
//         << "}";
// }

// ----------------------------------------------------------------------------

}  // namespace grid
}  // namespace atlas
