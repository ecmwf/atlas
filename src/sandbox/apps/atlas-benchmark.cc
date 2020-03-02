/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

/**
 * @file atlas-benchmark.cc
 * @author Willem Deconinck
 *
 * Benchmark testing parallel performance of gradient computation using the
 * Green-Gauss Theorem on an edge-based median-dual mesh.
 *
 * Configurable is
 *   - Horizontal mesh resolution, which is unstructured and
 *     domain-decomposed,
 *   - Vertical resolution, which is structured, and is beneficial for caching
 *   - Number of iterations, so caches can warm up, and timings can be averaged
 *   - Number of OpenMP threads per MPI task
 *
 * Results should be bit-identical when changing number of OpenMP threads or MPI
 * tasks.
 * A checksum on all bits is used to verify between scaling runs.
 *
 *
 */

#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <vector>

#include "atlas/functionspace.h"
#include "atlas/grid.h"
#include "atlas/library.h"
#include "atlas/mesh.h"
#include "atlas/mesh/IsGhostNode.h"
#include "atlas/mesh/actions/BuildDualMesh.h"
#include "atlas/mesh/actions/BuildEdges.h"
#include "atlas/mesh/actions/BuildHalo.h"
#include "atlas/mesh/actions/BuildParallelFields.h"
#include "atlas/mesh/actions/BuildPeriodicBoundaries.h"
#include "atlas/mesh/actions/Reorder.h"
#include "atlas/meshgenerator.h"
#include "atlas/output/Gmsh.h"
#include "atlas/parallel/Checksum.h"
#include "atlas/parallel/HaloExchange.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/parallel/omp/omp.h"
#include "atlas/runtime/AtlasTool.h"
#include "atlas/runtime/Exception.h"
#include "atlas/runtime/Trace.h"
#include "atlas/util/CoordinateEnums.h"
#include "atlas/util/Earth.h"

//----------------------------------------------------------------------------------------------------------------------

using std::cout;
using std::endl;
using std::fixed;
using std::max;
using std::min;
using std::numeric_limits;
using std::scientific;
using std::setprecision;
using std::setw;
using std::string;
using std::stringstream;
using std::unique_ptr;
using std::vector;

using Topology = atlas::mesh::Nodes::Topology;
using atlas::mesh::IsGhostNode;

using namespace eckit::option;
using namespace atlas;
using namespace atlas::array;
using namespace atlas::grid;
using namespace atlas::mesh::actions;
using namespace atlas::functionspace;
using namespace atlas::meshgenerator;
using atlas::AtlasTool;

//----------------------------------------------------------------------------------------------------------------------

struct TimerStats {
    TimerStats( const string& _name = "timer" ) {
        max  = -1;
        min  = -1;
        avg  = 0;
        cnt  = 0;
        name = _name;
    }
    void update( Trace& timer ) {
        double t = timer.elapsed();
        if ( min < 0 ) {
            min = t;
        }
        if ( max < 0 ) {
            max = t;
        }
        min = std::min( min, t );
        max = std::max( max, t );
        avg = ( avg * cnt + t ) / ( cnt + 1 );
        ++cnt;
    }
    string str() {
        stringstream stream;
        stream << name << ": min, max, avg -- " << fixed << setprecision( 5 ) << min << ", " << fixed
               << setprecision( 5 ) << max << ", " << fixed << setprecision( 5 ) << avg;
        return stream.str();
    }
    string name;
    double max;
    double min;
    double avg;
    int cnt;
};

//----------------------------------------------------------------------------------------------------------------------

class AtlasBenchmark : public AtlasTool {
    int execute( const Args& args ) override;

public:
    AtlasBenchmark( int argc, char** argv ) : AtlasTool( argc, argv ) {
        add_option( new SimpleOption<std::string>( "grid", "Grid unique identifier" ) );
        add_option( new SimpleOption<size_t>( "nlev", "Vertical resolution: Number of levels" ) );
        add_option( new SimpleOption<size_t>( "niter", "Number of iterations" ) );
        add_option( new SimpleOption<size_t>( "omp", "Number of OpenMP threads per MPI task" ) );
        add_option( new SimpleOption<bool>( "progress", "Show progress bar instead of intermediate timings" ) );
        add_option( new SimpleOption<bool>( "output", "Write output in gmsh format" ) );
        add_option( new SimpleOption<long>( "exclude", "Exclude number of iterations in statistics (default=1)" ) );
        add_option( new SimpleOption<bool>( "details", "Show detailed timers (default=false)" ) );
        add_option( new SimpleOption<std::string>( "reorder", "Reorder mesh (default=none)" ) );
        add_option( new SimpleOption<bool>( "sort_edges", "Sort edges by lowest node local index" ) );
    }

    void setup();

    void iteration();

    double result();

    int verify( const double& );

    void initial_condition( Field& field, const double& beta );

private:
    Mesh mesh;
    functionspace::NodeColumns nodes_fs;
    functionspace::NodeColumns edges_fs;
    Field scalar_field;
    Field grad_field;

    unique_ptr<array::Array> avgS_arr;

    std::vector<idx_t> pole_edges;
    std::vector<bool> is_ghost;

    idx_t nnodes;
    idx_t nedges;
    idx_t nlev;
    idx_t niter;
    idx_t exclude;
    bool output;
    long omp_threads;
    double dz;
    std::string gridname;
    std::string reorder{"none"};
    bool sort_edges{false};

    TimerStats iteration_timer;
    TimerStats haloexchange_timer;
    idx_t iter;
    bool progress;
};

//----------------------------------------------------------------------------------------------------------------------

int AtlasBenchmark::execute( const Args& args ) {
    Trace timer( Here(), "atlas-benchmark" );
    // Timer::Logging set_channel( Log::info() );

    nlev = 137;
    args.get( "nlev", nlev );
    gridname = "N64";
    args.get( "grid", gridname );
    niter = 100;
    args.get( "niter", niter );
    omp_threads = -1;
    args.get( "omp", omp_threads );
    progress = false;
    args.get( "progress", progress );
    exclude = niter == 1 ? 0 : 1;
    args.get( "exclude", exclude );
    output = false;
    args.get( "output", output );
    args.get( "reorder", reorder );
    args.get( "sort_edges", sort_edges );
    bool help( false );
    args.get( "help", help );

    iteration_timer    = TimerStats( "iteration" );
    haloexchange_timer = TimerStats( "halo-exchange" );

    if ( omp_threads > 0 ) {
        atlas_omp_set_num_threads( omp_threads );
    }

    Log::info() << "atlas-benchmark\n" << endl;
    Log::info() << Library::instance().information() << endl;
    Log::info() << "Configuration:" << endl;
    Log::info() << "  grid: " << gridname << endl;
    Log::info() << "  nlev: " << nlev << endl;
    Log::info() << "  niter: " << niter << endl;
    Log::info() << endl;
    Log::info() << "  MPI tasks: " << mpi::comm().size() << endl;
    Log::info() << "  OpenMP threads per MPI task: " << atlas_omp_get_max_threads() << endl;
    Log::info() << endl;

    Log::info() << "Timings:" << endl;

    ATLAS_TRACE_SCOPE( "setup", {"atlas-benchmark-setup"} ) { setup(); }

    Log::info() << "  Executing " << niter << " iterations: \n";
    if ( progress ) {
        Log::info() << "      0%   10   20   30   40   50   60   70   80   90   100%\n";
        Log::info() << "      |----|----|----|----|----|----|----|----|----|----|\n";
        Log::info() << "      " << std::flush;
    }
    unsigned int tic = 0;
    for ( iter = 0; iter < niter; ++iter ) {
        if ( progress ) {
            unsigned int tics_needed =
                static_cast<unsigned int>( static_cast<double>( iter ) / static_cast<double>( niter - 1 ) * 50.0 );
            while ( tic <= tics_needed ) {
                Log::info() << '*' << std::flush;
                ++tic;
            }
            if ( iter == niter - 1 ) {
                if ( tic < 51 ) {
                    Log::info() << '*';
                }
                Log::info() << endl;
            }
        }
        iteration();
    }
    timer.stop();

    Log::info() << "Iteration timer Statistics:\n"
                << "  min: " << setprecision( 5 ) << fixed << iteration_timer.min << "  max: " << setprecision( 5 )
                << fixed << iteration_timer.max << "  avg: " << setprecision( 5 ) << fixed << iteration_timer.avg
                << endl;
    Log::info() << "Communication timer Statistics:\n"
                << "  min: " << setprecision( 5 ) << fixed << haloexchange_timer.min << "  max: " << setprecision( 5 )
                << fixed << haloexchange_timer.max << "  avg: " << setprecision( 5 ) << fixed << haloexchange_timer.avg
                << " ( " << setprecision( 2 ) << haloexchange_timer.avg / iteration_timer.avg * 100. << "% )" << endl;

    util::Config report_config;
    report_config.set( "indent", 4 );
    if ( not args.getBool( "details", false ) ) {
        report_config.set( "exclude", std::vector<std::string>{"halo-exchange", "atlas-benchmark-setup/*"} );
    }
    Log::info() << timer.report( report_config ) << std::endl;
    Log::info() << endl;

    mpi::comm().barrier();

    Log::info() << "Results:" << endl;

    double res = result();

    Log::info() << endl;
    return verify( res );
}

//----------------------------------------------------------------------------------------------------------------------

void AtlasBenchmark::initial_condition( Field& field, const double& beta ) {
    const double radius  = util::Earth::radius();
    const double USCAL   = 20.;
    const double pvel    = USCAL / radius;
    const double deg2rad = M_PI / 180.;

    auto lonlat_deg = array::make_view<double, 2>( mesh.nodes().lonlat() );
    auto var        = array::make_view<double, 2>( field );

    idx_t nnodes = mesh.nodes().size();
    for ( idx_t jnode = 0; jnode < nnodes; ++jnode ) {
        double x = lonlat_deg( jnode, LON ) * deg2rad;
        double y = lonlat_deg( jnode, LAT ) * deg2rad;
        double Ux =
            pvel * ( std::cos( beta ) + std::tan( y ) * std::cos( x ) * std::sin( beta ) ) * radius * std::cos( y );
        double Uy = -pvel * std::sin( x ) * std::sin( beta ) * radius;
        for ( idx_t jlev = 0; jlev < nlev; ++jlev ) {
            var( jnode, jlev ) = std::sqrt( Ux * Ux + Uy * Uy );
        }
    }
}

//----------------------------------------------------------------------------------------------------------------------

void AtlasBenchmark::setup() {
    idx_t halo = 1;

    StructuredGrid grid;
    ATLAS_TRACE_SCOPE( "Create grid" ) { grid = Grid( gridname ); }
    ATLAS_TRACE_SCOPE( "Create mesh" ) {
        mesh = MeshGenerator( "structured", util::Config( "partitioner", "equal_regions" ) ).generate( grid );
    }
    mesh::actions::Reorder{option::type( reorder )}( mesh );

    ATLAS_TRACE_SCOPE( "Create node_fs" ) { nodes_fs = functionspace::NodeColumns( mesh, option::halo( halo ) ); }
    ATLAS_TRACE_SCOPE( "Create edges_fs" ) {
        edges_fs = functionspace::EdgeColumns( mesh, option::halo( halo ) | util::Config( "sort_edges", sort_edges ) );
    }

    // mesh.polygon(0).outputPythonScript("plot_polygon.py");
    //  atlas::output::Output gmsh = atlas::output::Gmsh( "edges.msh",
    //  util::Config("ghost",true)("edges",true)("elements",false) );
    //  gmsh.write( mesh );

    //  gmsh = atlas::output::Gmsh( "elements.msh",
    //  util::Config("ghost",true)("edges",false)("elements",true) );
    //  gmsh.write( mesh );

    ATLAS_TRACE_SCOPE( "build_median_dual_mesh" ) { build_median_dual_mesh( mesh ); }
    ATLAS_TRACE_SCOPE( "build_node_to_edge_connectivity" ) { build_node_to_edge_connectivity( mesh ); }

    scalar_field = nodes_fs.createField<double>( option::name( "field" ) | option::levels( nlev ) );
    grad_field =
        nodes_fs.createField<double>( option::name( "grad" ) | option::levels( nlev ) | option::variables( 3 ) );

    nnodes = mesh.nodes().size();
    nedges = mesh.edges().size();

    auto lonlat = array::make_view<double, 2>( mesh.nodes().xy() );
    auto V      = array::make_view<double, 1>( mesh.nodes().field( "dual_volumes" ) );
    auto S      = array::make_view<double, 2>( mesh.edges().field( "dual_normals" ) );

    initial_condition( scalar_field, 0. );

    double radius  = 6371.22e+03;  // Earth's radius
    double height  = 80.e+03;      // Height of atmosphere
    double deg2rad = M_PI / 180.;
    atlas_omp_parallel_for( idx_t jnode = 0; jnode < nnodes; ++jnode ) {
        lonlat( jnode, LON ) = lonlat( jnode, LON ) * deg2rad;
        lonlat( jnode, LAT ) = lonlat( jnode, LAT ) * deg2rad;
        double y             = lonlat( jnode, LAT );
        double hx            = radius * std::cos( y );
        double hy            = radius;
        double G             = hx * hy;
        V( jnode ) *= std::pow( deg2rad, 2 ) * G;
    }
    atlas_omp_parallel_for( idx_t jedge = 0; jedge < nedges; ++jedge ) {
        S( jedge, LON ) *= deg2rad;
        S( jedge, LAT ) *= deg2rad;
    }
    dz = height / static_cast<double>( nlev );

    const mesh::Connectivity& node2edge           = mesh.nodes().edge_connectivity();
    const mesh::MultiBlockConnectivity& edge2node = mesh.edges().node_connectivity();
    auto node2edge_sign                           = array::make_view<double, 2>( mesh.nodes().add(
        Field( "to_edge_sign", array::make_datatype<double>(), array::make_shape( nnodes, node2edge.maxcols() ) ) ) );

    atlas_omp_parallel_for( idx_t jnode = 0; jnode < nnodes; ++jnode ) {
        for ( idx_t jedge = 0; jedge < node2edge.cols( jnode ); ++jedge ) {
            idx_t iedge = node2edge( jnode, jedge );
            idx_t ip1   = edge2node( iedge, 0 );
            if ( jnode == ip1 ) {
                node2edge_sign( jnode, jedge ) = 1.;
            }
            else {
                node2edge_sign( jnode, jedge ) = -1.;
            }
        }
    }

    auto edge_flags   = array::make_view<int, 1>( mesh.edges().flags() );
    auto is_pole_edge = [&]( size_t e ) { return Topology::check( edge_flags( e ), Topology::POLE ); };

    std::vector<idx_t> tmp( nedges );
    int c( 0 );
    for ( idx_t jedge = 0; jedge < nedges; ++jedge ) {
        if ( is_pole_edge( jedge ) ) {
            tmp[c++] = jedge;
        }
    }
    pole_edges.reserve( c );
    for ( idx_t jedge = 0; jedge < c; ++jedge ) {
        pole_edges.push_back( tmp[jedge] );
    }

    auto flags = array::make_view<int, 1>( mesh.nodes().flags() );
    is_ghost.reserve( nnodes );
    for ( idx_t jnode = 0; jnode < nnodes; ++jnode ) {
        is_ghost.push_back( Topology::check( flags( jnode ), Topology::GHOST ) );
    }
}

//----------------------------------------------------------------------------------------------------------------------

void AtlasBenchmark::iteration() {
    Trace t( Here() );
    Trace compute( Here(), "compute" );
    if ( !avgS_arr ) {
        avgS_arr.reset( array::Array::create<double>( nedges, nlev, 2ul ) );
    }
    const auto& node2edge     = mesh.nodes().edge_connectivity();
    const auto& edge2node     = mesh.edges().node_connectivity();
    const auto field          = array::make_view<double, 2>( scalar_field );
    const auto S              = array::make_view<double, 2>( mesh.edges().field( "dual_normals" ) );
    const auto V              = array::make_view<double, 1>( mesh.nodes().field( "dual_volumes" ) );
    const auto node2edge_sign = array::make_view<double, 2>( mesh.nodes().field( "to_edge_sign" ) );

    auto grad = array::make_view<double, 3>( grad_field );
    auto avgS = array::make_view<double, 3>( *avgS_arr );

    atlas_omp_parallel_for( idx_t jedge = 0; jedge < nedges; ++jedge ) {
        idx_t ip1 = edge2node( jedge, 0 );
        idx_t ip2 = edge2node( jedge, 1 );

        for ( idx_t jlev = 0; jlev < nlev; ++jlev ) {
            double avg               = ( field( ip1, jlev ) + field( ip2, jlev ) ) * 0.5;
            avgS( jedge, jlev, LON ) = S( jedge, LON ) * avg;
            avgS( jedge, jlev, LAT ) = S( jedge, LAT ) * avg;
        }
    }

    atlas_omp_parallel_for( idx_t jnode = 0; jnode < nnodes; ++jnode ) {
        for ( idx_t jlev = 0; jlev < nlev; ++jlev ) {
            grad( jnode, jlev, LON ) = 0.;
            grad( jnode, jlev, LAT ) = 0.;
        }
        for ( idx_t jedge = 0; jedge < node2edge.cols( jnode ); ++jedge ) {
            idx_t iedge = node2edge( jnode, jedge );
            double add  = node2edge_sign( jnode, jedge );
            for ( idx_t jlev = 0; jlev < nlev; ++jlev ) {
                grad( jnode, jlev, LON ) += add * avgS( iedge, jlev, LON );
                grad( jnode, jlev, LAT ) += add * avgS( iedge, jlev, LAT );
            }
        }
        for ( idx_t jlev = 0; jlev < nlev; ++jlev ) {
            grad( jnode, jlev, LON ) /= V( jnode );
            grad( jnode, jlev, LAT ) /= V( jnode );
        }
    }
    // special treatment for the north & south pole cell faces
    // Sx == 0 at pole, and Sy has same sign at both sides of pole
    for ( size_t jedge = 0; jedge < pole_edges.size(); ++jedge ) {
        idx_t iedge = pole_edges[jedge];
        idx_t ip2   = edge2node( iedge, 1 );
        // correct for wrong Y-derivatives in previous loop
        for ( idx_t jlev = 0; jlev < nlev; ++jlev ) {
            grad( ip2, jlev, LAT ) += 2. * avgS( iedge, jlev, LAT ) / V( ip2 );
        }
    }

    double dzi   = 1. / dz;
    double dzi_2 = 0.5 * dzi;

    atlas_omp_parallel_for( idx_t jnode = 0; jnode < nnodes; ++jnode ) {
        if ( nlev > 2 ) {
            for ( idx_t jlev = 1; jlev < nlev - 1; ++jlev ) {
                grad( jnode, jlev, ZZ ) = ( field( jnode, jlev + 1 ) - field( jnode, jlev - 1 ) ) * dzi_2;
            }
        }
        if ( nlev > 1 ) {
            grad( jnode, 0, ZZ )        = ( field( jnode, 1 ) - field( jnode, 0 ) ) * dzi;
            grad( jnode, nlev - 1, ZZ ) = ( field( jnode, nlev - 2 ) - field( jnode, nlev - 1 ) ) * dzi;
        }
        if ( nlev == 1 ) {
            grad( jnode, 0, ZZ ) = 0.;
        }
    }
    compute.stop();

    // halo-exchange
    Trace halo( Here(), "halo-exchange" );
    nodes_fs.halo_exchange().execute<double, 3>( grad_field.array() );
    halo.stop();

    t.stop();

    if ( iter >= exclude ) {
        haloexchange_timer.update( halo );
        iteration_timer.update( t );
    }

    if ( !progress ) {
        Log::info() << setw( 6 ) << iter + 1 << "    total: " << fixed << setprecision( 5 ) << t.elapsed()
                    << "    communication: " << setprecision( 5 ) << halo.elapsed() << " ( " << setprecision( 2 )
                    << fixed << setw( 3 ) << halo.elapsed() / t.elapsed() * 100 << "% )" << endl;
    }
}

//----------------------------------------------------------------------------------------------------------------------

template <typename DATA_TYPE>
DATA_TYPE vecnorm( const DATA_TYPE vec[], size_t size ) {
    DATA_TYPE norm = 0;
    for ( size_t j = 0; j < size; ++j ) {
        norm += vec[j] * vec[j];
    }
    return norm;
}

double AtlasBenchmark::result() {
    auto grad     = array::make_view<double, 3>( grad_field );
    double maxval = -std::numeric_limits<double>::max();
    double minval = std::numeric_limits<double>::max();
    ;
    double norm = 0.;

    nodes_fs.haloExchange( grad_field );
    for ( idx_t jnode = 0; jnode < nnodes; ++jnode ) {
        if ( !is_ghost[jnode] ) {
            for ( idx_t jlev = 0; jlev < 1; ++jlev ) {
                const double scaling = 1.e12;
                grad( jnode, jlev, LON ) *= scaling;
                grad( jnode, jlev, LAT ) *= scaling;
                grad( jnode, jlev, ZZ ) *= scaling;

                std::array<double, 3> v;
                v[0] = grad( jnode, jlev, LON );
                v[1] = grad( jnode, jlev, LAT );
                v[2] = grad( jnode, jlev, ZZ );

                maxval = std::max( maxval, v[0] );
                maxval = std::max( maxval, v[1] );
                maxval = std::max( maxval, v[2] );
                minval = std::min( minval, v[0] );
                minval = std::min( minval, v[1] );
                minval = std::min( minval, v[2] );

                // if( mpi::comm().rank() == 478 ) {
                //  std::cout << "    " << jnode << "   part " << part(jnode) << "
                //  glb_idx " <<  glb_idx(jnode) << "   x,y,z " << v[0] << "," << v[1]
                //  << ","<< v[2] <<  std::endl;
                //}

                norm += v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
            }
        }
    }

    if ( output ) {
        std::vector<long> levels( 1, 0 );
        atlas::output::Output gmsh =
            atlas::output::Gmsh( "benchmark.msh", util::Config( "levels", levels )( "ghost", true ) );
        gmsh.write( mesh );
        gmsh.write( scalar_field );
        gmsh.write( grad_field );
    }

    ATLAS_TRACE_MPI( ALLREDUCE ) {
        mpi::comm().allReduceInPlace( maxval, eckit::mpi::max() );
        mpi::comm().allReduceInPlace( minval, eckit::mpi::min() );
        mpi::comm().allReduceInPlace( norm, eckit::mpi::sum() );
    }

    norm = std::sqrt( norm );

    Log::info() << "  maxval: " << setw( 13 ) << setprecision( 6 ) << scientific << maxval << endl;
    Log::info() << "  minval: " << setw( 13 ) << setprecision( 6 ) << scientific << minval << endl;
    Log::info() << "  norm:   " << setw( 13 ) << setprecision( 6 ) << scientific << norm << endl;

    Log::info() << "  checksum: " << nodes_fs.checksum().execute( grad ) << endl;

    return norm;
}

int AtlasBenchmark::verify( const double& norm ) {
    Log::warning() << "Verification is not yet implemented" << endl;
    return success();
}

//----------------------------------------------------------------------------------------------------------------------

int main( int argc, char** argv ) {
    AtlasBenchmark tool( argc, argv );
    return tool.start();
}
