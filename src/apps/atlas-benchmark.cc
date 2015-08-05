/*
 * (C) Copyright 1996-2014 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
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
 * Results should be bit-identical when changing number of OpenMP threads or MPI tasks.
 * A checksum on all bits is used to verify between scaling runs.
 *
 *
 */

#include <limits>
#include <cassert>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <vector>
#include <memory>

#include "eckit/exception/Exceptions.h"
#include "eckit/config/Resource.h"
#include "eckit/runtime/Tool.h"
#include "eckit/runtime/Context.h"
#include "eckit/filesystem/PathName.h"
#include "eckit/memory/Factory.h"
#include "eckit/memory/Builder.h"
#include "eckit/parser/JSON.h"
#include "eckit/log/Timer.h"


#include "atlas/atlas.h"
#include "atlas/util/Array.h"
#include "atlas/io/Gmsh.h"
#include "atlas/actions/GenerateMesh.h"
#include "atlas/actions/BuildEdges.h"
#include "atlas/actions/BuildPeriodicBoundaries.h"
#include "atlas/actions/BuildHalo.h"
#include "atlas/actions/BuildParallelFields.h"
#include "atlas/actions/BuildDualMesh.h"
#include "atlas/Mesh.h"
#include "atlas/Nodes.h"
#include "atlas/grids/grids.h"
#include "atlas/io/Gmsh.h"
#include "atlas/atlas_omp.h"
#include "atlas/mpi/mpi.h"
#include "atlas/util/Bitflags.h"
#include "atlas/util/IndexView.h"
#include "atlas/util/IsGhost.h"

//------------------------------------------------------------------------------------------------------

using std::string;
using std::stringstream;
using std::min;
using std::max;
using std::vector;
using std::setw;
using std::setprecision;
using std::scientific;
using std::fixed;
using std::cout;
using std::endl;
using std::numeric_limits;

using atlas::util::Topology;
using atlas::util::IsGhost;

using namespace eckit;
using namespace atlas;
using namespace atlas::grids;
using namespace atlas::actions;


//------------------------------------------------------------------------------------------------------

struct TimerStats
{
  TimerStats(const string& _name = "timer")
  {
    max = -1;
    min = -1;
    avg = 0;
    cnt = 0;
    name = _name;
  }
  void update(Timer& timer)
  {
    double t = timer.elapsed();
    if( min < 0 ) min = t;
    if( max < 0 ) max = t;
    min = std::min(min, t);
    max = std::max(max, t);
    avg = (avg*cnt+t)/(cnt+1);
    ++cnt;
  }
  string str()
  {
    stringstream stream;
    stream << name << ": min, max, avg -- " << min << ", " << max << ", " << avg;
    return stream.str();
  }
  string name;
  double max;
  double min;
  double avg;
  int cnt;
};

class AtlasBenchmark : public eckit::Tool {

  virtual void run();

public:

  AtlasBenchmark(int argc,char **argv): eckit::Tool(argc,argv), do_run(true)
  {
    N     = Resource<size_t>("-N",1280);
    nlev  = Resource<size_t>("-nlev",137);
    niter = Resource<size_t>("-niter",100);
    omp_threads = Resource<long>("-omp",-1);
    progress = Resource<bool>("-progress",false);
    iteration_timer = TimerStats("iteration");
    haloexchange_timer = TimerStats("halo-exchange");
    exclude = Resource<int>("-exclude", niter==1?0:1);
    output = Resource<bool>("-output", false);
    bool help = Resource<bool>("-h",false);
    if( help )
    {
      string help_str =
          "NAME\n"
          "       atlas-benchmark - Benchmark parallel performance\n"
          "\n"
          "SYNOPSIS\n"
          "       atlas-benchmark [OPTIONS]...\n"
          "\n"
          "DESCRIPTION\n"
          "       Parallel performance of gradient computation using\n"
          "       Green-Gauss theorem on median-dual mesh based on\n"
          "       IFS reduced Gaussian grid\n"
          "\n"
          "       -N         Horizontal resolution: Number of\n"
          "                  latitudes between pole and equator\n"
          "\n"
          "       -nlev      Vertical resolution: Number of levels\n"
          "\n"
          "       -niter     Number of iterations to run\n"
          "\n"
          "       -omp       Number of threads per MPI task\n"
          "\n"
          "       -exclude   Exclude number of iterations in statistics (default=1)\n"
          "\n"
          "       -progress  Show progress bar instead of intermediate timings\n"
          "\n"
          "       -output    Write output in gmsh format\n"
          "\n"
          "AUTHOR\n"
          "       Written by Willem Deconinck.\n"
          "\n"
          "ECMWF                        December 2014"
          ;

      eckit::mpi::init();
      if( eckit::mpi::rank()==0 )
      {
        std::cout << help_str << endl;
      }
      eckit::mpi::finalize();
      do_run = false;
    }
  }

  void setup();

  void iteration();

  double result();

  int verify(const double&);

private:

  Mesh::Ptr mesh;
  IndexView<int,2> edge2node;


  ArrayView<double,2> lonlat;
  ArrayView<double,1> V;
  ArrayView<double,2> S;

  ArrayView<double,2> field;
  ArrayView<double,3> grad;

  IndexView<int,   2> node2edge;
  ArrayView<int,   1> node2edge_size;
  ArrayView<double,2> node2edge_sign;
  ArrayView<int,   1> edge_is_pole;
  vector<int> pole_edges;
  vector<bool> is_ghost;

  size_t nnodes;
  size_t nedges;
  size_t nlev;
  size_t N;
  size_t niter;
  long exclude;
  bool output;
  long omp_threads;
  double dz;

  TimerStats iteration_timer;
  TimerStats haloexchange_timer;
  size_t iter;
  bool progress;
  bool do_run;

public:
  int exit_code;

};

//------------------------------------------------------------------------------------------------------

void AtlasBenchmark::run()
{
  if( !do_run )
    return;

  atlas_init();

  if( omp_threads > 0 )
    omp_set_num_threads(omp_threads);

  Log::info() << "atlas-benchmark\n" << endl;
  Log::info() << "Atlas:" << endl;
  Log::info() << "  version:  ["<< atlas_version() << "]" << endl;
  Log::info() << "  git:      ["<< atlas_git_sha1() << "]" << endl;
  Log::info() << endl;
  Log::info() << "Configuration:" << endl;
  Log::info() << "  N: " << N << endl;
  Log::info() << "  nlev: " << nlev << endl;
  Log::info() << "  niter: " << niter << endl;
  Log::info() << endl;
  Log::info() << "  MPI tasks: "<<eckit::mpi::size()<<endl;
  Log::info() << "  OpenMP threads per MPI task: " << omp_get_max_threads() << endl;
  Log::info() << endl;

  Log::info() << "Timings:" << endl;

  setup();

  Log::info() << "  Executing " << niter << " iterations: \n";
  if( progress )
  {
    Log::info() << "      0%   10   20   30   40   50   60   70   80   90   100%\n";
    Log::info() << "      |----|----|----|----|----|----|----|----|----|----|\n";
    Log::info() << "      " << std::flush;
  }
  int tic=0;
  for( iter=0; iter<niter; ++iter )
  {
    if( progress )
    {
      unsigned int tics_needed = static_cast<unsigned int>(static_cast<double>(iter)/static_cast<double>(niter-1)*50.0);
      while( tic <= tics_needed )
      {
        Log::info() << '*' << std::flush;
        ++tic;
      }
      if ( iter == niter-1 )
      {
        if ( tic < 51 ) Log::info() << '*';
          Log::info() << endl;
      }
    }
    iteration();
  }


  Log::info() << "Iteration timer Statistics:\n"
              << "  min: " << setprecision(5) << fixed << iteration_timer.min
              << "  max: " << setprecision(5) << fixed << iteration_timer.max
              << "  avg: " << setprecision(5) << fixed << iteration_timer.avg << endl;
  Log::info() << "Communication timer Statistics:\n"
              << "  min: " << setprecision(5) << fixed << haloexchange_timer.min
              << "  max: " << setprecision(5) << fixed << haloexchange_timer.max
              << "  avg: " << setprecision(5) << fixed << haloexchange_timer.avg
              << " ( "<< setprecision(2) << haloexchange_timer.avg/iteration_timer.avg*100. << "% )" << endl;

  Log::info() << endl;
  Log::info() << "Results:" << endl;

  double res = result();

  Log::info() << endl;
  exit_code = verify( res );

  atlas_finalize();
}

//------------------------------------------------------------------------------------------------------

void AtlasBenchmark::setup()
{
  Timer timer( "setup", Log::debug());

  grids::load();

  stringstream gridname; gridname << "rgg.N"<<N;
  ReducedGrid::Ptr grid( ReducedGrid::create(gridname.str()) );
  mesh = Mesh::Ptr( generate_mesh(*grid) );
  // mesh = Mesh::Ptr( generate_regular_grid( 2*N, N) );
  build_nodes_parallel_fields(mesh->nodes());
  build_periodic_boundaries(*mesh);
  build_halo(*mesh,1);
  renumber_nodes_glb_idx(mesh->nodes());
  build_edges(*mesh);
  build_pole_edges(*mesh);
  build_edges_parallel_fields(mesh->function_space("edges"),mesh->nodes());
  build_median_dual_mesh(*mesh);
  build_node_to_edge_connectivity(*mesh);

  ArrayView<double,2> lonlat ( mesh->nodes().lonlat() );

  nnodes = mesh->nodes().shape(0);
  nedges = mesh->function_space("edges").shape(0);

  edge2node  = IndexView<int,   2> ( mesh->function_space("edges").field("nodes") );
  lonlat = ArrayView<double,2> ( mesh->nodes().create_field<double>("lonlat",2) );
  V      = ArrayView<double,1> ( mesh->nodes().field("dual_volumes") );
  S      = ArrayView<double,2> ( mesh->function_space("edges").field("dual_normals") );
  field  = ArrayView<double,2> ( mesh->nodes().create_field<double>("field",nlev)  );
  Field& gradfield = ( mesh->nodes().create_field<double>("grad",nlev*3) );
  grad   = ArrayView<double,3> ( gradfield.data<double>(), make_shape(nnodes,nlev,3) );
  mesh->nodes().field("field").metadata().set("nb_levels",nlev);
  mesh->nodes().field("grad").metadata().set("nb_levels",nlev);

  double radius = 6371.22e+03; // Earth's radius
  double height = 80.e+03;     // Height of atmosphere
  double deg2rad = M_PI/180.;
  for( int jnode=0; jnode<nnodes; ++jnode )
  {
    lonlat(jnode,LON) = lonlat(jnode,LON) * deg2rad;
    lonlat(jnode,LAT) = lonlat(jnode,LAT) * deg2rad;
    double y  = lonlat(jnode,LAT);
    double hx = radius*std::cos(y);
    double hy = radius;
    double G  = hx*hy;
    V(jnode) *= std::pow(deg2rad,2) * G;

    for( int jlev=0; jlev<nlev; ++jlev )
      field(jnode,jlev) = 100.+50.*std::cos(2*y);
  }
  for( int jedge=0; jedge<nedges; ++jedge )
  {
    S(jedge,LON) *= deg2rad;
    S(jedge,LAT) *= deg2rad;
  }
  dz = height/static_cast<double>(nlev);

  edge_is_pole   = ArrayView<int,1> ( mesh->function_space("edges").field("is_pole_edge") );
  node2edge      = IndexView<int,2> ( mesh->nodes().field("to_edge") );
  node2edge_size = ArrayView<int,1> ( mesh->nodes().field("to_edge_size") );

  node2edge_sign = ArrayView<double,2> ( mesh->nodes().
                                         create_field<double>("to_edge_sign",node2edge.shape(1)) );

  for( int jnode=0; jnode<nnodes; ++jnode )
  {
    for( int jedge=0; jedge<node2edge_size(jnode); ++jedge )
    {
      int iedge = node2edge(jnode,jedge);
      int ip1 = edge2node(iedge,0);
      int ip2 = edge2node(iedge,1);
      if( jnode == ip1 )
        node2edge_sign(jnode,jedge) = 1.;
      else
        node2edge_sign(jnode,jedge) = -1.;
    }
  }

  vector<int> tmp(nedges);
  int c(0);
  for( int jedge=0; jedge<nedges; ++jedge )
  {
    if( edge_is_pole(jedge) )
      tmp[c++] = jedge;
  }
  pole_edges.reserve(c);
  for( int jedge=0; jedge<c; ++jedge )
    pole_edges.push_back(tmp[jedge]);

  ArrayView<int,1> flags( mesh->nodes().field("flags") );
  is_ghost.reserve(nnodes);
  for( int jnode=0; jnode<nnodes; ++jnode )
  {
    is_ghost.push_back( Topology::check(flags(jnode),Topology::GHOST) );
  }


  Log::info() << "  setup: " << timer.elapsed() << endl;


  // Check bit-reproducibility after setup()
  // ---------------------------------------
  //ArrayView<double,1> V ( mesh->nodes().field("dual_volumes") );
  //ArrayView<double,2> S ( mesh->function_space("edges").field("dual_normals") );
  //Log::info() << "  checksum coordinates : " << mesh->nodes().checksum().execute( lonlat ) << endl;
  //Log::info() << "  checksum dual_volumes: " << mesh->nodes().checksum().execute( V ) << endl;
  //Log::info() << "  checksum dual_normals: " << mesh->function_space("edges").checksum().execute( S ) << endl;
  //Log::info() << "  checksum field       : " << mesh->nodes().checksum().execute( field ) << endl;
}

//------------------------------------------------------------------------------------------------------

void AtlasBenchmark::iteration()
{
  Timer t("iteration", Log::debug(5));

  Array<double> avgS_arr(nedges,nlev,2);
  ArrayView<double,3> avgS(avgS_arr);

#ifdef ATLAS_HAVE_OMP
  #pragma omp parallel for
#endif
  for( int jedge=0; jedge<nedges; ++jedge )
  {
    int ip1 = edge2node(jedge,0);
    int ip2 = edge2node(jedge,1);

    for( int jlev=0; jlev<nlev; ++jlev )
    {
      double avg = ( field(ip1,jlev) + field(ip2,jlev) ) * 0.5;
      avgS(jedge,jlev,LON) = S(jedge,LON)*avg;
      avgS(jedge,jlev,LAT) = S(jedge,LAT)*avg;
    }
  }

#ifdef ATLAS_HAVE_OMP
  #pragma omp parallel for
#endif
  for( int jnode=0; jnode<nnodes; ++jnode )
  {
    for( int jlev=0; jlev<nlev; ++jlev )
    {
      grad(jnode,jlev,LON) = 0.;
      grad(jnode,jlev,LAT) = 0.;
    }
    for( int jedge=0; jedge<node2edge_size(jnode); ++jedge )
    {
      int iedge = node2edge(jnode,jedge);
      double add = node2edge_sign(jnode,jedge);
      for( int jlev=0; jlev<nlev; ++jlev )
      {
        grad(jnode,jlev,LON) += add*avgS(iedge,jlev,LON);
        grad(jnode,jlev,LAT) += add*avgS(iedge,jlev,LAT);
      }
    }
    for( int jlev=0; jlev<nlev; ++jlev )
    {
      grad(jnode,jlev,LON) /= V(jnode);
      grad(jnode,jlev,LAT) /= V(jnode);
    }
  }
  // special treatment for the north & south pole cell faces
  // Sx == 0 at pole, and Sy has same sign at both sides of pole
  for( int jedge=0; jedge<pole_edges.size();++jedge )
  {
    int iedge = pole_edges[jedge];
    int ip1 = edge2node(iedge,0);
    int ip2 = edge2node(iedge,1);
    // correct for wrong Y-derivatives in previous loop
    for( int jlev=0; jlev<nlev; ++jlev )
      grad(ip2,jlev,LAT) += 2.*avgS(iedge,jlev,LAT)/V(ip2);
  }

  double dzi = 1./dz;
  double dzi_2 = 0.5*dzi;

#ifdef ATLAS_HAVE_OMP
  #pragma omp parallel for
#endif
  for( int jnode=0; jnode<nnodes; ++jnode )
  {
    if( nlev > 2 )
    {
      for( int jlev=1; jlev<nlev-1; ++jlev )
      {
        grad(jnode,jlev,ZZ)   = (field(jnode,jlev+1)     - field(jnode,jlev-1))*dzi_2;
      }
    }
    if( nlev > 1 )
    {
      grad(jnode,  0   ,ZZ) = (field(jnode,  1   ) - field(jnode,  0   ))*dzi;
      grad(jnode,nlev-1,ZZ) = (field(jnode,nlev-2) - field(jnode,nlev-1))*dzi;
    }
    if( nlev == 1 )
      grad(jnode,0,ZZ) = 0.;
  }

  // halo-exchange
  eckit::mpi::barrier();
  Timer halo("halo-exchange", Log::debug(5));
  mesh->nodes().halo_exchange().execute(grad);
  eckit::mpi::barrier();
  t.stop();
  halo.stop();

  if( iter >= exclude )
  {
    haloexchange_timer.update(halo);
    iteration_timer.update(t);
  }

  if( !progress )
  {
    Log::info() << setw(6) << iter+1
                << "    total: " << fixed << setprecision(5) << t.elapsed()
                << "    communication: " << setprecision(5) << halo.elapsed()
                << " ( "<< setprecision(2) << fixed << setw(3)
                << halo.elapsed()/t.elapsed()*100 << "% )" << endl;
  }
}

//------------------------------------------------------------------------------------------------------

template< typename DATA_TYPE >
DATA_TYPE vecnorm( DATA_TYPE vec[], size_t size )
{
  DATA_TYPE norm=0;
  for( int j=0; j<size; ++j )
    norm += std::pow(vec[j],2);
  return std::sqrt(norm);
}

double AtlasBenchmark::result()
{
  double maxval = numeric_limits<double>::min();
  double minval = numeric_limits<double>::max();;
  double norm = 0.;
  for( int jnode=0; jnode<nnodes; ++jnode )
  {
    if( !is_ghost[jnode] )
    {
      for( int jlev=0; jlev<nlev; ++jlev )
      {
        maxval = max(maxval,grad(jnode,jlev,LON));
        maxval = max(maxval,grad(jnode,jlev,LAT));
        maxval = max(maxval,grad(jnode,jlev,ZZ));
        minval = min(minval,grad(jnode,jlev,LON));
        minval = min(minval,grad(jnode,jlev,LAT));
        minval = min(minval,grad(jnode,jlev,ZZ));
        norm += std::pow(vecnorm(grad[jnode][jlev].data(),3),2);
      }
    }
  }

  ECKIT_MPI_CHECK_RESULT( MPI_Allreduce(MPI_IN_PLACE,&maxval,1,eckit::mpi::datatype<double>(),MPI_MAX,eckit::mpi::comm()) );
  ECKIT_MPI_CHECK_RESULT( MPI_Allreduce(MPI_IN_PLACE,&minval,1,eckit::mpi::datatype<double>(),MPI_MIN,eckit::mpi::comm()) );
  ECKIT_MPI_CHECK_RESULT( MPI_Allreduce(MPI_IN_PLACE,&norm  ,1,eckit::mpi::datatype<double>(),MPI_SUM,eckit::mpi::comm()) );
  norm = std::sqrt(norm);

  Log::info() << "  checksum: " << mesh->nodes().checksum().execute( grad ) << endl;
  Log::info() << "  maxval: " << setw(13) << setprecision(6) << scientific << maxval << endl;
  Log::info() << "  minval: " << setw(13) << setprecision(6) << scientific << minval << endl;
  Log::info() << "  norm:   " << setw(13) << setprecision(6) << scientific << norm << endl;

  if( output )
  {
    //io::Gmsh().write(mesh->nodes().field("dual_volumes"),"benchmark.gmsh",std::ios_base::app);
    //io::Gmsh().write(mesh->nodes().field("field"),"benchmark.gmsh",std::ios_base::out);
    io::Gmsh().write(mesh->nodes().field("grad"),"benchmark.gmsh",std::ios_base::app);
  }
  return norm;
}

int AtlasBenchmark::verify(const double& norm)
{
  if( nlev != 137 )
  {
    Log::warning() << "Results cannot be verified with nlev != 137" << endl;
    return 1;
  }
  std::map<int,double> norms;
  norms[  16] = 1.473937e-09;
  norms[  24] = 2.090045e-09;
  norms[  32] = 2.736576e-09;
  norms[  48] = 3.980306e-09;
  norms[  64] = 5.219642e-09;
  norms[  80] = 6.451913e-09;
  norms[  96] = 7.647690e-09;
  norms[ 128] = 1.009042e-08;
  norms[ 160] = 1.254571e-08;
  norms[ 200] = 1.557589e-08;
  norms[ 256] = 1.983944e-08;
  norms[ 320] = 2.469347e-08;
  norms[ 400] = 3.076775e-08;
  norms[ 512] = 3.924470e-08;
  norms[ 576] = 4.409003e-08;
  norms[ 640] = 4.894316e-08;
  norms[ 800] = 6.104009e-08;
  norms[1024] = 7.796900e-08;
  norms[1280] = 9.733947e-08;
  norms[1600] = 1.215222e-07;
  norms[2000] = 1.517164e-07;
  norms[4000] = 2.939562e-07;

  if( norms.count(N) == 0 )
  {
    Log::warning() << "Results cannot be verified with resolution N="<< N << endl;
    return 1;
  }

  double diff = (norm-norms[N])/norms[N];
  if( diff < 0.01 )
  {
    Log::info() << "Results are verified and correct.\n  difference = " << setprecision(6) << fixed << diff*100 << "%" << endl;
    return 0;
  }

  Log::info() << "Results are wrong.\n  difference = " << setprecision(6) << fixed << diff*100 << "%" << endl;
  return 1;
}

//------------------------------------------------------------------------------------------------------

int main( int argc, char **argv )
{
  AtlasBenchmark tool(argc,argv);
  tool.start();
  return tool.exit_code;
}
