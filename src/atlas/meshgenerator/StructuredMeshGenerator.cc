/*
 * (C) Copyright 1996-2017 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include <numeric>
#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

#include "eckit/runtime/Main.h"

#include "atlas/library/config.h"
#include "atlas/grid/Partitioner.h"
#include "atlas/grid/Grid.h"
#include "atlas/grid/Distribution.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/mesh/Nodes.h"
#include "atlas/mesh/HybridElements.h"
#include "atlas/mesh/ElementType.h"
#include "atlas/mesh/Elements.h"
#include "atlas/meshgenerator/StructuredMeshGenerator.h"
#include "atlas/field/Field.h"
#include "atlas/util/CoordinateEnums.h"
#include "atlas/runtime/Log.h"
#include "atlas/array/Array.h"
#include "atlas/array/ArrayView.h"
#include "atlas/array/IndexView.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/runtime/Debug.h"

#define DEBUG_OUTPUT 0

using namespace eckit;
using atlas::mesh::Mesh;
using Topology = atlas::mesh::Nodes::Topology;

namespace atlas {
namespace meshgenerator {

namespace {
static double to_rad = M_PI/180.;
static double to_deg = 180.*M_1_PI;
}

struct Region
{
  int north;
  int south;
  array::ArrayT<int> elems;
  int ntriags;
  int nquads;
  int nnodes;
  std::vector<int> lat_begin;
  std::vector<int> lat_end;
  std::vector<int> nb_lat_elems;
};

StructuredMeshGenerator::StructuredMeshGenerator(const eckit::Parametrisation& p)
{
  configure_defaults();

  bool include_pole;
  if( p.get("include_pole",include_pole) )
    options.set("include_pole",include_pole);

  bool patch_pole;
  if( p.get("patch_pole",patch_pole) )
    options.set("patch_pole",patch_pole);

  bool unique_pole;
  if( p.get("unique_pole",unique_pole) )
    options.set("unique_pole",unique_pole);

  bool three_dimensional;
  if( p.get("three_dimensional",three_dimensional) || p.get("3d",three_dimensional) )
    options.set("3d",three_dimensional);

  size_t nb_parts;
  if( p.get("nb_parts",nb_parts) )
    options.set("nb_parts",nb_parts);

  size_t part;
  if( p.get("part",part) )
    options.set("part",part);

  double angle;
  if( p.get("angle",angle) )
    options.set("angle",angle);

  bool triangulate;
  if( p.get("triangulate",triangulate) )
    options.set("triangulate",triangulate);

  bool ghost_at_end;
  if( p.get("ghost_at_end",ghost_at_end) )
    options.set("ghost_at_end",ghost_at_end);

  std::string partitioner;
  if( grid::Partitioner::exists("trans") )
    partitioner = "trans";
  else
    partitioner = "equal_regions";
  options.set("partitioner",partitioner);

  if( p.get("partitioner",partitioner) )
  {
    if( not grid::Partitioner::exists(partitioner) ) {
      Log::warning() << "Atlas does not have support for partitioner " << partitioner << ". "
                     << "Defaulting to use partitioner EqualRegions" << std::endl;
      partitioner = "equal_regions";
    }
    options.set("partitioner",partitioner);
  }
}


void StructuredMeshGenerator::configure_defaults()
{
  // This option creates a point at the pole when true
  options.set( "include_pole", false );

  // This option sets the part that will be generated
  options.set( "patch_pole", false );

  // This option disregards multiple poles in grid (e.g. lonlat up to poles) and connects elements
  // to the first node only. Note this option will only be looked at in case other option
  // "3d"==true
  options.set( "unique_pole", true );

  // This option creates elements that connect east to west at greenwich meridian
  // when true, instead of creating periodic ghost-points at east boundary when false
  options.set( "3d", false );

  // This option sets number of parts the mesh will be split in
  options.set( "nb_parts", parallel::mpi::comm().size() );

  // This option sets the part that will be generated
  options.set( "part", parallel::mpi::comm().rank() );

  // Experimental option. The result is a non-standard Reduced Gaussian Grid, with a ragged Greenwich line
  options.set("stagger", false );

  // This option sets the maximum angle deviation for a quadrilateral element
  // angle = 30  -->  minimises number of triangles
  // angle = 0   -->  maximises number of triangles
  options.set<double>("angle", 0. );

  options.set<bool>("triangulate", false );

  options.set<bool>("ghost_at_end", true );

}

void StructuredMeshGenerator::generate(const grid::Grid& grid, Mesh& mesh ) const
{
    ASSERT(!mesh.generated());

  const grid::StructuredGrid rg = grid::StructuredGrid(grid);
  if( !rg )
    throw eckit::BadCast("Structured can only work with a Structured",Here());

  size_t nb_parts = options.get<size_t>("nb_parts");

  std::string partitioner_type = "trans";
  options.get("partitioner",partitioner_type);

  if ( rg.ny()%2 == 1 ) partitioner_type = "equal_regions"; // Odd number of latitudes
  if ( nb_parts == 1 || parallel::mpi::comm().size() == 1 ) partitioner_type = "equal_regions"; // Only one part --> Trans is slower

  grid::Partitioner partitioner( partitioner_type, grid, nb_parts );
  grid::Distribution distribution( partitioner.distribution() );
  generate( grid, distribution, mesh );
}

void StructuredMeshGenerator::hash(MD5& md5) const
{
    md5.add("Structured");
    options.hash(md5);
}

void StructuredMeshGenerator::generate(const grid::Grid& grid, const grid::Distribution& distribution, Mesh& mesh ) const
{
  const grid::StructuredGrid rg = grid::StructuredGrid(grid);
  if( !rg )
    throw eckit::BadCast("Grid could not be cast to a Structured",Here());

  ASSERT(!mesh.generated());

  if( grid.npts() != distribution.partition().size() )
  {
    std::stringstream msg;
    msg << "Number of points in grid ("<<grid.npts()<<") different from "
           "number of points in grid distribution ("<<distribution.partition().size()<<")";
    throw eckit::AssertionFailed(msg.str(),Here());
  }

  int mypart   = options.get<size_t>("part");


  // show distribution
#if DEBUG_OUTPUT
    int inode=0;
    std::vector<int> parts=distribution;
    Log::info(Here()) << "Partition : " << std::endl;
    for (int ilat=0; ilat<rg.ny(); ilat++) {
      for (int ilon=0; ilon<rg.nx(ilat); ilon++ ) {
        Log::info(Here()) << std::setw(3) << parts[inode];
        inode++;
      }
      Log::info(Here()) << "\n";
    }
#endif

  // clone some grid properties
  mesh.setProjection(rg.projection());

  Region region;
  generate_region(rg,distribution,mypart,region);

  generate_mesh(rg,distribution,region,mesh);
}


void StructuredMeshGenerator::generate_region(const grid::StructuredGrid& rg, const std::vector<int>& parts, int mypart, Region& region) const
{
  double max_angle          = options.get<double>("angle");
  bool   triangulate_quads  = options.get<bool>("triangulate");
  bool   three_dimensional  = options.get<bool>("3d");
  bool   has_north_pole = rg.y().front() == 90;
  bool   has_south_pole = rg.y().back()  == -90;
  bool   unique_pole        = options.get<bool>("unique_pole") && three_dimensional && has_north_pole && has_south_pole;
  bool   periodic_east_west = rg.periodic();


  int n;
  /*
  Find min and max latitudes used by this part.
  */
  n=0;
  int lat_north=-1;
  for(size_t jlat = 0; jlat < rg.ny(); ++jlat) {
    for(size_t jlon = 0; jlon < rg.nx(jlat); ++jlon) {
      if( parts.at(n) == mypart ) {
        lat_north=jlat;
        goto end_north;
      }
      ++n;
    }
  } end_north:

  n=rg.npts()-1;
  int lat_south=-1;
  for( int jlat=rg.ny()-1; jlat>=0; --jlat) {
    for( int jlon=rg.nx(jlat)-1; jlon>=0; --jlon) {
      if( parts.at(n) == mypart ) {
        lat_south=jlat;
        goto end_south;
      }
      --n;
    }
  } end_south:


  std::vector<int> offset(rg.ny(),0);

  n=0;
  for(size_t jlat = 0; jlat < rg.ny(); ++jlat)
  {
    offset.at(jlat)=n;
    n+=rg.nx(jlat);
  };

  /*
  We need to connect to next region
  */
  if( lat_north-1 >=0             && rg.nx(lat_north-1) > 0 )
    --lat_north;
  if(size_t(lat_south+1) <= rg.ny()-1 && rg.nx(lat_south+1) > 0 )
    ++lat_south;
  region.lat_begin.resize(rg.ny(),-1);
  region.lat_end.resize(rg.ny(),-1);
  region.nb_lat_elems.resize(rg.ny(),0);
  region.north = lat_north;
  region.south = lat_south;

  array::ArrayShape shape = array::make_shape(region.south-region.north, 4*rg.nxmax(), 4);

  region.elems.resize(shape);
  region.elems = -1;

  int nelems=0;
  region.nquads=0;
  region.ntriags=0;

  array::ArrayView<int,3> elemview(region.elems);

  bool stagger = options.get<bool>("stagger");
  for (int jlat=lat_north; jlat<lat_south; ++jlat)
  {
//    std::stringstream filename; filename << "/tmp/debug/"<<jlat;

    size_t ilat, latN, latS;
    size_t ipN1, ipN2, ipS1, ipS2;
    double xN1, xN2, yN, xS1, xS2, yS;
    double dN1S2, dS1N2;  // dN2S2;
    bool try_make_triangle_up, try_make_triangle_down, try_make_quad;
    bool add_triag, add_quad;

    ilat = jlat-region.north;

    array::ArrayView<int,2> lat_elems_view = elemview.at(ilat);

    latN = jlat;
    latS = jlat+1;
    yN = rg.y(latN);
    yS = rg.y(latS);

    size_t beginN, beginS, endN, endS;

    beginN = 0;
    endN   = rg.nx(latN) - (periodic_east_west ? 0 : 1);
    if( yN == 90 && unique_pole )
      endN = beginN;


    beginS = 0;
    endS   = rg.nx(latS) - (periodic_east_west ? 0 : 1);
    if( yS == -90 && unique_pole )
      endS = beginS;

    ipN1 = beginN;
    ipS1 = beginS;
    ipN2 = std::min(ipN1+1,endN);
    ipS2 = std::min(ipS1+1,endS);

    int jelem=0;

#if DEBUG_OUTPUT
    Log::info()  << "=================\n";
    Log::info()  << "latN, latS : " << latN << ", " << latS << '\n';
#endif

    while (true)
    {

      if( ipN1 == endN && ipS1 == endS ) break;

#if DEBUG_OUTPUT
      Log::info(Here())  << "-------\n";
#endif

      //ASSERT(offset.at(latN)+ipN1 < parts.size());
      //ASSERT(offset.at(latS)+ipS1 < parts.size());


      int pN1, pS1, pN2, pS2;
      if( ipN1 != rg.nx(latN) )
        pN1 = parts.at(offset.at(latN)+ipN1);
      else
        pN1 = parts.at(offset.at(latN)+ipN1-1);
      if( ipS1 != rg.nx(latS) )
        pS1 = parts.at(offset.at(latS)+ipS1);
      else
        pS1 = parts.at(offset.at(latS)+ipS1-1);


      if( ipN2 == rg.nx(latN) )
        pN2 = pN1;
      else
        pN2 = parts.at(offset.at(latN)+ipN2);
      if( ipS2 == rg.nx(latS) )
        pS2 = pS1;
      else
        pS2 = parts.at(offset.at(latS)+ipS2);

      //Log::info()  << ipN1 << "("<<pN1<<") " << ipN2 <<"("<<pN2<<")" <<  std::endl;
      //Log::info()  << ipS1 << "("<<pS2<<") " << ipS2 <<"("<<pS2<<")" <<  std::endl;

#if DEBUG_OUTPUT
      Log::info()  << ipN1 << "("<<pN1<<") " << ipN2 <<"("<<pN2<<")" <<  std::endl;
      Log::info()  << ipS1 << "("<<pS2<<") " << ipS2 <<"("<<pS2<<")" <<  std::endl;
#endif

      xN1 = rg.x(ipN1,latN) * to_rad;
      xN2 = rg.x(ipN2,latN) * to_rad;
      xS1 = rg.x(ipS1,latS) * to_rad;
      xS2 = rg.x(ipS2,latS) * to_rad;

      if( stagger && (latN+1)%2==0 ) xN1 += M_PI/static_cast<double>(rg.nx(latN));
      if( stagger && (latN+1)%2==0 ) xN2 += M_PI/static_cast<double>(rg.nx(latN));
      if( stagger && (latS+1)%2==0 ) xS1 += M_PI/static_cast<double>(rg.nx(latS));
      if( stagger && (latS+1)%2==0 ) xS2 += M_PI/static_cast<double>(rg.nx(latS));

#if DEBUG_OUTPUT
      Log::info()  << "-------\n";
#endif
      // Log::info()  << "  access  " << region.elems.stride(0)*(jlat-region.north) + region.elems.stride(1)*jelem + 5 << std::endl;
//      Log::info()  << ipN1 << "("<< xN1 << ")  " << ipN2 <<  "("<< xN2 << ")  " << std::endl;
//      Log::info()  << ipS1 << "("<< xS1 << ")  " << ipS2 <<  "("<< xS2 << ")  " << std::endl;
      try_make_triangle_up   = false;
      try_make_triangle_down = false;
      try_make_quad = false;


// ------------------------------------------------
// START RULES
// ------------------------------------------------

      const double dxN = std::abs(xN2-xN1);
      const double dxS = std::abs(xS2-xS1);
      const double dx  = std::min(dxN,dxS);
      const double alpha1 = ( dx==0. ? 0. : std::atan2((xN1-xS1)/dx,1.) * to_deg );
      const double alpha2 = ( dx==0. ? 0. : std::atan2((xN2-xS2)/dx,1.) * to_deg );
      if( std::abs(alpha1) <= max_angle && std::abs(alpha2) <= max_angle )
      {
        if( triangulate_quads )
        {
          if( false ) //std::abs(alpha1) < 1 && std::abs(alpha2) < 1)
          {
            try_make_triangle_up   = (jlat+ipN1) % 2;
            try_make_triangle_down = (jlat+ipN1+1) % 2;
          }
          else
          {
            dN1S2 = std::abs(xN1-xS2);
            dS1N2 = std::abs(xS1-xN2);
            // dN2S2 = std::abs(xN2-xS2);
            // Log::info()  << "  dN1S2 " << dN1S2 << "   dS1N2 " << dS1N2 << "   dN2S2 " << dN2S2 << std::endl;
            if (dN1S2 == dS1N2)
            {
              try_make_triangle_up   = (jlat+ipN1) % 2;
              try_make_triangle_down = (jlat+ipN1+1) % 2;
            }
            else if (dN1S2 < dS1N2)
            {
              if (ipS1 != ipS2) { try_make_triangle_up = true; }
              else { try_make_triangle_down = true ; }
            }
            else if (dN1S2 > dS1N2)
            {
              if (ipN1 != ipN2) { try_make_triangle_down = true;}
              else { try_make_triangle_up = true; }
            }
            else
            {
              throw eckit::Exception("Should not be here", Here());
            }
          }

        }
        else
        {
          if     ( ipN1 == ipN2 ) try_make_triangle_up   = true;
          else if( ipS1 == ipS2 ) try_make_triangle_down = true;
          else                    try_make_quad          = true;

//          try_make_quad          = true;
        }
      }
      else
      {
        dN1S2 = std::abs(xN1-xS2);
        dS1N2 = std::abs(xS1-xN2);
        // dN2S2 = std::abs(xN2-xS2);
        // Log::info()  << "  dN1S2 " << dN1S2 << "   dS1N2 " << dS1N2 << "   dN2S2 " << dN2S2 << std::endl;
        if( (dN1S2 <= dS1N2) && (ipS1 != ipS2) ) { try_make_triangle_up = true;}
        else if( (dN1S2 >= dS1N2) && (ipN1 != ipN2) ) { try_make_triangle_down = true;}
        else Exception("Should not try to make a quadrilateral!",Here());
      }
// ------------------------------------------------
// END RULES
// ------------------------------------------------


#if DEBUG_OUTPUT
      DEBUG_VAR(jelem);
#endif

      array::ArrayView<int,1> elem = lat_elems_view.at(jelem);

      if( try_make_quad )
      {
        // add quadrilateral
#if DEBUG_OUTPUT
        Log::info()  << "          " << ipN1 << "  " << ipN2 << '\n';
        Log::info()  << "          " << ipS1 << "  " << ipS2 << '\n';
#endif
        elem.at(0) = ipN1;
        elem.at(1) = ipS1;
        elem.at(2) = ipS2;
        elem.at(3) = ipN2;
        add_quad = false;
        int np[] = {pN1, pN2, pS1, pS2};
        int cnt_mypart = std::count(np, np+4, mypart);
        if( cnt_mypart > 0 )
        {
          int pcnts[4];
          for( int j=0; j<4; ++j )
            pcnts[j] = std::count(np, np+4, np[j]);
          int cnt_max = *std::max_element(pcnts,pcnts+4);

          if( latN == 0 )
          {
            if( pN1 == mypart )
            {
              add_quad = true;
            }
          }
          else if ( latS == rg.ny()-1 )
          {
            if( pS2 == mypart )
            {
              add_quad = true;
            }
          }
          else if( cnt_mypart > 2 ) // 3 or 4 points belong to mypart
          {
            add_quad = true;
          }
          else if( cnt_max < 3 ) // 3 or 4 points don't belong to mypart
          {
            if( 0.5*(yN+yS) > 1e-6)
            {
              if ( pS1 == mypart )  add_quad = true;
            }
            else
            {
              if ( pN2 == mypart )  add_quad = true;
            }
          }
        }
        if (add_quad)
        {
          ++region.nquads;
          ++jelem;
          ++nelems;

          if( region.lat_begin.at(latN) == -1 ) region.lat_begin.at(latN) = ipN1;
          if( region.lat_begin.at(latS) == -1 ) region.lat_begin.at(latS) = ipS1;
          region.lat_begin.at(latN) = std::min<int>(region.lat_begin.at(latN), ipN1);
          region.lat_begin.at(latS) = std::min<int>(region.lat_begin.at(latS), ipS1);
          region.lat_end.at(latN) = std::max<int>(region.lat_end.at(latN), ipN2);
          region.lat_end.at(latS) = std::max<int>(region.lat_end.at(latS), ipS2);
        }
        else
        {
#if DEBUG_OUTPUT
          Log::info() << "Quad belongs to other partition" << std::endl;
#endif
        }
        ipN1=ipN2;
        ipS1=ipS2;
      }
      else if( try_make_triangle_down  ) // make triangle down
      {
        // triangle without ip3
#if DEBUG_OUTPUT
        Log::info()  << "          " << ipN1 << "  " << ipN2 << '\n';
        Log::info()  << "          " << ipS1 << '\n';
#endif
        elem.at(0) = ipN1;
        elem.at(1) = ipS1;
        elem.at(2) = -1;
        elem.at(3) = ipN2;

        add_triag = false;

        int cnt_mypart = 0;
        int np[3] = {pN1, pN2, pS1};
        for( int j=0; j<3; ++j ) {
          if (np[j]==mypart) ++cnt_mypart;
        }
        if( cnt_mypart > 1 )
        {
          add_triag=true;
        }
        else if( cnt_mypart == 1)
        {
          int pcnts[3];
          for( int j=0; j<3; ++j )
            pcnts[j] = std::count(np, np+3, np[j]);
          int cnt_max = *std::max_element(pcnts,pcnts+3);
          if( cnt_max == 1)
          {
            if( 0.5*(yN+yS) > 1e-6 )
            {
              if( pS1 == mypart )
                add_triag = true;
            }
            else
            {
              if( pN1 == mypart )
                add_triag = true;
            }
          }
        }
        if (add_triag)
        {
          ++region.ntriags;
          ++jelem;
          ++nelems;

          if( region.lat_begin.at(latN) == -1 ) region.lat_begin.at(latN) = ipN1;
          if( region.lat_begin.at(latS) == -1 ) region.lat_begin.at(latS) = ipS1;
          region.lat_begin.at(latN) = std::min<int>(region.lat_begin.at(latN), ipN1);
          region.lat_begin.at(latS) = std::min<int>(region.lat_begin.at(latS), ipS1);
          region.lat_end.at(latN) = std::max<int>(region.lat_end.at(latN), ipN2);
          region.lat_end.at(latS) = std::max<int>(region.lat_end.at(latS), ipS1);
        }
        else
        {
#if DEBUG_OUTPUT
          Log::info() << "Downward Triag belongs to other partition" << std::endl;
#endif
        }
        ipN1=ipN2;
        // and ipS1=ipS1;

      }
      else if( try_make_triangle_up ) // make triangle up
      {
        // triangle without ip4
#if DEBUG_OUTPUT
        Log::info()  << "          " << ipN1 << " ("<<pN1<<")" << '\n';
        Log::info()  << "          " << ipS1 << " ("<<pS1<<")" << "  " << ipS2 << " ("<<pS2<<")" << '\n';
#endif
        elem.at(0) = ipN1;
        elem.at(1) = ipS1;
        elem.at(2) = ipS2;
        elem.at(3) = -1;

        add_triag = false;


        int cnt_mypart = 0;
        int np[3] = {pN1, pS1, pS2};
        for( int j=0; j<3; ++j ) {
          if (np[j]==mypart) ++cnt_mypart;
        }

        if( cnt_mypart > 1 )
        {
          add_triag=true;
        }
        else if( cnt_mypart == 1)
        {
          int pcnts[3];
          for( int j=0; j<3; ++j )
            pcnts[j] = std::count(np, np+3, np[j]);
          int cnt_max = *std::max_element(pcnts,pcnts+3);
          if( cnt_max == 1)
          {
            // if( latN == 0 && mypart = pN1 ) || latS == rg.nlat()-1 )
            // {
            //   add_triag = true;
            // }
            // else
            if( 0.5*(yN+yS) > 1e-6 )
            {
              if( pS2 == mypart )
                add_triag = true;
            }
            else
            {
              if( pN1 == mypart )
                add_triag = true;
            }
          }
        }

        if (add_triag)
        {
          ++region.ntriags;
          ++jelem;
          ++nelems;

          if( region.lat_begin.at(latN) == -1 ) region.lat_begin.at(latN) = ipN1;
          if( region.lat_begin.at(latS) == -1 ) region.lat_begin.at(latS) = ipS1;
          region.lat_begin.at(latN) = std::min<int>(region.lat_begin.at(latN), ipN1);
          region.lat_begin.at(latS) = std::min<int>(region.lat_begin.at(latS), ipS1);
          region.lat_end.at(latN) = std::max<int>(region.lat_end.at(latN), ipN1);
          region.lat_end.at(latS) = std::max<int>(region.lat_end.at(latS), ipS2);
        }
        else
        {
#if DEBUG_OUTPUT
          Log::info() << "Upward Triag belongs to other partition" << std::endl;
#endif
        }
        ipS1=ipS2;
        // and ipN1=ipN1;

      }
      else
      {
        throw eckit::SeriousBug("Could not detect which element to create", Here());
      }
      ipN2 = std::min(endN,ipN1+1);
      ipS2 = std::min(endS,ipS1+1);
    }
    region.nb_lat_elems.at(jlat) = jelem;
#if DEBUG_OUTPUT
    DEBUG_VAR(region.nb_lat_elems.at(jlat));
#endif
    if( region.nb_lat_elems.at(jlat) == 0 && latN == size_t(region.north) ) {
      ++region.north;
    }
    if( region.nb_lat_elems.at(jlat) == 0 && latS == size_t(region.south) ) {
      --region.south;
    }
    region.lat_end.at(latN) = std::min(region.lat_end.at(latN), int(rg.nx(latN)-1));
    region.lat_end.at(latS) = std::min(region.lat_end.at(latS), int(rg.nx(latS)-1));
    if( yN == 90 && unique_pole )
      region.lat_end.at(latN) = rg.nx(latN)-1;
    if( yS == -90 && unique_pole )
      region.lat_end.at(latS) = rg.nx(latS)-1;

    if( region.nb_lat_elems.at(jlat) > 0 )
    {
      region.lat_end.at(latN) = std::max(region.lat_end.at(latN), region.lat_begin.at(latN));
      region.lat_end.at(latS) = std::max(region.lat_end.at(latS), region.lat_begin.at(latS));
    }
  } // for jlat

//  Log::info()  << "nb_triags = " << region.ntriags << std::endl;
//  Log::info()  << "nb_quads = " << region.nquads << std::endl;
//  Log::info()  << "nb_elems = " << nelems << std::endl;

  int nb_region_nodes = 0;
  for(int jlat = region.north; jlat <= region.south; ++jlat) {
    n = offset.at(jlat);
    region.lat_begin.at(jlat) = std::max( 0, region.lat_begin.at(jlat) );
    for(size_t jlon = 0; jlon < rg.nx(jlat); ++jlon) {
      if( parts.at(n) == mypart ) {
        region.lat_begin.at(jlat) = std::min( region.lat_begin.at(jlat), int(jlon) );
        region.lat_end.  at(jlat) = std::max( region.lat_end  .at(jlat), int(jlon) );
      }
      ++n;
    }
    nb_region_nodes += region.lat_end.at(jlat)-region.lat_begin.at(jlat)+1;

    // Count extra periodic node
    if( periodic_east_west && size_t(region.lat_end.at(jlat)) == rg.nx(jlat) - 1) ++nb_region_nodes;
  }

  region.nnodes = nb_region_nodes;
  if (region.nnodes == 0)
  {
    throw Exception("Trying to generate mesh with too many partitions. Reduce the number of partitions.",Here());
  }
#if DEBUG_OUTPUT
  DEBUG("End of generate_region()");
#endif
}


namespace {
struct GhostNode {
  GhostNode(int _jlat, int _jlon, int _jnode) { jlat = _jlat; jlon=_jlon; jnode=_jnode; }
  int jlat;
  int jlon;
  int jnode;
};
}

void StructuredMeshGenerator::generate_mesh(const grid::StructuredGrid& rg, const std::vector<int>& parts, const Region& region, Mesh& mesh) const
{


  ASSERT(!mesh.generated());

  int mypart = options.get<size_t>("part");
  int nparts = options.get<size_t>("nb_parts");
  int n, l;

  bool has_north_pole = rg.y().front() ==  90 && rg.nx().front() > 0;
  bool has_south_pole = rg.y().back()  == -90 && rg.nx().back()  > 0;
  bool three_dimensional  = options.get<bool>("3d");
  bool periodic_east_west = rg.periodic();
  bool include_periodic_ghost_points = periodic_east_west  && !three_dimensional;
  bool remove_periodic_ghost_points =  periodic_east_west  &&  three_dimensional;

  bool include_north_pole = (mypart == 0 )
          && options.get<bool>("include_pole")
          && !has_north_pole
          && rg.domain().includesNorthPole(rg.projection());

  bool include_south_pole = (mypart == nparts-1)
          && options.get<bool>("include_pole")
          && !has_south_pole
          && rg.domain().includesSouthPole(rg.projection());

  bool patch_north_pole = (mypart == 0)
          && options.get<bool>("patch_pole")
          && !has_north_pole
          && rg.domain().includesNorthPole(rg.projection())
          && rg.nx(1) > 0;

  bool patch_south_pole = (mypart == nparts-1)
          && options.get<bool>("patch_pole")
          && !has_south_pole
          && rg.domain().includesSouthPole(rg.projection())
          && rg.nx(rg.ny()-2) > 0;

  if( three_dimensional && nparts != 1 )
    throw BadParameter("Cannot generate three_dimensional mesh in parallel",Here());
  int nnodes  = region.nnodes;
  int ntriags = region.ntriags;
  int nquads  = region.nquads;

  if (include_north_pole) {
    ++nnodes;
    ntriags += rg.nx(0) - (periodic_east_west ? 0 : 1);
  }
  else if (patch_north_pole) {
    ntriags += rg.nx(0)-2;
  }
  if (include_south_pole) {
    ++nnodes;
    ntriags += rg.nx(rg.ny()-1) - (periodic_east_west ? 0 : 1);
  }
  else if (patch_south_pole) {
    ntriags += rg.nx(rg.ny()-1)-2;
  }

  size_t node_numbering_size = nnodes;
  if ( remove_periodic_ghost_points ) {
    for(size_t jlat = 0; jlat < rg.ny(); ++jlat)
    {
      if( rg.nx(jlat) > 0 )
        --nnodes;
    }
  }

#if DEBUG_OUTPUT
  DEBUG_VAR(include_periodic_ghost_points);
  DEBUG_VAR(include_north_pole);
  DEBUG_VAR(include_south_pole);
  DEBUG_VAR(three_dimensional);
  DEBUG_VAR(patch_north_pole);
  DEBUG_VAR(patch_south_pole);
  DEBUG_VAR(rg.npts());
  DEBUG_VAR(nnodes);
  DEBUG_VAR(ntriags);
  DEBUG_VAR(nquads);
  DEBUG_VAR(options.get<bool>("ghost_at_end"));
#endif


  std::vector<int> offset_glb(rg.ny());
  std::vector<int> offset_loc(region.south-region.north+1,0);

  n=0;
  for(size_t jlat = 0; jlat < rg.ny(); ++jlat)
  {
    offset_glb.at(jlat)=n;
    n+=rg.nx(jlat);
  };

  std::vector<int> periodic_glb(rg.ny());

  if( include_periodic_ghost_points )
  {
    for(size_t jlat = 0; jlat < rg.ny(); ++jlat)
    {
      if( rg.nx(jlat) > 0 )
      {
        periodic_glb.at(jlat) = n;
        ++n;
      }
    }
  }
  else
  {
    for(size_t jlat = 0; jlat < rg.ny(); ++jlat)
    {
      if( rg.nx(jlat) > 0 )
      {
        periodic_glb.at(jlat) = offset_glb.at(jlat) + rg.nx(jlat)-1;
      }
    }
  }


  mesh.nodes().resize(nnodes);
  mesh::Nodes& nodes = mesh.nodes();

  array::ArrayView<double,2> lonlat        ( nodes.lonlat() );
  array::ArrayView<double,2> geolonlat     ( nodes.geolonlat() );
  array::ArrayView<gidx_t,1> glb_idx       ( nodes.global_index() );
  array::ArrayView<int,   1> part          ( nodes.partition() );
  array::ArrayView<int,   1> ghost         ( nodes.ghost() );
  array::ArrayView<int,   1> flags         ( nodes.field("flags") );

  bool stagger = options.get<bool>("stagger");

  std::vector<int> node_numbering(node_numbering_size,-1);
  if( options.get<bool>("ghost_at_end") )
  {
    std::vector< GhostNode > ghost_nodes;
    ghost_nodes.reserve( nnodes );
    int node_number=0;
    int jnode=0;
    l=0;
    ASSERT( region.south >= region.north );
    for( int jlat = region.north; jlat <= region.south; ++jlat )
    {
      int ilat=jlat-region.north;
      offset_loc.at(ilat)=l;
      l+=region.lat_end.at(jlat)-region.lat_begin.at(jlat)+1;


      if( region.lat_end.at(jlat) < region.lat_begin.at(jlat) )
      {
        DEBUG_VAR(jlat);
        DEBUG_VAR(region.lat_begin[jlat]);
        DEBUG_VAR(region.lat_end[jlat]);
      }
      for( int jlon=region.lat_begin.at(jlat); jlon<=region.lat_end.at(jlat); ++jlon )
      {
        n = offset_glb.at(jlat) + jlon;
        if( parts.at(n) == mypart ) {
          node_numbering.at(jnode) = node_number;
          ++node_number;
        }
        else {
          ghost_nodes.push_back( GhostNode(jlat,jlon,jnode));
        }
        ++jnode;
      }
      if(include_periodic_ghost_points && size_t(region.lat_end.at(jlat)) == rg.nx(jlat) - 1) // add periodic point
      {
        ++l;
        part(jnode)      = part(jnode-1);
        ghost(jnode)     = 1;
        ghost_nodes.push_back( GhostNode(jlat,rg.nx(jlat),jnode) );
        ++jnode;
      }
    }
    for( size_t jghost=0; jghost<ghost_nodes.size(); ++jghost )
    {
      node_numbering.at( ghost_nodes.at(jghost).jnode ) = node_number;
      ++node_number;
    }
    if (include_north_pole)
    {
      node_numbering.at( jnode ) = jnode;
      ++jnode;
    }
    if (include_south_pole)
    {
      node_numbering.at( jnode ) = jnode;
      ++jnode;
    }
    ASSERT( jnode == nnodes );
  }
  else // No renumbering
  {
    for( int jnode=0; jnode<nnodes; ++jnode )
      node_numbering.at(jnode) = jnode;
  }

  int jnode=0;
  l=0;
  for(int jlat = region.north; jlat <= region.south; ++jlat)
  {
    int ilat=jlat-region.north;
    offset_loc.at(ilat)=l;
    l+=region.lat_end.at(jlat)-region.lat_begin.at(jlat)+1;

    double y = rg.y(jlat);
    for( int jlon=region.lat_begin.at(jlat); jlon<=region.lat_end.at(jlat); ++jlon )
    {
      int inode = node_numbering.at(jnode);
      n = offset_glb.at(jlat) + jlon;

      double x = rg.x(jlon,jlat);
      //std::cout << "jlat = " << jlat << "; jlon = " << jlon << "; x = " << x << std::endl;
      if( stagger && (jlat+1)%2==0 ) x += 180./static_cast<double>(rg.nx(jlat));

      lonlat(inode,LON) = x;
      lonlat(inode,LAT) = y;

      // geographic coordinates by using projection
      double crd[] = {x,y};
      rg.projection().xy2lonlat(crd);
      geolonlat(inode,LON) = crd[LON];
      geolonlat(inode,LAT) = crd[LAT];

      glb_idx(inode)   = n+1;
      part(inode) = parts.at(n);
      ghost(inode) = 0;
      Topology::reset(flags(inode));
      if( jlat == 0 && !include_north_pole) {
        Topology::set(flags(inode),Topology::BC|Topology::NORTH);
      }
      if( size_t(jlat) == rg.ny()-1 && !include_south_pole) {
        Topology::set(flags(inode),Topology::BC|Topology::SOUTH);
      }
      if( jlon == 0 && include_periodic_ghost_points) {
        Topology::set(flags(inode),Topology::BC|Topology::WEST);
      }
      if( part(inode) != mypart ) {
        Topology::set(flags(inode),Topology::GHOST);
        ghost(inode) = 1;
      }
      ++jnode;
    }
    if(include_periodic_ghost_points && size_t(region.lat_end.at(jlat)) == rg.nx(jlat) - 1) // add periodic point
    {
      int inode = node_numbering.at(jnode);
      int inode_left = node_numbering.at(jnode-1);
      ++l;
      double x = rg.x(rg.nx(jlat),jlat);
      if( stagger && (jlat+1)%2==0 ) x += 180./static_cast<double>(rg.nx(jlat));

      lonlat(inode,LON) = x;
      lonlat(inode,LAT) = y;

      // geographic coordinates by using projection
      double crd[] = {x,y};
      rg.projection().xy2lonlat(crd);
      geolonlat(inode,LON) = crd[LON];
      geolonlat(inode,LAT) = crd[LAT];


      glb_idx(inode)   = periodic_glb.at(jlat)+1;
      part(inode)      = part(inode_left);
      ghost(inode)     = 1;
      Topology::reset(flags(inode));
      Topology::set(flags(inode),Topology::BC|Topology::EAST);
      Topology::set(flags(inode),Topology::GHOST);
      ++jnode;
    }
  };

  int jnorth=-1;
  if (include_north_pole)
  {
    int inode = node_numbering.at(jnode);
    jnorth = jnode;
    double y = 90.;
    double x = 180.;
    lonlat(inode,LON) = x;
    lonlat(inode,LAT) = y;

    // geographic coordinates by using projection
    double crd[] = {x,y};
    rg.projection().xy2lonlat(crd);
    geolonlat(inode,LON) = crd[LON];
    geolonlat(inode,LAT) = crd[LAT];

    glb_idx(inode)   = periodic_glb.at(rg.ny()-1)+2;
    part(inode)      = mypart;
    ghost(inode)     = 0;
    Topology::reset(flags(inode));
    Topology::set(flags(inode),Topology::NORTH);
    ++jnode;
  }

  int jsouth=-1;
  if (include_south_pole)
  {
    int inode = node_numbering.at(jnode);
    jsouth = jnode;
    double y = -90.;
    double x =  180.;
    lonlat(inode,LON) = x;
    lonlat(inode,LAT) = y;

    // geographic coordinates by using projection
    double crd[] = {x,y};
    rg.projection().xy2lonlat(crd);
    geolonlat(inode,LON) = crd[LON];
    geolonlat(inode,LAT) = crd[LAT];

    glb_idx(inode)   = periodic_glb.at(rg.ny()-1)+3;
    part(inode)      = mypart;
    ghost(inode)     = 0;
    Topology::reset(flags(inode));
    Topology::set(flags(inode),Topology::SOUTH);
    ++jnode;
  }

  mesh.cells().add( new mesh::temporary::Quadrilateral(), nquads  );
  mesh.cells().add( new mesh::temporary::Triangle(),      ntriags );

  mesh::HybridElements::Connectivity& node_connectivity = mesh.cells().node_connectivity();
  array::ArrayView<gidx_t,1> cells_glb_idx( mesh.cells().global_index() );
  array::ArrayView<int,1>    cells_part(    mesh.cells().partition() );
  array::ArrayView<int,1>    cells_patch(   mesh.cells().field("patch") );

  /*
   * label all patch cells a non-patch
   */
  cells_patch = 0;

  /*
  Fill in connectivity tables with global node indices first
  */
  int jcell;
  int jquad=0;
  int jtriag=0;
  int quad_begin  = mesh.cells().elements(0).begin();
  int triag_begin = mesh.cells().elements(1).begin();
  int quad_nodes[4];
  int triag_nodes[3];

  for (int jlat=region.north; jlat<region.south; ++jlat)
  {

    int ilat = jlat-region.north;
    int jlatN = jlat;
    int jlatS = jlat+1;
    int ilatN = ilat;
    int ilatS = ilat+1;
    for (int jelem=0; jelem<region.nb_lat_elems.at(jlat); ++jelem)
    {
      const array::ArrayView<int,1> elem = array::ArrayView<int,3>(region.elems).at(ilat).at(jelem);

      if(elem.at(2)>=0 && elem.at(3)>=0) // This is a quad
      {
        quad_nodes[0] = node_numbering.at( offset_loc.at(ilatN) + elem.at(0) - region.lat_begin.at(jlatN) );
        quad_nodes[1] = node_numbering.at( offset_loc.at(ilatS) + elem.at(1) - region.lat_begin.at(jlatS) );
        quad_nodes[2] = node_numbering.at( offset_loc.at(ilatS) + elem.at(2) - region.lat_begin.at(jlatS) );
        quad_nodes[3] = node_numbering.at( offset_loc.at(ilatN) + elem.at(3) - region.lat_begin.at(jlatN) );

        if( three_dimensional && periodic_east_west )
        {
          if (size_t(elem.at(2)) == rg.nx(jlatS)) quad_nodes[2] = node_numbering.at( offset_loc.at(ilatS) );
          if (size_t(elem.at(3)) == rg.nx(jlatN)) quad_nodes[3] = node_numbering.at( offset_loc.at(ilatN) );
        }

        jcell = quad_begin + jquad++;
        node_connectivity.set( jcell, quad_nodes );
        cells_glb_idx(jcell) = jcell+1;
        cells_part(jcell)    = mypart;
      }
      else // This is a triag
      {
        if(elem.at(3)<0) // This is a triangle pointing up
        {
          triag_nodes[0] = node_numbering.at( offset_loc.at(ilatN) + elem.at(0) - region.lat_begin.at(jlatN) );
          triag_nodes[1] = node_numbering.at( offset_loc.at(ilatS) + elem.at(1) - region.lat_begin.at(jlatS) );
          triag_nodes[2] = node_numbering.at( offset_loc.at(ilatS) + elem.at(2) - region.lat_begin.at(jlatS) );
          if( three_dimensional && periodic_east_west )
          {
            if (size_t(elem.at(0)) == rg.nx(jlatN)) triag_nodes[0] = node_numbering.at( offset_loc.at(ilatN) );
            if (size_t(elem.at(2)) == rg.nx(jlatS)) triag_nodes[2] = node_numbering.at( offset_loc.at(ilatS) );
          }
        }
        else // This is a triangle pointing down
        {
          triag_nodes[0] = node_numbering.at( offset_loc.at(ilatN) + elem.at(0) - region.lat_begin.at(jlatN) );
          triag_nodes[1] = node_numbering.at( offset_loc.at(ilatS) + elem.at(1) - region.lat_begin.at(jlatS) );
          triag_nodes[2] = node_numbering.at( offset_loc.at(ilatN) + elem.at(3) - region.lat_begin.at(jlatN) );
          if( three_dimensional && periodic_east_west )
          {
            if (size_t(elem.at(1)) == rg.nx(jlatS)) triag_nodes[1] = node_numbering.at( offset_loc.at(ilatS) );
            if (size_t(elem.at(3)) == rg.nx(jlatN)) triag_nodes[2] = node_numbering.at( offset_loc.at(ilatN) );
          }
        }
        jcell = triag_begin + jtriag++;
        node_connectivity.set( jcell, triag_nodes );
        cells_glb_idx(jcell) = jcell+1;
        cells_part(jcell) = mypart;
      }
    }
  }

  if (include_north_pole)
  {
    int ilat = 0;
    int ip1  = 0;
    size_t nlon = rg.nx(0) - (periodic_east_west ? 0 : 1);
    for (size_t ip2 = 0; ip2 < nlon; ++ip2)
    {
      jcell = triag_begin + jtriag++;
      size_t ip3 = ip2+1;
      if( three_dimensional && ip3 == rg.nx(0) ) ip3 = 0;
      triag_nodes[0] = node_numbering.at( jnorth           + ip1 );
      triag_nodes[1] = node_numbering.at( offset_loc.at(ilat) + ip2 );
      triag_nodes[2] = node_numbering.at( offset_loc.at(ilat) + ip3 );
      node_connectivity.set( jcell, triag_nodes );
      cells_glb_idx(jcell) = jcell+1;
      cells_part(jcell) = mypart;
    }
  }
  else if (patch_north_pole )
  {
    int jlat = 0;
    int ilat = 0;
    int ip1, ip2, ip3;

    int jforward = 0;
    int jbackward = rg.nx(jlat) - 1;
    bool forward = true;

    while( true )
    {
      if( forward )
      {
        ip1 = jforward;
        ip2 = jforward+1;
        ip3 = jbackward;
      }
      else
      {
        ip1 = jforward;
        ip2 = jbackward-1;
        ip3 = jbackward;
      }

      triag_nodes[0] = node_numbering.at( offset_loc.at(ilat) + ip1 );
      triag_nodes[1] = node_numbering.at( offset_loc.at(ilat) + ip2 );
      triag_nodes[2] = node_numbering.at( offset_loc.at(ilat) + ip3 );

      jcell = triag_begin + jtriag++;
      node_connectivity.set( jcell, triag_nodes );

      cells_glb_idx (jcell) = jcell+1;
      cells_part    (jcell) = mypart;
      cells_patch   (jcell) = 1;  // mark cell as "patch"

      if (jbackward == jforward+2 ) break;

      if( forward )
      {
        ++jforward;
        forward = false;
      }
      else
      {
        --jbackward;
        forward = true;
      }
    }
  }

  if (include_south_pole)
  {
    int jlat = rg.ny()-1;
    int ilat = region.south-region.north;
    int ip1 = 0;
    size_t nlon = rg.nx(jlat)+1 - (periodic_east_west ? 0 : 1);
    for (size_t ip2 = 1; ip2 < nlon; ++ip2)
    {
      jcell = triag_begin + jtriag++;
      int ip3 = ip2-1;
      triag_nodes[0] = node_numbering.at( jsouth           + ip1 );
      triag_nodes[1] = node_numbering.at( offset_loc.at(ilat) + ip2 );
      triag_nodes[2] = node_numbering.at( offset_loc.at(ilat) + ip3 );
      if( three_dimensional && ip2 == rg.nx(jlat) )
        triag_nodes[1] = node_numbering.at( offset_loc.at(ilat) + 0 );
      node_connectivity.set( jcell, triag_nodes );
      cells_glb_idx(jcell) = jcell+1;
      cells_part(jcell) = mypart;
    }
  }
  else if (patch_south_pole)
  {
    int jlat = rg.ny()-1;
    int ilat = region.south-region.north;
    int ip1, ip2, ip3;

    int jforward = 0;
    int jbackward = rg.nx(jlat) - 1;
    bool forward = true;

    while( true )
    {
      if( forward )
      {
        ip1 = jforward;
        ip2 = jforward+1;
        ip3 = jbackward;
      }
      else
      {
        ip1 = jforward;
        ip2 = jbackward-1;
        ip3 = jbackward;
      }

      triag_nodes[0] = node_numbering.at( offset_loc.at(ilat) + ip1 );
      triag_nodes[1] = node_numbering.at( offset_loc.at(ilat) + ip2 );
      triag_nodes[2] = node_numbering.at( offset_loc.at(ilat) + ip3 );

      jcell = triag_begin + jtriag++;
      node_connectivity.set( jcell, triag_nodes );

      cells_glb_idx (jcell) = jcell+1;
      cells_part    (jcell) = mypart;
      cells_patch   (jcell) = 1;  // mark cell as "patch"

      if (jbackward == jforward+2 ) break;

      if( forward )
      {
        ++jforward;
        forward = false;
      }
      else
      {
        --jbackward;
        forward = true;
      }
    }
  }

  generate_global_element_numbering( mesh );
}

namespace {
static MeshGeneratorBuilder< StructuredMeshGenerator > __Structured("structured");
}

} // namespace meshgenerator
} // namespace atlas
