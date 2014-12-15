/*
 * (C) Copyright 1996-2014 ECMWF.
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

#include <eckit/runtime/Context.h>
#include <eckit/config/Configurable.h>
#include <eckit/config/Resource.h>

#include "atlas/util/Array.h"
#include "atlas/util/ArrayView.h"
#include "atlas/util/IndexView.h"
#include "atlas/Field.h"
#include "atlas/FunctionSpace.h"
#include "atlas/Mesh.h"
#include "atlas/Util.h"
#include "atlas/Parameters.h"
#include "atlas/meshgen/EqualAreaPartitioner.h"
#include "atlas/grids/ReducedGrid.h"
#include "atlas/meshgen/RGGMeshGenerator.h"

#define DEBUG_OUTPUT 0

using namespace eckit;
using namespace atlas::grids;

namespace atlas {
namespace meshgen {

namespace {
static double to_rad = M_PI/180.;
static double to_deg = 180.*M_1_PI;
}

struct Region
{
  int north;
  int south;
  Array<int> elems;
  int ntriags;
  int nquads;
  int nnodes;
  std::vector<int> lat_begin;
  std::vector<int> lat_end;
  std::vector<int> nb_lat_elems;
};

ReducedGridMeshGenerator::ReducedGridMeshGenerator()
{
  // This option creates a point at the pole when true
  options.set("include_pole",false);

  // This option creates elements that connect east to west at greenwich meridian
  // when true, instead of creating periodic ghost-points at east boundary when false
  options.set("three_dimensional",false);

  // This option sets number of parts the mesh will be split in
  options.set("nb_parts",1);

  // This option sets the part that will be generated
  options.set("part",0);

  // Experimental option. The result is a non-standard Reduced Gaussian Grid, with a ragged Greenwich line
  options.set("stagger", bool( Resource<bool>("-stagger;meshgen.stagger", false ) ) );

  // This option sets the maximum angle deviation for a quadrilateral element
  // max_angle = 30  -->  minimises number of triangles
  // max_angle = 0   -->  maximises number of triangles
  //Resource<double>( &Context::instance(), "meshgen.angle", 27 ).dump(Log::warning());
  //Log::warning() << "\n\n Getting max_angle_ ... " << std::endl;
  max_angle_ = Resource< double > ( "meshgen.angle", 29.5 );
  //Log::warning() << "max_angle_ = " << max_angle_ << std::endl;

  triangulate_quads_ = Resource< bool > ( "-triangulate;meshgen.triangulate", false );
  //Log::error() << "triangulate =" << triangulate_quads_ << std::endl;
}

Mesh* ReducedGridMeshGenerator::operator()( const ReducedGrid& grid )
{
  return generate(grid);
}


//std::vector<int> RGGMeshGenerator::partition(const RGG& rgg) const
//{
//  eckit::Log::info(Here())  << "partition start" << std::endl;
//  int nb_parts = options.get<int>("nb_parts");
//  EqualAreaPartitioner partitioner(nb_parts);
//  int ngptot = rgg.ngptot();
//  int n;
//  int p;
//  int i;
//  int begin;
//  int end;
//  int band;

//  double mem_per_node_in_MB = (32+2*32+32)/8./1024./1024.;
//  eckit::Log::info(Here())  << "required memory = " << static_cast<double>(ngptot)*mem_per_node_in_MB << " MB" << std::endl;
//  eckit::Log::info(Here())  << "required memory T8K ~ " << static_cast<double>(16e3*8e3*2./3.)*mem_per_node_in_MB << " MB" << std::endl;
//  eckit::Log::info(Here())  << "required memory T16K ~ " << static_cast<double>(32e3*16e3*2./3.)*mem_per_node_in_MB << " MB" << std::endl;

//  // Output
//  std::vector<int> part(ngptot);

//  // Create structure which we can sort with multiple keys (lat and lon)
//  std::vector<NodeInt> nodes(ngptot);
//  n=0;
//  for( int jlat=0; jlat<rgg.nlat(); ++jlat)
//  {
//    for( int jlon=0; jlon<rgg.nlon(jlat); ++jlon)
//    {
//      nodes[n].x = static_cast<int>(rgg.lon(jlat,jlon)*1e6);
//      nodes[n].y = static_cast<int>(rgg.lat(jlat)*1e6);
//      nodes[n].n = n;
//      ++n;
//    }
//  }
//  partitioner.partition(ngptot,nodes.data(),part.data());
//  return part;
//}

Mesh* ReducedGridMeshGenerator::generate(const ReducedGrid& rgg)
{
  int mypart   = options.get<int>("part");
  int nb_parts = options.get<int>("nb_parts");
  EqualAreaPartitioner partitioner(nb_parts);
  int n;
  int ngptot = rgg.npts();
  std::vector<int> part(ngptot);
  bool stagger = options.get<bool>("stagger");

  /*
  Create structure which we can partition with multiple keys (lat and lon)
  */
  std::vector<NodeInt> nodes(ngptot);
  n=0;
  for( int jlat=0; jlat<rgg.nlat(); ++jlat)
  {
    for( int jlon=0; jlon<rgg.nlon(jlat); ++jlon)
    {
      nodes[n].x = microdeg(rgg.lon(jlat,jlon));
      if( stagger ) nodes[n].x += static_cast<int>(1e6*180./static_cast<double>(rgg.nlon(jlat)));
      nodes[n].y = microdeg(rgg.lat(jlat));
      nodes[n].n = n;
      ++n;
    }
  }
  partitioner.partition(ngptot,nodes.data(),part.data());
  std::vector<NodeInt>().swap(nodes); // Deallocate completely

  Region region;
  generate_region(rgg,part,mypart,region);

  Mesh* mesh = generate_mesh(rgg,part,region);
  mesh->grid(rgg);
  return mesh;
}


void ReducedGridMeshGenerator::generate_region(const ReducedGrid& rgg, const std::vector<int>& parts, int mypart, Region& region)
{
  double max_angle = max_angle_;

  int n;
  /*
  Find min and max latitudes used by this part.
  */
  n=0;
  int lat_north=-1;
  for( int jlat=0; jlat<rgg.nlat(); ++jlat) {
    for( int jlon=0; jlon<rgg.nlon(jlat); ++jlon) {
      if( parts[n] == mypart ) {
        lat_north=jlat;
        goto end_north;
      }
      ++n;
    }
  } end_north:

  n=rgg.npts()-1;
  int lat_south=-1;
  for( int jlat=rgg.nlat()-1; jlat>=0; --jlat) {
    for( int jlon=rgg.nlon(jlat)-1; jlon>=0; --jlon) {
      if( parts[n] == mypart ) {
        lat_south=jlat;
        goto end_south;
      }
      --n;
    }
  } end_south:

  std::vector<int> offset(rgg.nlat(),0);

  n=0;
  for( int jlat=0; jlat<rgg.nlat(); ++jlat )
  {
    offset[jlat]=n;
    n+=rgg.nlon(jlat);
  };

  /*
  We need to connect to next region
  */
  int add_lat = 1;
  //if( rgg.lat(std::max(0,lat_north-add_lat)) > 0. )
    lat_north = std::max(0,lat_north-add_lat);
  //if( rgg.lat(std::min(rgg.nlat()-1,lat_south+add_lat)) < 0. )
    lat_south = std::min(rgg.nlat()-1,lat_south+add_lat);
  region.lat_begin.resize(rgg.nlat(),-1);
  region.lat_end.resize(rgg.nlat(),-1);
  region.nb_lat_elems.resize(rgg.nlat(),-1);
  region.north = lat_north;
  region.south = lat_south;

  ArrayShape shape = make_shape(region.south-region.north, 4*rgg.nlonmax(), 4);
  // eckit::Log::info(Here())  << "allocating elems" <<  "(" << extents[0] << "," << extents[1] << "," << extents[2] << ")" << std::endl;
  region.elems.resize(shape);
  region.elems = -1;

  int nelems=0;
  region.nquads=0;
  region.ntriags=0;

  ArrayView<int,3> elemview(region.elems);

  bool stagger = options.get<bool>("stagger");
  for (int jlat=lat_north; jlat<lat_south; ++jlat)
  {

    int ilat, latN, latS;
    int ipN1, ipN2, ipS1, ipS2;
    double xN1, xN2, yN, xS1, xS2, yS;
    double dN1S2, dS1N2, dN2S2;
    bool try_make_triangle_up, try_make_triangle_down, try_make_quad;
    bool add_triag, add_quad;

    ilat = jlat-region.north;

    ArrayView<int,2> lat_elems_view = elemview[ilat];

    latN = jlat;
    latS = jlat+1;
    yN = rgg.lat(latN);
    yS = rgg.lat(latS);

    int beginN, beginS, endN, endS;

    beginN = 0;
    endN   = rgg.nlon(latN); // include periodic point
    beginS = 0;
    endS   = rgg.nlon(latS); // include periodic point

    ipN1 = beginN;
    ipS1 = beginS;
    ipN2 = ipN1+1;
    ipS2 = ipS1+1;

    int jelem=0;

#if DEBUG_OUTPUT
    eckit::Log::info(Here())  << "=================" << std::endl;
#endif

    while (true)
    {
      if( ipN1 == endN && ipS1 == endS ) break;

      //ASSERT(offset[latN]+ipN1 < parts.size());
      //ASSERT(offset[latS]+ipS1 < parts.size());

      int pN1, pS1, pN2, pS2;
      if( ipN1 != rgg.nlon(latN) )
        pN1 = parts[offset[latN]+ipN1];
      else
        pN1 = parts[offset[latN]+ipN1-1];
      if( ipS1 != rgg.nlon(latS) )
        pS1 = parts[offset[latS]+ipS1];
      else
        pS1 = parts[offset[latS]+ipS1-1];

      if( ipN2 == rgg.nlon(latN) )
        pN2 = pN1;
      else
        pN2 = parts[offset[latN]+ipN2];
      if( ipS2 == rgg.nlon(latS) )
        pS2 = pS1;
      else
        pS2 = parts[offset[latS]+ipS2];

      //eckit::Log::info(Here())  << ipN1 << "("<<pN1<<") " << ipN2 <<"("<<pN2<<")" <<  std::endl;
      //eckit::Log::info(Here())  << ipS1 << "("<<pS2<<") " << ipS2 <<"("<<pS2<<")" <<  std::endl;


      xN1 = rgg.lon(latN,ipN1) * to_rad;
      xN2 = rgg.lon(latN,ipN2) * to_rad;
      xS1 = rgg.lon(latS,ipS1) * to_rad;
      xS2 = rgg.lon(latS,ipS2) * to_rad;

      if( stagger && (latN+1)%2==0 ) xN1 += M_PI/static_cast<double>(rgg.nlon(latN));
      if( stagger && (latN+1)%2==0 ) xN2 += M_PI/static_cast<double>(rgg.nlon(latN));
      if( stagger && (latS+1)%2==0 ) xS1 += M_PI/static_cast<double>(rgg.nlon(latS));
      if( stagger && (latS+1)%2==0 ) xS2 += M_PI/static_cast<double>(rgg.nlon(latS));

#if DEBUG_OUTPUT
      eckit::Log::info(Here())  << "-------" << std::endl;
#endif
      // eckit::Log::info(Here())  << "  access  " << region.elems.stride(0)*(jlat-region.north) + region.elems.stride(1)*jelem + 5 << std::endl;
//      eckit::Log::info(Here())  << ipN1 << "("<< xN1 << ")  " << ipN2 <<  "("<< xN2 << ")  " << std::endl;
//      eckit::Log::info(Here())  << ipS1 << "("<< xS1 << ")  " << ipS2 <<  "("<< xS2 << ")  " << std::endl;
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
        if( triangulate_quads_ )
        {
          if( false) //std::abs(alpha1) < 1 && std::abs(alpha2) < 1)
          {
            try_make_triangle_up   = (jlat+ipN1) % 2;
            try_make_triangle_down = (jlat+ipN1+1) % 2;
          }
          else
          {
            dN1S2 = std::abs(xN1-xS2);
            dS1N2 = std::abs(xS1-xN2);
            dN2S2 = std::abs(xN2-xS2);
            // eckit::Log::info(Here())  << "  dN1S2 " << dN1S2 << "   dS1N2 " << dS1N2 << "   dN2S2 " << dN2S2 << std::endl;
            if (dN1S2 <= dS1N2)
            {
              if (ipS1 != ipS2) { try_make_triangle_up = true; }
              else { try_make_triangle_down = true ; }
            }
            else if (dN1S2 >= dS1N2)
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
        dN2S2 = std::abs(xN2-xS2);
        // eckit::Log::info(Here())  << "  dN1S2 " << dN1S2 << "   dS1N2 " << dS1N2 << "   dN2S2 " << dN2S2 << std::endl;
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

      ArrayView<int,1> elem = lat_elems_view[jelem];

      if( try_make_quad )
      {
        // add quadrilateral
#if DEBUG_OUTPUT
        eckit::Log::info(Here())  << "          " << ipN1 << "  " << ipN2 << std::endl;
        eckit::Log::info(Here())  << "          " << ipS1 << "  " << ipS2 << std::endl;
#endif
        elem[0] = ipN1;
        elem[1] = ipS1;
        elem[2] = ipS2;
        elem[3] = ipN2;
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
          else if ( latS == rgg.nlat()-1 )
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

          if( region.lat_begin[latN] == -1 ) region.lat_begin[latN] = ipN1;
          if( region.lat_begin[latS] == -1 ) region.lat_begin[latS] = ipS1;
          region.lat_begin[latN] = std::min(region.lat_begin[latN], ipN1);
          region.lat_begin[latS] = std::min(region.lat_begin[latS], ipS1);
          region.lat_end[latN] = std::max(region.lat_end[latN], ipN2);
          region.lat_end[latS] = std::max(region.lat_end[latS], ipS2);
        }
        ipN1=ipN2;
        ipS1=ipS2;
      }
      else if( try_make_triangle_down  ) // make triangle down
      {
        // triangle without ip3
#if DEBUG_OUTPUT
        eckit::Log::info(Here())  << "          " << ipN1 << "  " << ipN2 << std::endl;
        eckit::Log::info(Here())  << "          " << ipS1 << std::endl;
#endif
        elem[0] = ipN1;
        elem[1] = ipS1;
        elem[2] = -1;
        elem[3] = ipN2;

        add_triag = false;

        int cnt_mypart = 0;
        int np[3] = {pN1, pN2, pS1};
        for( int j=0; j<3; ++j )
          if (np[j]==mypart) ++cnt_mypart;

        if( latN == 0 )
        {
          if( pN1 == mypart )
          {
            add_triag = true;
          }
        }
        else if ( latS == rgg.nlat()-1 )
        {
          if( pS1 == mypart )
          {
            add_triag = true;
          }
        }
        else if( cnt_mypart > 1 )
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

          if( region.lat_begin[latN] == -1 ) region.lat_begin[latN] = ipN1;
          if( region.lat_begin[latS] == -1 ) region.lat_begin[latS] = ipS1;
          region.lat_begin[latN] = std::min(region.lat_begin[latN], ipN1);
          region.lat_begin[latS] = std::min(region.lat_begin[latS], ipS1);
          region.lat_end[latN] = std::max(region.lat_end[latN], ipN2);
          region.lat_end[latS] = std::max(region.lat_end[latS], ipS1);
        }
        ipN1=ipN2;
        // and ipS1=ipS1;

      }
      else // make triangle up
      {
        // triangle without ip4
#if DEBUG_OUTPUT
        eckit::Log::info(Here())  << "          " << ipN1 << std::endl;
        eckit::Log::info(Here())  << "          " << ipS1 << "  " << ipS2 << std::endl;
#endif
        elem[0] = ipN1;
        elem[1] = ipS1;
        elem[2] = ipS2;
        elem[3] = -1;

        add_triag = false;


        int cnt_mypart = 0;
        int np[3] = {pN1, pS1, pS2};
        for( int j=0; j<3; ++j )
          if (np[j]==mypart) ++cnt_mypart;

        if( latN == 0 )
        {
          if( pN1 == mypart )
          {
            add_triag = true;
          }
        }
        else if ( latS == rgg.nlat()-1 )
        {
          if( pS2 == mypart )
          {
            add_triag = true;
          }
        }
        else if( cnt_mypart > 1 )
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
            if( latN == 0 || latS == rgg.nlat()-1 )
            {
              add_triag = true;
            }
            else if( 0.5*(yN+yS) > 1e-6 )
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

          if( region.lat_begin[latN] == -1 ) region.lat_begin[latN] = ipN1;
          if( region.lat_begin[latS] == -1 ) region.lat_begin[latS] = ipS1;
          region.lat_begin[latN] = std::min(region.lat_begin[latN], ipN1);
          region.lat_begin[latS] = std::min(region.lat_begin[latS], ipS1);
          region.lat_end[latN] = std::max(region.lat_end[latN], ipN1);
          region.lat_end[latS] = std::max(region.lat_end[latS], ipS2);
        }
        ipS1=ipS2;
        // and ipN1=ipN1;

      }
      ipN2 = std::min(endN,ipN1+1);
      ipS2 = std::min(endS,ipS1+1);
    }
    region.nb_lat_elems[jlat] = jelem;
    region.lat_end[latN] = std::min(region.lat_end[latN], rgg.nlon(latN)-1);
    region.lat_end[latS] = std::min(region.lat_end[latS], rgg.nlon(latS)-1);
  }

//  eckit::Log::info(Here())  << "nb_triags = " << region.ntriags << std::endl;
//  eckit::Log::info(Here())  << "nb_quads = " << region.nquads << std::endl;
//  eckit::Log::info(Here())  << "nb_elems = " << nelems << std::endl;

  int nb_region_nodes = 0;
  for( int jlat=lat_north; jlat<=lat_south; ++jlat ) {
    region.lat_begin[jlat] = std::max( 0, region.lat_begin[jlat] );
    nb_region_nodes += region.lat_end[jlat]-region.lat_begin[jlat]+1;

    // Count extra periodic node to be added in this case
    if( region.lat_end[jlat] == rgg.nlon(jlat)-1 ) ++nb_region_nodes;
  }
  region.nnodes = nb_region_nodes;
  if (region.nnodes == 0)
  {
    throw Exception("Trying to generate mesh with too many partitions. Reduce the number of partitions.",Here());
  }
}

Mesh* ReducedGridMeshGenerator::generate_mesh(const ReducedGrid& rgg,
                                      const std::vector<int>& parts,
                                      const Region& region)
{
  double tol = 1e-3;
  Mesh* mesh = new Mesh();
  int mypart = options.get<int>("part");
  int nparts = options.get<int>("nb_parts");
  int n, l;

  bool include_north_pole = (mypart == 0       ) && options.get<bool>("include_pole");
  bool include_south_pole = (mypart == nparts-1) && options.get<bool>("include_pole");
  bool three_dimensional  = options.get<bool>("three_dimensional");
  if( three_dimensional && nparts != 1 )
    throw BadParameter("Cannot generate three_dimensional mesh in parallel",Here());
  int nnodes  = region.nnodes;
  int ntriags = region.ntriags;
  int nquads  = region.nquads;
  if (include_north_pole) {
    ++nnodes;
    ntriags += rgg.nlon(0);
  }
  if (include_south_pole) {
    ++nnodes;
    ntriags += rgg.nlon(rgg.nlat()-1);
  }
  if (three_dimensional) {
    nnodes -= rgg.nlat();
  }

  std::vector<int> offset_glb(rgg.nlat());
  std::vector<int> offset_loc(region.south-region.north+1,0);

  n=0;
  for( int jlat=0; jlat<rgg.nlat(); ++jlat )
  {
    offset_glb[jlat]=n;
    n+=rgg.nlon(jlat);
  };

  std::vector<int> periodic_glb(rgg.nlat());

  if( !three_dimensional )
  {
    for( int jlat=0; jlat<rgg.nlat(); ++jlat )
    {
      periodic_glb[jlat] = n;
      ++n;
    }
  }
  else
  {
    for( int jlat=0; jlat<rgg.nlat(); ++jlat )
    {
      periodic_glb[jlat] = offset_glb[jlat] + rgg.nlon(jlat)-1;
    }
  }


  ArrayShape shape = make_shape(nnodes,Field::UNDEF_VARS);

  FunctionSpace& nodes = mesh->create_function_space( "nodes","LagrangeP1",shape );

  nodes.metadata().set("type",static_cast<int>(Entity::NODES));

  ArrayView<double,2> coords        ( nodes.create_field<double>("coordinates",   2) );
  ArrayView<gidx_t,1> glb_idx       ( nodes.create_field<gidx_t>("glb_idx",       1) );
  ArrayView<int,   1> part          ( nodes.create_field<int   >("partition",     1) );
  ArrayView<int,   1> flags         ( nodes.create_field<int   >("flags",         1) );

  bool stagger = options.get<bool>("stagger");

  int jnode=0;
  l=0;
  for( int jlat=region.north; jlat<=region.south; ++jlat )
  {
    int ilat=jlat-region.north;
    offset_loc[ilat]=l;
    l+=region.lat_end[jlat]-region.lat_begin[jlat]+1;
    double y = rgg.lat(jlat);

    for( int jlon=region.lat_begin[jlat]; jlon<=region.lat_end[jlat]; ++jlon )
    {
      n = offset_glb[jlat] + jlon;
      double x = rgg.lon(jlat,jlon);
      if( stagger && (jlat+1)%2==0 ) x += 180./static_cast<double>(rgg.nlon(jlat));

      coords(jnode,XX) = x;
      coords(jnode,YY) = y;
      glb_idx(jnode)   = n+1;
      part(jnode) = parts[n];
      Topology::reset(flags(jnode));
      if( jlat == 0 && !include_north_pole)
        Topology::set(flags(jnode),Topology::BC|Topology::NORTH);
      if( jlat == rgg.nlat()-1 && !include_south_pole)
        Topology::set(flags(jnode),Topology::BC|Topology::SOUTH);
      if( jlon == 0 && !three_dimensional)
        Topology::set(flags(jnode),Topology::BC|Topology::WEST);
      if( part(jnode) != mypart )
        Topology::set(flags(jnode),Topology::GHOST);

      //Log::info() << "meshgen " << std::setw(2) << glb_idx(jnode) << " ghost = " << Topology::check(flags(jnode),Topology::GHOST) << std::endl;
      ++jnode;
    }
    if( !three_dimensional &&  region.lat_end[jlat]==rgg.nlon(jlat)-1 ) // add periodic point
    {
      ++l;
      double x = rgg.lon(jlat,rgg.nlon(jlat));
      if( stagger && (jlat+1)%2==0 ) x += 180./static_cast<double>(rgg.nlon(jlat));

      coords(jnode,XX) = x;
      coords(jnode,YY) = y;
      glb_idx(jnode)   = periodic_glb[jlat]+1;
      part(jnode)      = part(jnode-1);
      Topology::reset(flags(jnode));
      Topology::set(flags(jnode),Topology::BC|Topology::EAST);
      if( part(jnode) != mypart )
        Topology::set(flags(jnode),Topology::GHOST);
      ++jnode;
    }
  };

  int jnorth=-1;
  if (include_north_pole)
  {
    jnorth = jnode;
    double y = 90.;
    double x = 180.;
    coords(jnode,XX) = x;
    coords(jnode,YY) = y;
    glb_idx(jnode)   = periodic_glb[rgg.nlat()-1]+2;
    part(jnode)      = mypart;
    Topology::reset(flags(jnode));
    Topology::set(flags(jnode),Topology::NORTH);
    ++jnode;
  }
  int jsouth=-1;
  if (include_south_pole)
  {
    jsouth = jnode;
    double y = -90.;
    double x =  180.;
    coords(jnode,XX) = x;
    coords(jnode,YY) = y;
    glb_idx(jnode)   = periodic_glb[rgg.nlat()-1]+3;
    part(jnode)      = mypart;
    Topology::reset(flags(jnode));
    Topology::set(flags(jnode),Topology::SOUTH);
    ++jnode;
  }

  shape = make_shape(nquads,Field::UNDEF_VARS);
  FunctionSpace& quads = mesh->create_function_space( "quads","LagrangeP1",shape );
  quads.metadata().set("type",static_cast<int>(Entity::ELEMS));
  IndexView<int,2> quad_nodes( quads.create_field<int>("nodes",4) );
  ArrayView<gidx_t,1> quad_glb_idx( quads.create_field<gidx_t>("glb_idx",1) );
  ArrayView<int,1> quad_part( quads.create_field<int>("partition",1) );

  shape = make_shape(ntriags,Field::UNDEF_VARS);
  FunctionSpace& triags = mesh->create_function_space( "triags","LagrangeP1",shape );
  triags.metadata().set("type",static_cast<int>(Entity::ELEMS));
  IndexView<int,2> triag_nodes( triags.create_field<int>("nodes",3) );
  ArrayView<gidx_t,1> triag_glb_idx( triags.create_field<gidx_t>("glb_idx",1) );
  ArrayView<int,1> triag_part( triags.create_field<int>("partition",1) );

  /*
  Fill in connectivity tables with global node indices first
  */
  int jquad = 0;
  int jtriag = 0;
  for (int jlat=region.north; jlat<region.south; ++jlat)
  {
    int ilat = jlat-region.north;
    int jlatN = jlat;
    int jlatS = jlat+1;
    int ilatN = ilat;
    int ilatS = ilat+1;
    for (int jelem=0; jelem<region.nb_lat_elems[jlat]; ++jelem)
    {
      const ArrayView<int,1> elem = ArrayView<int,3>(region.elems)[ilat][jelem];

      if(elem[2]>0 && elem[3]>0) // This is a quad
      {
        quad_nodes(jquad,0) = offset_loc[ilatN] + elem[0] - region.lat_begin[jlatN];
        quad_nodes(jquad,1) = offset_loc[ilatS] + elem[1] - region.lat_begin[jlatS];
        quad_nodes(jquad,2) = offset_loc[ilatS] + elem[2] - region.lat_begin[jlatS];
        quad_nodes(jquad,3) = offset_loc[ilatN] + elem[3] - region.lat_begin[jlatN];

        if( three_dimensional )
        {
          if (elem[2] == rgg.nlon(jlatS)) quad_nodes(jquad,2) = offset_loc[ilatS];
          if (elem[3] == rgg.nlon(jlatN)) quad_nodes(jquad,3) = offset_loc[ilatN];
        }

        // eckit::Log::info(Here())  << quad_nodes(0,jquad) << " " << quad_nodes(1,jquad) << " " << quad_nodes(2,jquad) << " " << quad_nodes(3,jquad) << std::endl;
        quad_glb_idx(jquad) = jquad+jtriag+1;
        quad_part(jquad) = mypart;
        ++jquad;
      }
      else // This is a triag
      {
        if(elem[3]<0) // This is a triangle pointing up
        {
          triag_nodes(jtriag,0) = offset_loc[ilatN] + elem[0] - region.lat_begin[jlatN];
          triag_nodes(jtriag,1) = offset_loc[ilatS] + elem[1] - region.lat_begin[jlatS];
          triag_nodes(jtriag,2) = offset_loc[ilatS] + elem[2] - region.lat_begin[jlatS];
          if( three_dimensional )
          {
            if (elem[0] == rgg.nlon(jlatN)) triag_nodes(jtriag,0) = offset_loc[ilatN];
            if (elem[2] == rgg.nlon(jlatS)) triag_nodes(jtriag,2) = offset_loc[ilatS];
          }
        }
        else // This is a triangle pointing down
        {
          triag_nodes(jtriag,0) = offset_loc[ilatN] + elem[0] - region.lat_begin[jlatN];
          triag_nodes(jtriag,1) = offset_loc[ilatS] + elem[1] - region.lat_begin[jlatS];
          triag_nodes(jtriag,2) = offset_loc[ilatN] + elem[3] - region.lat_begin[jlatN];
          if( three_dimensional )
          {
            if (elem[1] == rgg.nlon(jlatS)) triag_nodes(jtriag,1) = offset_loc[ilatS];
            if (elem[3] == rgg.nlon(jlatN)) triag_nodes(jtriag,2) = offset_loc[ilatN];
          }
        }
        triag_glb_idx(jtriag) = jquad+jtriag+1;
        triag_part(jtriag) = mypart;
        ++jtriag;
      }
    }
  }

  if (include_north_pole)
  {
    int jlat = 0;
    int ilat = 0;
    int ip1  = 0;
    for (int ip2=0; ip2<rgg.nlon(0); ++ip2)
    {
      int ip3 = ip2+1;
      if( three_dimensional && ip3 == rgg.nlon(0) ) ip3 = 0;
      triag_nodes(jtriag,0) = jnorth           + ip1;
      triag_nodes(jtriag,1) = offset_loc[ilat] + ip2;
      triag_nodes(jtriag,2) = offset_loc[ilat] + ip3;
      triag_glb_idx(jtriag) = jquad+jtriag+1;
      triag_part(jtriag) = mypart;
      ++jtriag;
    }
  }

  if (include_south_pole)
  {
    int jlat = rgg.nlat()-1;
    int ilat = region.south-region.north;
    int ip1 = 0;
    for (int ip2=1; ip2<rgg.nlon(jlat)+1; ++ip2)
    {
      int ip3 = ip2-1;
      triag_nodes(jtriag,0) = jsouth           + ip1;
      triag_nodes(jtriag,1) = offset_loc[ilat] + ip2;
      triag_nodes(jtriag,2) = offset_loc[ilat] + ip3;
      if( three_dimensional && ip2 == rgg.nlon(jlat) )
        triag_nodes(jtriag,1) = offset_loc[ilat] + 0;

      triag_glb_idx(jtriag) = jquad+jtriag+1;
      triag_part(jtriag) = mypart;
      ++jtriag;
    }
  }


  nodes.metadata().set("nb_owned",nnodes);
  quads.metadata().set("nb_owned",nquads);
  triags.metadata().set("nb_owned",ntriags);

  gidx_t max_glb_idx = rgg.npts()+rgg.nlat();
  if( three_dimensional ) max_glb_idx -= rgg.nlat();
  if( include_north_pole ) max_glb_idx += 1;
  if( include_south_pole ) max_glb_idx += 1;
  nodes.metadata().set("max_glb_idx",max_glb_idx);
  quads.metadata().set("max_glb_idx", nquads+ntriags);
  triags.metadata().set("max_glb_idx",nquads+ntriags);

  generate_global_element_numbering( *mesh );

  return mesh;
}

void ReducedGridMeshGenerator::generate_global_element_numbering( Mesh& mesh )
{
  int loc_nb_elems = 0;
  std::vector<int> elem_counts( MPL::size() );
  std::vector<int> elem_displs( MPL::size() );
  for( int f=0; f<mesh.nb_function_spaces(); ++f )
  {
    FunctionSpace& elements = mesh.function_space(f);
    if( elements.metadata().get<int>("type") == Entity::ELEMS )
    {
      loc_nb_elems += elements.shape(0);
    }
  }
  MPL_CHECK_RESULT( MPI_Allgather( &loc_nb_elems, 1, MPI_INT,
                                   elem_counts.data(), 1, MPI_INT, MPI_COMM_WORLD) );
  elem_displs[0] = 0;
  for( int jpart=1; jpart<MPL::size(); ++jpart )
  {
    elem_displs[jpart] = elem_displs[jpart-1] + elem_counts[jpart-1];
  }

  gidx_t gid = 1+elem_displs[ MPL::rank() ];

  for( int f=0; f<mesh.nb_function_spaces(); ++f )
  {
    FunctionSpace& elements = mesh.function_space(f);
    if( elements.metadata().get<int>("type") == Entity::ELEMS )
    {
      ArrayView<gidx_t,1> glb_idx( elements.field("glb_idx") );
      int nb_elems = elements.shape(0);
      for( int e=0; e<nb_elems; ++e )
      {
        glb_idx(e) = gid++;
      }
    }
  }
}

} // namespace meshgen
} // namespace atlas

