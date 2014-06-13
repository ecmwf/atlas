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

#include "atlas/mesh/Array.hpp"
#include "atlas/mesh/ArrayView.hpp"
#include "atlas/mesh/IndexView.hpp"
#include "atlas/mesh/Field.hpp"
#include "atlas/mesh/FunctionSpace.hpp"
#include "atlas/mesh/Mesh.hpp"
#include "atlas/mesh/Util.hpp"
#include "atlas/mesh/Parameters.hpp"
#include "atlas/meshgen/EqualAreaPartitioner.hpp"
#include "atlas/meshgen/RGG.hpp"

#define DEBUG_OUTPUT 0
namespace atlas {
namespace meshgen {

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

double spherical_distance(double x1,double y1,double x2,double y2,double rad)
{
  using namespace std;
  if((std::abs(x2-x1))<1e-8) return (y2-y1)*rad;
  else if((std::abs(y2-y1))<1e-8) return (x2-x1)*rad;
  else return acos(cos(x1)*cos(y1)*cos(x2)*cos(y2)+cos(x1)*sin(y1)*cos(x2)*sin(y2)+sin(x1)*sin(x2))*rad;
}

int RGG::ngptot() const
{
  return std::accumulate(lon_.data(),lon_.data()+lon_.size(),0);
}  

RGGMeshGenerator::RGGMeshGenerator()
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
  
}

std::vector<int> RGGMeshGenerator::partition(const RGG& rgg) const
{
  std::cout << "partition start" << std::endl;
  int nb_parts = options.get<int>("nb_parts");
  EqualAreaPartitioner partitioner(nb_parts);
  int ngptot = rgg.ngptot();
  int n;
  int p;
  int i;
  int begin;
  int end;
  int band;
  
  double mem_per_node_in_MB = (32+2*32+32)/8./1024./1024.;
  std::cout << "required memory = " << static_cast<double>(ngptot)*mem_per_node_in_MB << " MB" << std::endl;
  std::cout << "required memory T8K ~ " << static_cast<double>(16e3*8e3*2./3.)*mem_per_node_in_MB << " MB" << std::endl;
  std::cout << "required memory T16K ~ " << static_cast<double>(32e3*16e3*2./3.)*mem_per_node_in_MB << " MB" << std::endl;
  
  // Output
  std::vector<int> part(ngptot);

  // Create structure which we can sort with multiple keys (lat and lon)
  std::vector<NodeInt> nodes(ngptot);
  n=0;
  for( int jlat=0; jlat<rgg.nlat(); ++jlat)
  {
    for( int jlon=0; jlon<rgg.nlon(jlat); ++jlon)
    {
      nodes[n].x = static_cast<int>(rgg.lon(jlon,jlat)*1e6);
      nodes[n].y = static_cast<int>(rgg.lat(jlat)*1e6);
      nodes[n].n = n;
      ++n;
    }
  }
  partitioner.partition(ngptot,nodes.data(),part.data());
  return part;
}

Mesh* RGGMeshGenerator::generate(const RGG& rgg)
{
  int mypart   = options.get<int>("part");
  int nb_parts = options.get<int>("nb_parts");
  EqualAreaPartitioner partitioner(nb_parts);
  int ngptot = rgg.ngptot();
  int n;
  
  std::vector<int> part(ngptot);

  /*
  Create structure which we can partition with multiple keys (lat and lon)
  */
  std::vector<NodeInt> nodes(ngptot);
  n=0;
  for( int jlat=0; jlat<rgg.nlat(); ++jlat)
  {
    for( int jlon=0; jlon<rgg.nlon(jlat); ++jlon)
    {
      nodes[n].x = static_cast<int>(rgg.lon(jlon,jlat)*1e6);
      nodes[n].y = static_cast<int>(rgg.lat(jlat)*1e6);
      nodes[n].n = n;
      ++n;
    }
  }
  partitioner.partition(ngptot,nodes.data(),part.data());
  std::vector<NodeInt>().swap(nodes); // Deallocate completely
  
  Region region;
  generate_region(rgg,part,mypart,region);
  
  Mesh* mesh = generate_mesh(rgg,part,region);
  return mesh;
}


void RGGMeshGenerator::generate_region(const RGG& rgg, const std::vector<int>& parts, int mypart, Region& region)
{
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
  
  n=rgg.ngptot()-1;
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
  if( rgg.lat(std::max(0,lat_north-1)) > 0. ) 
    lat_north = std::max(0,lat_north-1);
  if( rgg.lat(std::min(rgg.nlat()-1,lat_south+1)) < 0. ) 
    lat_south = std::min(rgg.nlat()-1,lat_south+1);
  region.lat_begin.resize(rgg.nlat(),-1);
  region.lat_end.resize(rgg.nlat(),-1);
  region.nb_lat_elems.resize(rgg.nlat(),-1);
  region.north = lat_north;
  region.south = lat_south;
  
  std::vector<int> extents = Extents(region.south-region.north, 2*rgg.nlonmax(), 4);
  // std::cout << "allocating elems" <<  "(" << extents[0] << "," << extents[1] << "," << extents[2] << ")" << std::endl;
  region.elems.resize(extents);
  region.elems = -1;
  
  int nelems=0;
  region.nquads=0;
  region.ntriags=0;

  ArrayView<int,3> elemview(region.elems);
  
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
    std::cout << "=================" << std::endl;
#endif
    
    while (true)
    {
      int pN1 = parts[offset[latN]+ipN1];
      int pS1 = parts[offset[latS]+ipS1];
      int pN2;
      int pS2;
      if( ipN2 == rgg.nlon(latN) )
        pN2 = pN1;
      else
        pN2 = parts[offset[latN]+ipN2];
      if( ipS2 == rgg.nlon(latS) )
        pS2 = pS1;
      else
        pS2 = parts[offset[latS]+ipS2];

      if( ipN1 == endN && ipS1 == endS ) break;
      
      xN1 = rgg.lon(ipN1,latN);
      xN2 = rgg.lon(ipN2,latN);
      xS1 = rgg.lon(ipS1,latS);
      xS2 = rgg.lon(ipS2,latS);
  
#if DEBUG_OUTPUT
      std::cout << "-------" << std::endl;
#endif
      // std::cout << "  access  " << region.elems.stride(0)*(jlat-region.north) + region.elems.stride(1)*jelem + 5 << std::endl;
      // std::cout << ipN1 << "("<< xN1 << ")  " << ipN2 <<  "("<< xN2 << ")  " << std::endl;
      // std::cout << ipS1 << "("<< xS1 << ")  " << ipS2 <<  "("<< xS2 << ")  " << std::endl;
      try_make_triangle_up   = false;
      try_make_triangle_down = false;
      try_make_quad = false;
    
      
      dN1S2 = std::abs(xN1-xS2);//std::abs(spherical_distance(xN1,yN,xS2,yS,1.));
      dS1N2 = std::abs(xS1-xN2);//std::abs(spherical_distance(xS1,yS,xN2,yN,1.));
      dN2S2 = std::abs(xN2-xS2);//std::abs(spherical_distance(xN2,yN,xS2,yS,1.));
      // std::cout << "  dN1S2 " << dN1S2 << "   dS1N2 " << dS1N2 << "   dN2S2 " << dN2S2 << std::endl;
      if ( (dN1S2 < dN2S2 && dN1S2 < dS1N2) && (ipS1 != ipS2) ) try_make_triangle_up = true;
      else if ( (dS1N2 < dN2S2 && dS1N2 < dN1S2) && (ipN1 != ipN2) ) try_make_triangle_down = true;
      else try_make_quad = true;
            
      // BlockT<int,2> lat_elems(region.elems, Idx(ilat));
      
      ArrayView<int,1> elem = lat_elems_view[jelem];
      // int* elem = lat_elems[jelem];
      // int* elem = view(region.elems,ilat,jelem);
      
      if( try_make_quad )
      {
        // add quadrilateral
#if DEBUG_OUTPUT
        std::cout << "          " << ipN1 << "  " << ipN2 << std::endl;
        std::cout << "          " << ipS1 << "  " << ipS2 << std::endl;
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

          if( cnt_mypart > 2 ) // 3 or 4 points belong to mypart
          {
            add_quad = true;
          }
          else if( cnt_max < 3 ) // 3 or 4 points don't belong to mypart
          {
            if( 0.5*(yN+yS) > 1e-6 )
            {
              if ( pS2 == mypart )  add_quad = true;
            }
            else
            {
              if ( pN1 == mypart )  add_quad = true;
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
        std::cout << "          " << ipN1 << "  " << ipN2 << std::endl;
        std::cout << "          " << ipS1 << std::endl;
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
        
          if( region.lat_begin[latN] == -1 ) region.lat_begin[latN] = ipN1;
          if( region.lat_begin[latS] == -1 ) region.lat_begin[latS] = ipS1;
          region.lat_begin[latN] = std::min(region.lat_begin[latN], ipN1);
          region.lat_begin[latS] = std::min(region.lat_begin[latS], ipS1);
          region.lat_end[latN] = std::max(region.lat_end[latN], ipN2);
          region.lat_end[latS] = std::max(region.lat_end[latS], ipS1);
        }
        ipN1=ipN2;
        ipS1=ipS1;
      }
      else // make triangle up
      {
        // triangle without ip4
#if DEBUG_OUTPUT
        std::cout << "          " << ipN1 << std::endl;
        std::cout << "          " << ipS1 << "  " << ipS2 << std::endl;
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
        ipN1=ipN1;
        ipS1=ipS2;
      }
      ipN2 = ipN1+1;
      ipS2 = ipS1+1;
    }
    region.nb_lat_elems[jlat] = jelem;
    region.lat_end[latN] = std::min(region.lat_end[latN], rgg.nlon(latN)-1);
    region.lat_end[latS] = std::min(region.lat_end[latS], rgg.nlon(latS)-1);
  }
  
  int nb_region_nodes = 0;
  for( int jlat=lat_north; jlat<=lat_south; ++jlat ) {
    region.lat_begin[jlat] = std::max( 0, region.lat_begin[jlat] );
    nb_region_nodes += region.lat_end[jlat]-region.lat_begin[jlat]+1;
    
    // Count extra periodic node to be added in this case
    if( region.lat_end[jlat] == rgg.nlon(jlat)-1 ) ++nb_region_nodes;
  }
  region.nnodes = nb_region_nodes;
}

Mesh* RGGMeshGenerator::generate_mesh(const RGG& rgg,
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
  bool three_dimensional  = (nparts == 1       ) && options.get<bool>("three_dimensional");
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

  
  std::vector<int> extents = Extents(nnodes,Field::UNDEF_VARS);
  FunctionSpace& nodes =
      mesh->add_function_space( new FunctionSpace("nodes","LagrangeP1",extents) );
  nodes.metadata().set("type",static_cast<int>(Entity::NODES));
  ArrayView<double,2> coords        ( nodes.create_field<double>("coordinates",   2) );
  ArrayView<int,   1> glb_idx       ( nodes.create_field<int   >("glb_idx",       1) );
  ArrayView<int,   1> part          ( nodes.create_field<int   >("partition",     1) );
  IndexView<int,   1> ridx          ( nodes.create_field<int   >("remote_idx",    1) );
  
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
      double x = rgg.lon(jlon,jlat);
      coords(jnode,XX) = x;
      coords(jnode,YY) = y;
      glb_idx(jnode)   = n+1;
      part(jnode) = parts[n];
      if( part(jnode) == mypart ) ridx(jnode) = jnode;
      else ridx(jnode) = -1;
      ++jnode;
    }
    if( !three_dimensional &&  region.lat_end[jlat]==rgg.nlon(jlat)-1 ) // add periodic point
    {
      ++l;
      double x = rgg.lon(rgg.nlon(jlat),jlat);
      coords(jnode,XX) = x;
      coords(jnode,YY) = y;
      glb_idx(jnode)   = periodic_glb[jlat]+1;
      part(jnode)      = parts[ offset_glb[jlat] ];
      if( part(jnode) == mypart ) ridx(jnode) = offset_glb[jlat];
      else ridx(jnode) = -1;
      ++jnode;
    }
  };
  
  int jnorth=-1;
  if (include_north_pole)
  {
    jnorth = jnode;
    double y = M_PI_2;
    double x = M_PI;
    coords(jnode,XX) = x;
    coords(jnode,YY) = y;
    glb_idx(jnode)   = periodic_glb[rgg.nlat()-1]+2;
    part(jnode)      = mypart;
    ridx(jnode)      = jnode;
    ++jnode;
  }
  int jsouth=-1;
  if (include_south_pole)
  {
    jsouth = jnode;
    double y = -M_PI_2;
    double x =  M_PI;
    coords(jnode,XX) = x;
    coords(jnode,YY) = y;
    glb_idx(jnode)   = periodic_glb[rgg.nlat()-1]+3;
    part(jnode)      = mypart;
    ridx(jnode)      = jnode;
    ++jnode;
  }
    
  extents = Extents(nquads,Field::UNDEF_VARS);
  FunctionSpace& quads = mesh->add_function_space( new FunctionSpace("quads","LagrangeP1",extents) );
  quads.metadata().set("type",static_cast<int>(Entity::ELEMS));
  IndexView<int,2> quad_nodes( quads.create_field<int>("nodes",4) );
  ArrayView<int,1> quad_glb_idx( quads.create_field<int>("glb_idx",1) );
  ArrayView<int,1> quad_part( quads.create_field<int>("partition",1) );
  
  extents = Extents(ntriags,Field::UNDEF_VARS);
  FunctionSpace& triags = mesh->add_function_space( new FunctionSpace("triags","LagrangeP1",extents) );
  triags.metadata().set("type",static_cast<int>(Entity::ELEMS));
  IndexView<int,2> triag_nodes( triags.create_field<int>("nodes",3) );
  ArrayView<int,1> triag_glb_idx( triags.create_field<int>("glb_idx",1) );
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
        
        // std::cout << quad_nodes(0,jquad) << " " << quad_nodes(1,jquad) << " " << quad_nodes(2,jquad) << " " << quad_nodes(3,jquad) << std::endl;
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
  
  int max_glb_idx = rgg.ngptot()+rgg.nlat();
  if( three_dimensional ) max_glb_idx -= rgg.nlat();
  if( include_north_pole ) max_glb_idx += 1;
  if( include_south_pole ) max_glb_idx += 1;
  nodes.metadata().set("max_glb_idx",max_glb_idx);
  quads.metadata().set("max_glb_idx", nquads+ntriags);
  triags.metadata().set("max_glb_idx",nquads+ntriags);

  generate_global_element_numbering( *mesh );

  return mesh;
}

void RGGMeshGenerator::generate_global_element_numbering( Mesh& mesh )
{
  int loc_nb_elems = 0;
  std::vector<int> elem_counts( MPL::size() );
  std::vector<int> elem_displs( MPL::size() );
  for( int f=0; f<mesh.nb_function_spaces(); ++f )
  {
    FunctionSpace& elements = mesh.function_space(f);
    if( elements.metadata().get<int>("type") == Entity::ELEMS )
    {
      loc_nb_elems += elements.extents()[0];
    }
  }
  MPL_CHECK_RESULT( MPI_Allgather( &loc_nb_elems, 1, MPI_INT,
                                   elem_counts.data(), 1, MPI_INT, MPI_COMM_WORLD) );
  elem_displs[0] = 0;
  for( int jpart=1; jpart<MPL::size(); ++jpart )
  {
    elem_displs[jpart] = elem_displs[jpart-1] + elem_counts[jpart-1];
  }

  int gid = 1+elem_displs[ MPL::rank() ];

  for( int f=0; f<mesh.nb_function_spaces(); ++f )
  {
    FunctionSpace& elements = mesh.function_space(f);
    if( elements.metadata().get<int>("type") == Entity::ELEMS )
    {
      ArrayView<int,1> glb_idx( elements.field("glb_idx") );
      int nb_elems = elements.extents()[0];
      for( int e=0; e<nb_elems; ++e )
      {
        glb_idx(e) = gid++;
      }
    }
  }
}

} // namespace meshgen
} // namespace atlas

