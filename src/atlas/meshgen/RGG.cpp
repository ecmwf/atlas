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
#include <cmath>

#include "atlas/Array.hpp"
#include "atlas/ArrayView.hpp"
#include "atlas/Field.hpp"
#include "atlas/FunctionSpace.hpp"
#include "atlas/Mesh.hpp"
#include "atlas/Parameters.hpp"
#include "atlas/meshgen/EqualAreaPartitioner.hpp"
#include "atlas/meshgen/RGG.hpp"

namespace atlas {
namespace meshgen {

  DebugMesh::DebugMesh()
  {
    int nlat=5;
    int lon[] = {
      6,
      12,
      16,
      18,
      18,
    };
    /*
    First prediction of colatitudes
    */
    std::vector<double> colat(nlat);
    double z;
    for( int i=0; i<nlat; ++i )
    {
      z = (4.*(i+1.)-1.)*M_PI/(4.*2.*nlat+2.);
      colat[i] = z+1./(tan(z)*(8.*(2.*nlat)*(2.*nlat)));
    }
    /*
    Fill in final structures
    */
    lat_.resize(2*nlat);
    lon_.resize(2*nlat);
    std::copy( lon, lon+nlat, lon_.begin() );
    std::reverse_copy( lon, lon+nlat, lon_.begin()+nlat );
    std::copy( colat.begin(), colat.begin()+nlat, lat_.begin() );
    std::reverse_copy( colat.begin(), colat.begin()+nlat, lat_.begin()+nlat );
    for (int i=0; i<nlat; ++i)
      lat_[i]=M_PI/2.-lat_[i];
    for (int i=nlat; i<2*nlat; ++i)
      lat_[i]=-M_PI/2.+lat_[i];
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
  options.set("include_pole",false);
  options.set("distributed",true);
  options.set("nb_parts",1);
  options.set("part",0);
  options.set("halo",1);
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
  
  std::cout << "latitudes = " << lat_north << " , " << lat_south << std::endl;
  
  std::vector<int> offset(rgg.nlat(),0);

  n=0;
  for( int jlat=0; jlat<rgg.nlat(); ++jlat )
  {
    offset[jlat]=n;
    n+=rgg.nlon(jlat);
  };

  /*
  Find begin and end points on latitudes
  Also count number of nodes belonging to this part
  */
  std::vector<int> lat_nodes_begin(rgg.nlat(),-1);
  std::vector<int> lat_nodes_end(rgg.nlat(),-1);
  for( int jlat=lat_north; jlat<=lat_south; ++jlat ) {
    for( int jlon=0; jlon<rgg.nlon(jlat); ++jlon) {
      if( lat_nodes_begin[jlat]<0 && parts[offset[jlat]+jlon] == mypart)
      {
        lat_nodes_begin[jlat]=jlon;
        break;
      }
    }
    for( int jlon=rgg.nlon(jlat)-1; jlon>=0; --jlon) {
      if( lat_nodes_end[jlat]<0 && parts[offset[jlat]+jlon] == mypart)
      {
        lat_nodes_end[jlat]=jlon;
        break;
      }
    }
    std::cout << "lat["<<jlat<<"] :  [" << lat_nodes_begin[jlat] << "("<< rgg.lon(lat_nodes_begin[jlat],jlat) << ") -- " << lat_nodes_end[jlat] << "("<< rgg.lon(lat_nodes_end[jlat],jlat) << ")]" << std::endl;
  }
  
  /*
  We need to connect to next region
  */
  for( int jlat=lat_north; jlat<=lat_south; ++jlat ) {
    lat_nodes_end[jlat] = std::min(rgg.nlon(jlat)-1, lat_nodes_end[jlat]+1);
  }
  if( (lat_north+lat_south)/2 <= rgg.nlat()/2 ) // This region is mostly on north hemisphere
  {
    if( lat_north > 0 ) 
    {
      int connect_lat = lat_north-1;
      double yN = rgg.lat(connect_lat);
      double yS = rgg.lat(lat_north);
      double xS;
      double xN;
      double dist;
      xS = rgg.lon(lat_nodes_begin[lat_north],lat_north);
      dist=2.*M_PI;
      for (int jlon=0; jlon<rgg.nlon(connect_lat); ++jlon )
      {
        xN = rgg.lon(jlon,connect_lat);
        double dist_new = std::abs(xN-xS);
        if (dist_new <= dist)
        {
          lat_nodes_begin[connect_lat] = jlon;
          dist = dist_new;
        }
        else
        {
          break;
        }
      }
      xS = rgg.lon(lat_nodes_end[lat_north],lat_north);
      dist=2.*M_PI;
      for (int jlon=rgg.nlon(connect_lat)-1; jlon>=0; --jlon )
      {
        xN = rgg.lon(jlon,connect_lat);
        double dist_new = std::abs(xN-xS);
        if (dist_new <= dist)
        {
          lat_nodes_end[connect_lat] = jlon;
          dist = dist_new;
        }
        else
        {
          break;
        }
      }
      lat_north = connect_lat;
      std::cout << "now add lat["<<lat_north<<"] : " << lat_nodes_begin[lat_north] << "  to " << lat_nodes_end[lat_north] << std::endl;
    }
  }
  
  for( int jlat=lat_north; jlat<=lat_south; ++jlat )
  {
    // std::cout << "lat["<<jlat<<"] :  [" << lat_nodes_begin[jlat] << "("<< rgg.lon(lat_nodes_begin[jlat],jlat) << ") -- " << lat_nodes_end[jlat] << "("<< rgg.lon(lat_nodes_end[jlat],jlat) << ")]" << std::endl;
  }
  
  region.lat_begin.resize(rgg.nlat(),-1);
  region.lat_end.resize(rgg.nlat(),-1);
  region.nb_lat_elems.resize(rgg.nlat(),-1);
  region.north = lat_north;
  region.south = lat_south;
  for( int jlat=region.north; jlat<=region.south; ++jlat )
  {
    region.lat_begin[jlat] = lat_nodes_begin[jlat];
    region.lat_end[jlat]   = lat_nodes_end[jlat];
  }
  
  std::vector<int> extents = Extents(region.south-region.north, 2*rgg.nlonmax(), 4);
  // std::cout << "allocating elems" <<  "(" << extents[0] << "," << extents[1] << "," << extents[2] << ")" << std::endl;
  region.elems.resize(extents);
  region.elems = -1;
  
  int nelems=0;
  region.nquads=0;
  region.ntriags=0;

  ArrayView<int,3> elemview(region.elems);
  
  std::cout << "generating..." << std::endl;
  for (int jlat=region.north; jlat<region.south; ++jlat)
  {
    
    int ilat, latN, latS;
    int ipN1, ipN2, ipS1, ipS2;
    double xN1, xN2, yN, xS1, xS2, yS;
    double dist, dist_new;
    double dN1S2, dS1N2, dN2S2;
    bool make_triangle_up, make_triangle_down, make_quad;
    
    ilat = jlat-region.north;
    
    ArrayView<int,2> lat_elems_view = elemview[ilat];
    
    latN = jlat;
    latS = jlat+1;
    yN = rgg.lat(latN);
    yS = rgg.lat(latS);

    int beginN, beginS, endN, endS;

    // Initial points
    {
      beginN = region.lat_begin[latN];
      xN1 = rgg.lon(beginN,latN);
      dist=2.*M_PI;
      for (int jlon=std::max(0,region.lat_begin[latS]-1); jlon<=region.lat_end[latS]; ++jlon )
      {
        xS1 = rgg.lon(jlon,latS);
        dist_new = std::abs(xN1-xS1);
        if (dist_new - 1e-6 < dist)
        {
          dist = dist_new;
          beginS = jlon;
        }
        else
        {
          // beginN connects beginS
          // But does it belong to this partition?
          // 1) if both beginN and beginS belong to mypart, then YES
          //    By construction beginN is already belonging to mypart.
          // 2) if beginS is not mypart, then increase it.
          //    Now we have to check if beginN needs to be increased, by
          //    checking spherical distances
          int partS = parts[offset[latS]+beginS];
          if( partS != mypart )
          {
            ++beginS;
            xS1 = rgg.lon(beginS,latS);
            double xN1_next = rgg.lon(beginN+1,latN);
            if( std::abs(xS1-xN1_next) < std::abs(xS1-xN1) ) ++beginN;
          } 
          break;
        }
      }
    }

    // End points
    {
      endS = region.lat_end[latS];
      std::cout << "endS = " << endS << std::endl;
      xS2 = rgg.lon(endS,latS);
      dist=2.*M_PI;
      for (int jlon=std::max(0,region.lat_begin[latN]-1); jlon<=region.lat_end[latN]+2; ++jlon )
      {
        xN2 = rgg.lon(jlon,latN);
        dist_new = std::abs(xS2-xN2);
        if (dist_new < dist)
        {
          std::cout << "dist = " << dist << std::endl;
          dist = dist_new;
          endN = jlon;
        }
        else // minimum distance found
        {
          std::cout << "larger dist = " << dist_new << std::endl;
          // ipN2 connects ipS2
          // But does it connect to next partition?
          // 1) if both ipN2 and ipS2 belong to other parts than mypart, then YES
          //    By construction ipS2 is already belonging to other part.
          // 2) if ipN2 is still mypart, then there are 2 options...
          //    a) if this is a latitude connecting to top, don't increase ipN2
          //    b) else, if it has to connect to side, increase it.
          //       Now we have to check if ipS2 needs to be increased, by
          //       checking spherical distances
          
          if( ilat == 0 ) break;
          
          int partN = parts[offset[latN]+endN];
          std::cout << "partN = " << partN << std::endl;
          if( partN == mypart )
          {
            // if( partN_next != mypart && ilat!=0)
            {
              std::cout << "dist = " << dist << "    endN : " << endN;
              endN = std::min(endN+1, rgg.nlon(latN)-1);
              std::cout << " --> " << endN << std::endl;
              xN2 = rgg.lon(endN,latN);
              double xS2_next = rgg.lon(endS+1,latS);
              if( std::abs(xN2-xS2_next) < std::abs(xN2-xS2) ){
                std::cout << "endS:  " << endS;
                endS = std::min(endS+1, rgg.nlon(latS)-1);    
                std::cout << " --> " << endS << std::endl;                        
              }
            }
          } 
          break;
        }
      }
    }

    std::cout << "begin strip["<<jlat<<"] : connect N"<< beginN << " with S" << beginS << std::endl;
    std::cout << "end strip  ["<<jlat<<"] : connect N"<< endN << " with S" << endS << std::endl;

    ipN1 = beginN;
    ipS1 = beginS;
    ipN2 = ipN1+1;
    ipS2 = ipS1+1;

    int jelem=0;
    std::cout << "=================" << std::endl;
    //continue until ipN1 and ipS1 are no longer of my part, or the east domain boundary
    while (true)
    {
      region.lat_begin[latN] = std::min(ipN1,region.lat_begin[latN]);
      region.lat_begin[latS] = std::min(ipS1,region.lat_begin[latS]);
      
      region.lat_end[latN] = std::max(ipN1,region.lat_end[latN]);
      region.lat_end[latS] = std::max(ipS1,region.lat_end[latS]);
      // int pN1 = parts[offset[latN]+ipN1];
      // int pS1 = parts[offset[latS]+ipS1];
      // if ( (ipN1>mypart && pS1>mypart) || (ipN1==rgg.nlon(latN)-1&&ipS1==rgg.nlon(latS)-1) )
        // break;
      if( ipN1 == endN && ipS1 == endS ) break;
      
      xN1 = rgg.lon(ipN1,latN);
      xN2 = rgg.lon(ipN2,latN);
      xS1 = rgg.lon(ipS1,latS);
      xS2 = rgg.lon(ipS2,latS);
  
      std::cout << "-------" << std::endl;
      // std::cout << "  access  " << region.elems.stride(0)*(jlat-region.north) + region.elems.stride(1)*jelem + 5 << std::endl;
      // std::cout << ipN1 << "("<< xN1 << ")  " << ipN2 <<  "("<< xN2 << ")  " << std::endl;
      // std::cout << ipS1 << "("<< xS1 << ")  " << ipS2 <<  "("<< xS2 << ")  " << std::endl;
      make_triangle_up   = false;
      make_triangle_down = false;
      make_quad = false;
    
      
      dN1S2 = std::abs(spherical_distance(xN1,yN,xS2,yS,1.));
      dS1N2 = std::abs(spherical_distance(xS1,yS,xN2,yN,1.));
      dN2S2 = std::abs(spherical_distance(xN2,yN,xS2,yS,1.));
      // std::cout << "  d31 " << d31 << "   d42 " << d42 << "   d43 " << d43 << std::endl;
      if ( (dN1S2 < dN2S2 && dN1S2 < dS1N2) && (ipS1 != ipS2) ) make_triangle_up = true;
      else if ( (dS1N2 < dN2S2 && dS1N2 < dN1S2) && (ipN1 != ipN2) ) make_triangle_down = true;
      else make_quad = true;
            
      // BlockT<int,2> lat_elems(region.elems, Idx(ilat));
      
      ArrayView<int,1> elem = lat_elems_view[jelem];
      // int* elem = lat_elems[jelem];
      // int* elem = view(region.elems,ilat,jelem);
      
      if( make_quad )
      {
        // add quadrilateral
        std::cout << "          " << ipN1 << "  " << ipN2 << std::endl;
        std::cout << "          " << ipS1 << "  " << ipS2 << std::endl;
        elem[0] = ipN1;
        elem[1] = ipS1;
        elem[2] = ipS2;
        elem[3] = ipN2;
        ipN1=std::min(rgg.nlon(latN)-1,ipN2);
        ipS1=std::min(rgg.nlon(latS)-1,ipS2);
        // std::cout << "  .  " << std::endl;
        // if( (pN1>=mypart && pS1>=mypart) )
          ++region.nquads;
      }
      else if( make_triangle_down  ) // make triangle down
      {
        // triangle without ip3
        std::cout << "          " << ipN1 << "  " << ipN2 << std::endl;
        std::cout << "          " << ipS1 << std::endl;
        elem[0] = ipN1;
        elem[1] = ipS1;
        elem[2] = -1;
        elem[3] = ipN2;
        ipN1=std::min(rgg.nlon(latN)-1,ipN2);
        ipS1=std::min(rgg.nlon(latS)-1,ipS1);
        // std::cout << "  ." << std::endl;
        // if( (pN1>=mypart && pS1>=mypart) )
          ++region.ntriags;
      }
      else // make triangle up
      {
        // triangle without ip4
        std::cout << "          " << ipN1 << std::endl;
        std::cout << "          " << ipS1 << "  " << ipS2 << std::endl;
        elem[0] = ipN1;
        elem[1] = ipS1;
        elem[2] = ipS2;
        elem[3] = -1;
        ipN1=std::min(rgg.nlon(latN)-1,ipN1);
        ipS1=std::min(rgg.nlon(latS)-1,ipS2);
        // std::cout << "  . "  << std::endl;
        // if( (pN1>=mypart && pS1>=mypart) )
          ++region.ntriags;
      }
      ipN2 = ipN1+1;
      ipS2 = ipS1+1;
      // if( (pN1>=mypart && pS1>=mypart) )
      {
        ++jelem;
        ++nelems;        
      }
    }
    region.nb_lat_elems[jlat] = jelem;
  }
  
  int nb_region_nodes = 0;
  for( int jlat=lat_north; jlat<=lat_south; ++jlat ) {
    nb_region_nodes += region.lat_end[jlat]-region.lat_begin[jlat]+1;
  }
  std::cout << "nb_region_nodes = " << nb_region_nodes << std::endl;
  region.nnodes = nb_region_nodes;
  
}

Mesh* RGGMeshGenerator::generate_mesh(const RGG& rgg,const std::vector<int>& parts, const Region& region)
{
  double tol = 1e-3;
  Mesh* mesh = new Mesh();
  int mypart = options.get<int>("part");
  int n, l;
  
  std::cout << "nnodes = "  << region.nnodes  << std::endl;
  std::cout << "nquads = "  << region.nquads  << std::endl;
  std::cout << "ntriags = " << region.ntriags << std::endl;

  std::vector<int> offset_glb(rgg.nlat(),0);
  std::vector<int> offset_loc(region.south-region.north+1,0);

  n=0;
  for( int jlat=0; jlat<rgg.nlat(); ++jlat )
  {
    offset_glb[jlat]=n;
    n+=rgg.nlon(jlat);
  };
  
  
  std::vector<int> extents = Extents(Field::UNDEF_VARS, region.nnodes);
  FunctionSpace& nodes = mesh->add_function_space( new FunctionSpace("nodes","LagrangeP1",extents) );
  nodes.metadata().set("type",static_cast<int>(Entity::NODES));
  ArrayView<double,2> coords( nodes.create_field<double>("coordinates",2) );
  ArrayView<int,1>    glb_idx( nodes.create_field<int>("glb_idx",1) );
  ArrayView<int,1>    master_glb_idx( nodes.create_field<int>("master_glb_idx",1) );
  ArrayView<int,1>    proc( nodes.create_field<int>("proc",1) );
  
  std::cout << "s1 = " << coords.strides()[0] << std::endl;
  std::cout << "s2 = " << coords.strides()[1] << std::endl;

  std::cout << "e1 = " << coords.extents()[0] << std::endl;
  std::cout << "e2 = " << coords.extents()[1] << std::endl;


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
      master_glb_idx(jnode) = glb_idx(jnode);
      proc(jnode) = parts[n];
      ++jnode;
    }
  };
  
  extents = Extents(Field::UNDEF_VARS, region.nquads);
  FunctionSpace& quads = mesh->add_function_space( new FunctionSpace("quads","LagrangeP1",extents) );
  quads.metadata().set("type",static_cast<int>(Entity::ELEMS));
  ArrayView<int,2> quad_nodes( quads.create_field<int>("nodes",4) );
  ArrayView<int,1> quad_glb_idx( quads.create_field<int>("glb_idx",1) );
  ArrayView<int,1> quad_master_glb_idx( quads.create_field<int>("master_glb_idx",1) );
  ArrayView<int,1> quad_proc( quads.create_field<int>("proc",1) );
  
  extents = Extents(Field::UNDEF_VARS, region.ntriags);
  FunctionSpace& triags = mesh->add_function_space( new FunctionSpace("triags","LagrangeP1",extents) );
  triags.metadata().set("type",static_cast<int>(Entity::ELEMS));
  ArrayView<int,2> triag_nodes( triags.create_field<int>("nodes",3) );
  ArrayView<int,1> triag_glb_idx( triags.create_field<int>("glb_idx",1) );
  ArrayView<int,1> triag_master_glb_idx( triags.create_field<int>("master_glb_idx",1) );
  ArrayView<int,1> triag_proc( triags.create_field<int>("proc",1) );

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
        quad_nodes(jquad,0) = F_IDX(offset_loc[ilatN] + elem[0] - region.lat_begin[jlatN]);
        quad_nodes(jquad,1) = F_IDX(offset_loc[ilatS] + elem[1] - region.lat_begin[jlatS]);
        quad_nodes(jquad,2) = F_IDX(offset_loc[ilatS] + elem[2] - region.lat_begin[jlatS]);
        quad_nodes(jquad,3) = F_IDX(offset_loc[ilatN] + elem[3] - region.lat_begin[jlatN]);
        // std::cout << quad_nodes(0,jquad) << " " << quad_nodes(1,jquad) << " " << quad_nodes(2,jquad) << " " << quad_nodes(3,jquad) << std::endl;
        quad_glb_idx(jquad) = jquad+jtriag+1;
        quad_master_glb_idx(jquad) = quad_glb_idx(jquad);
        int p0 = proc(C_IDX(quad_nodes(jquad,0)));
        int p1 = proc(C_IDX(quad_nodes(jquad,1)));
        int p2 = proc(C_IDX(quad_nodes(jquad,2)));
        int p3 = proc(C_IDX(quad_nodes(jquad,3)));
        quad_proc(jquad) = mypart; //std::min(p0,std::min(p1,std::max(p2,p3)));
        ++jquad;
      }
      else // This is a triag
      {
        if(elem[3]<0) // This is a triangle pointing up
        {
          triag_nodes(jtriag,0) = F_IDX(offset_loc[ilatN] + elem[0] - region.lat_begin[jlatN]);
          triag_nodes(jtriag,1) = F_IDX(offset_loc[ilatS] + elem[1] - region.lat_begin[jlatS]);
          triag_nodes(jtriag,2) = F_IDX(offset_loc[ilatS] + elem[2] - region.lat_begin[jlatS]);
        }
        else // This is a triangle pointing down
        {
          triag_nodes(jtriag,0) = F_IDX(offset_loc[ilatN] + elem[0] - region.lat_begin[jlatN]);
          triag_nodes(jtriag,1) = F_IDX(offset_loc[ilatS] + elem[1] - region.lat_begin[jlatS]);
          triag_nodes(jtriag,2) = F_IDX(offset_loc[ilatN] + elem[3] - region.lat_begin[jlatN]);
        }
        triag_glb_idx(jtriag) = jquad+jtriag+1;
        triag_master_glb_idx(jtriag) = triag_glb_idx(jtriag);
        int p0 = proc(C_IDX(triag_nodes(jtriag,0)));
        int p1 = proc(C_IDX(triag_nodes(jtriag,1)));
        int p2 = proc(C_IDX(triag_nodes(jtriag,2)));
        triag_proc(jtriag) = mypart; //std::max(p0,std::max(p1,p2));
        ++jtriag;
      }
    }
  }

  extents[1]=0;
  FunctionSpace& edges = mesh->add_function_space( new FunctionSpace("edges","LagrangeP1",extents) );
  edges.metadata().set("type",static_cast<int>(Entity::FACES));
  edges.create_field<int>("nodes",2);
  edges.create_field<int>("glb_idx",1);
  edges.create_field<int>("master_glb_idx",1);
  edges.create_field<int>("proc",1);
  
  nodes.metadata().set("nb_owned",rgg.ngptot());
  quads.metadata().set("nb_owned",region.nquads);
  triags.metadata().set("nb_owned",region.ntriags);
  edges.metadata().set("nb_owned",0);
  
  nodes.metadata().set("max_glb_idx",rgg.ngptot()+2*rgg.nlat()+1);
  quads.metadata().set("max_glb_idx",region.nquads+region.ntriags+1);
  triags.metadata().set("max_glb_idx",region.nquads+region.ntriags+1);
  edges.metadata().set("max_glb_idx",region.nquads+region.ntriags+1);
  
  return mesh;
}

} // namespace meshgen
} // namespace atlas

