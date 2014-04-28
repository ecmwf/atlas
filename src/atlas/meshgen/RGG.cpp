// (C) Copyright 1996-2014 ECMWF.

#include "atlas/meshgen/RGG.hpp"
#include "atlas/Mesh.hpp"
#include "atlas/FunctionSpace.hpp"
#include "atlas/Parameters.hpp"
#include "atlas/Field.hpp"
#include "atlas/Array.hpp"
#include "atlas/ArrayView.hpp"
#include "atlas/meshgen/EqualAreaPartitioner.hpp"
#include <numeric>
#include <cmath>

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
  region.lat_begin.resize(rgg.nlat(),-1);
  region.lat_end.resize(rgg.nlat(),-1);
  region.nb_lat_elems.resize(rgg.nlat(),-1);
  for( int jlat=lat_north; jlat<=lat_south; ++jlat ) {
    for( int jlon=0; jlon<rgg.nlon(jlat); ++jlon) {
      if( region.lat_begin[jlat]<0 && parts[offset[jlat]+jlon] == mypart)
      {
        region.lat_begin[jlat]=jlon;
        break;
      }
    }
    for( int jlon=rgg.nlon(jlat)-1; jlon>=0; --jlon) {
      if( region.lat_end[jlat]<0 && parts[offset[jlat]+jlon] == mypart)
      {
        region.lat_end[jlat]=jlon;
        break;
      }
    }
    std::cout << "lat["<<jlat<<"] :  [" << region.lat_begin[jlat] << "("<< rgg.lon(region.lat_begin[jlat],jlat) << ") -- " << region.lat_end[jlat] << "("<< rgg.lon(region.lat_end[jlat],jlat) << ")]" << std::endl;
  }
  
  /*
  We need to connect to next region
  */
  for( int jlat=lat_north; jlat<=lat_south; ++jlat ) {
    region.lat_end[jlat] = std::min(rgg.nlon(jlat), region.lat_end[jlat]+1);
  }

  region.north = lat_north;
  region.south = lat_south;
  
  std::vector<int> extents = Extents(region.south-region.north, 2*rgg.nlonmax(), 4);
  std::cout << "allocating elems" <<  "(" << extents[0] << "," << extents[1] << "," << extents[2] << ")" << std::endl;
  region.elems.resize(extents);
  region.elems = -1;
  
  int nelems=0;
  region.nquads=0;
  region.ntriags=0;

  ArrayView<int,3> elemview(region.elems);
  
  for (int jlat=region.north; jlat<region.south; ++jlat)
  {
    int ilat=jlat-region.north;
    
    ArrayView<int,2> lat_elems_view = elemview[ilat];
    
    
    int latN = jlat;
    int latS = jlat+1;
    double yN = rgg.lat(latN);
    double yS = rgg.lat(latS);
    

    // Initial points
    int ipN1 = region.lat_begin[latN];
    double xN1 = rgg.lon(ipN1,latN);
    int ipS1;
    double dist=2.*M_PI;
    for (int jlon=std::max(0,region.lat_begin[latS]-1); jlon<=region.lat_end[latS]; ++jlon )
    {
      double xS = rgg.lon(jlon,latS);
      double dist_new = std::abs(xN1-xS);
      if (dist_new < dist)
      {
        ipS1 = jlon;
        dist = dist_new;
      }
      else
      {
        break;
      }
    }

    std::cout << "begin strip["<<jlat<<"] : connect "<< ipN1 << " with " << ipS1 << std::endl;

    int ipN2 = ipN1+1;
    int ipS2 = ipS1+1;

    int jelem=0;
    std::cout << "=================" << std::endl;
    // continue until ipN1 and ipS1 are no longer of my part, or the east domain boundary
    while ( true )
    {
      region.lat_end[latN] = std::max(ipN1,region.lat_end[latN]);
      region.lat_end[latS] = std::max(ipS1,region.lat_end[latS]);
      int pN1 = parts[offset[latN]+ipN1];
      int pS1 = parts[offset[latS]+ipS1];
      if ( (pN1>mypart && pS1>mypart) || (ipN1==rgg.nlon(latN)&&ipS1==rgg.nlon(latS)) )
        break;
      
      if( (pN1==mypart && pS1==mypart) )
      {
        // region.lat_begin[latN] = std::min(ipN1,region.lat_begin[latN]);
        // region.lat_begin[latS] = std::min(ipS1,region.lat_begin[latS]);
      }
      
      double xN1 = rgg.lon(ipN1,latN);
      double xN2 = rgg.lon(ipN2,latN);
      double xS1 = rgg.lon(ipS1,latS);
      double xS2 = rgg.lon(ipS2,latS);
  
      std::cout << "-------" << std::endl;
      // std::cout << "  access  " << region.elems.stride(0)*(jlat-region.north) + region.elems.stride(1)*jelem + 5 << std::endl;
      std::cout << ipN1 << "("<< xN1 << ")  " << ipN2 <<  "("<< xN2 << ")  " << std::endl;
      std::cout << ipS1 << "("<< xS1 << ")  " << ipS2 <<  "("<< xS2 << ")  " << std::endl;
      bool make_triangle_up   = false;
      bool make_triangle_down = false;
      bool make_quad = false;
    
      
      double dN1S2 = std::abs(spherical_distance(xN1,yN,xS2,yS,1.));
      double dS1N2 = std::abs(spherical_distance(xS1,yS,xN2,yN,1.));
      double dN2S2 = std::abs(spherical_distance(xN2,yN,xS2,yS,1.));
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
        // std::cout << "          " << ip1 << "  " << ip4 << std::endl;
        // std::cout << "          " << ip2 << "  " << ip3;
        elem[0] = ipN1;
        elem[1] = ipS1;
        elem[2] = ipS2;
        elem[3] = ipN2;
        ipN1=std::min(rgg.nlon(latN),ipN2);
        ipS1=std::min(rgg.nlon(latS),ipS2);
        // std::cout << "  .  " << std::endl;
        if( (pN1>=mypart && pS1>=mypart) )
          ++region.nquads;
      }
      else if( make_triangle_down  ) // make triangle down
      {
        // triangle without ip3
        // std::cout << "          " << ip1 << "  " << ip4 << std::endl;
        // std::cout << "          " << ip2;
        elem[0] = ipN1;
        elem[1] = ipS1;
        elem[2] = -1;
        elem[3] = ipN2;
        ipN1=std::min(rgg.nlon(latN),ipN2);
        ipS1=std::min(rgg.nlon(latS),ipS1);
        // std::cout << "  ." << std::endl;
        if( (pN1>=mypart && pS1>=mypart) )
          ++region.ntriags;
      }
      else // make triangle up
      {
        // triangle without ip4
        // std::cout << "          " << ip1 << std::endl;
        // std::cout << "          " << ip2 << "  " << ip3;
        elem[0] = ipN1;
        elem[1] = ipS1;
        elem[2] = ipS2;
        elem[3] = -1;
        ipN1=std::min(rgg.nlon(latN),ipN1);
        ipS1=std::min(rgg.nlon(latS),ipS2);
        // std::cout << "  . "  << std::endl;
        if( (pN1>=mypart && pS1>=mypart) )
          ++region.ntriags;
      }
      ipN2 = ipN1+1;
      ipS2 = ipS1+1;
      if( (pN1>=mypart && pS1>=mypart) )
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
        quad_proc(jquad) = std::min(p0,std::min(p1,std::min(p2,p3)));
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
        triag_proc(jtriag) = std::min(p0,std::min(p1,p2));
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

