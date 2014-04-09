// (C) Copyright 1996-2014 ECMWF.

#include "atlas/meshgen/RGG.hpp"
#include "atlas/Mesh.hpp"
#include "atlas/FunctionSpace.hpp"
#include "atlas/Parameters.hpp"
#include "atlas/Field.hpp"
#include <numeric>
#include <cmath>

namespace atlas {
namespace meshgen {

double spherical_distance(double x1,double y1,double x2,double y2,double rad)
{
  using namespace std;
  if((std::abs(x2-x1))<1e-8) return (y2-y1)*rad;
  else if((std::abs(y2-y1))<1e-8) return (x2-x1)*rad;
  else return acos(cos(x1)*cos(y1)*cos(x2)*cos(y2)+cos(x1)*sin(y1)*cos(x2)*sin(y2)+sin(x1)*sin(x2))*rad;
}

int RGG::ngptot() const
{
  return 2*std::accumulate(lon_.data(),lon_.data()+lon_.size(),0);
}    

RGGMeshGenerator::RGGMeshGenerator()
{
  options.set("include_pole",false);
}

Mesh* RGGMeshGenerator::generate(const RGG& rgg)
{
  double beta = M_PI/6.;
  double tol = 1e-3;
  Mesh* mesh = new Mesh();
  
  
  std::vector<int> bounds(2); 
  bounds[0] = Field::UNDEF_VARS;
  bounds[1] = rgg.ngptot()+2*rgg.nlat();

  FunctionSpace& nodes = mesh->add_function_space( new FunctionSpace("nodes","LagrangeP1",bounds) );
  nodes.metadata().set("type",static_cast<int>(Entity::NODES));
  FieldT<double>& coords  = nodes.create_field<double>("coordinates",2);
  FieldT<int>&    glb_idx = nodes.create_field<int>("glb_idx",1);
  FieldT<int>&    master_glb_idx = nodes.create_field<int>("master_glb_idx",1);
  FieldT<int>&    proc = nodes.create_field<int>("proc",1);

  FunctionSpace& quads = mesh->add_function_space( new FunctionSpace("quads","LagrangeP1",bounds) );
  quads.metadata().set("type",static_cast<int>(Entity::ELEMS));
  FieldT<int>& quad_nodes = quads.create_field<int>("nodes",4);
  FieldT<int>& quad_glb_idx = quads.create_field<int>("glb_idx",1);
  FieldT<int>& quad_master_glb_idx = quads.create_field<int>("master_glb_idx",1);
  FieldT<int>& quad_proc = quads.create_field<int>("proc",1);
  
  FunctionSpace& triags = mesh->add_function_space( new FunctionSpace("triags","LagrangeP1",bounds) );
  triags.metadata().set("type",static_cast<int>(Entity::ELEMS));
  FieldT<int>& triag_nodes = triags.create_field<int>("nodes",3);
  FieldT<int>& triag_glb_idx = triags.create_field<int>("glb_idx",1);
  FieldT<int>& triag_master_glb_idx = triags.create_field<int>("master_glb_idx",1);
  FieldT<int>& triag_proc = triags.create_field<int>("proc",1);


  std::vector<int> offset; offset.reserve(2*rgg.nlat());
  std::vector<int> length; length.reserve(2*rgg.nlat());
  
  std::vector<int> periodic; periodic.reserve(2*rgg.nlat());
  
  int cnt=0;
  for (int jlat=0; jlat<rgg.nlat(); ++jlat)
  {
    int ilat = std::abs(rgg.nlat()-1-jlat);
    offset.push_back(cnt);
    length.push_back(rgg.nlon(ilat));
    for (int jlon=0; jlon<rgg.nlon(ilat); ++jlon)
    {
      coords(XX,cnt) = rgg.lon(jlon,ilat);
      coords(YY,cnt) = rgg.lat(ilat);
      glb_idx(cnt)   = cnt+1;
      master_glb_idx(cnt) = glb_idx(cnt);
      proc(cnt) = 0;
      ++cnt;
    }
  }
  for (int jlat=0; jlat<rgg.nlat(); ++jlat)
  {
    int ilat = jlat;
    offset.push_back(cnt);
    length.push_back(rgg.nlon(ilat));
    for (int jlon=0; jlon<rgg.nlon(ilat); ++jlon)
    {
      coords(XX,cnt) = rgg.lon(jlon,ilat);
      coords(YY,cnt) = -rgg.lat(ilat);
      glb_idx(cnt)   = cnt+1;
      master_glb_idx(cnt) = glb_idx(cnt);
      proc(cnt) = 0;
      ++cnt;
    }
  }
  
  for (int jlat=0; jlat<rgg.nlat(); ++jlat)
  {
    periodic.push_back(cnt);
    int ilat = std::abs(rgg.nlat()-1-jlat);
    coords(XX,cnt) = 2.*M_PI;
    coords(YY,cnt) = rgg.lat(ilat);
    glb_idx(cnt)   = cnt+1;
    master_glb_idx(cnt) = glb_idx(offset[ilat]);
    proc(cnt) = 0;
    ++cnt;
  }
  for (int jlat=0; jlat<rgg.nlat(); ++jlat)
  {
    periodic.push_back(cnt);
    int ilat = jlat;
    coords(XX,cnt) = 2.*M_PI;
    coords(YY,cnt) = -rgg.lat(ilat);
    glb_idx(cnt)   = cnt+1;
    master_glb_idx(cnt) = glb_idx(offset[ilat]);
    proc(cnt) = 0;
    ++cnt;
  }
    
  int nquads=0;
  int ntriags=0;
  int nelems=0;

  for (int jlat=0; jlat<2*rgg.nlat()-1; ++jlat)
  {
    int il2 = jlat;
    int il1 = jlat+1;
    int ip1 = offset[il2];
    int ip2 = offset[il1];
    int ip3 = ip2+1;
    int ip4 = ip1+1;
    
    int lon2 = offset[il2] + length[il2]-1;
    int lon1 = offset[il1] + length[il1]-1;
    int last_triag=-1;
    
    int ilon=0;
    while (ip2 != periodic[il1] && ip1 != periodic[il2])      
    {
      bool make_triangle_up   = false;
      bool make_triangle_down = false;
      bool make_quad = false;
      
      double d31 = spherical_distance(coords(XX,ip3),coords(YY,ip3),coords(XX,ip1),coords(YY,ip1),1.);
      double d42 = spherical_distance(coords(XX,ip4),coords(YY,ip4),coords(XX,ip2),coords(YY,ip2),1.);
      double d43 = spherical_distance(coords(XX,ip4),coords(YY,ip4),coords(XX,ip3),coords(YY,ip3),1.);
      if (d31 < d43 && d31 < d42) make_triangle_up = true;
      else if (d42 < d43 && d42 < d31) make_triangle_down = true;
      else make_quad = true;
      // double dx = spherical_distance(coords(XX,ip3),coords(YY,ip4),coords(XX,ip4),coords(YY,ip4),1.);
      // double dy = spherical_distance(coords(XX,ip4),coords(YY,ip3),coords(XX,ip4),coords(YY,ip4),1.);
      // double alpha = std::atan2(dx,std::abs(dy));
      // il2  1--------+------4
      //      |        |alph/
      //      |        |->/
      // il1  2--------3
      // bool make_triangle_up   = ( alpha > beta  || (alpha >= tol && last_triag==1 ) );
      // bool make_triangle_down = ( alpha < -beta || (alpha <= tol && last_triag==0 ) );
      // bool make_quad = ( ! make_triangle_up && ! make_triangle_down );
      if( make_quad )
      {
        // add quadrilateral
        quad_nodes(0,nquads) = F_IDX(ip1);
        quad_nodes(1,nquads) = F_IDX(ip2);
        quad_nodes(2,nquads) = F_IDX(ip3);
        quad_nodes(3,nquads) = F_IDX(ip4);
        quad_glb_idx(nquads) = nelems+1;
        quad_master_glb_idx(nquads) = quad_glb_idx(nquads);
        quad_proc(nquads) = 0;
        ++nquads;
        ip1=ip4;
        ip2=ip3;
        last_triag=-1;
        // std::cout << "quad " << ip1 << "  " << ip2 << " "<< ip3 << " " << ip4 <<std::endl;
      }
      else if( make_triangle_down  ) // make triangle down
      {
        // triangle without ip3
        triag_nodes(0,ntriags) = F_IDX(ip1);
        triag_nodes(1,ntriags) = F_IDX(ip2);
        triag_nodes(2,ntriags) = F_IDX(ip4);
        triag_glb_idx(ntriags) = nelems+1;
        triag_master_glb_idx(ntriags) = triag_glb_idx(ntriags);
        triag_proc(ntriags) = 0;
        
        ++ntriags;
        ip1=ip4;
        ip2=ip2;
        last_triag=0;
      }
      else // make triangle up
      {
        // triangle without ip4
        triag_nodes(0,ntriags) = F_IDX(ip1);
        triag_nodes(1,ntriags) = F_IDX(ip2);
        triag_nodes(2,ntriags) = F_IDX(ip3);
        triag_glb_idx(ntriags) = nelems+1;
        ++ntriags;
        ip1=ip1;
        ip2=ip3;
        last_triag=1;
      }
      ++nelems;
      ip3 = ip2+1;
      ip4 = ip1+1;
      if( ip1 == lon2 ){
        ip4 = periodic[il2];
      }
      if( ip2 == lon1 ){
        ip3 = periodic[il1];
      } 
    }
  }
  
  std::cout << "nquads = " << nquads << std::endl;
  std::cout << "ntriags = " << ntriags << std::endl;
  std::cout << "nelems = " << nelems << std::endl;
  
  bounds[1] = nquads;
  quads.resize(bounds);
  bounds[1] = ntriags;
  triags.resize(bounds);
  
  bounds[1] = 0;
  FunctionSpace& edges = mesh->add_function_space( new FunctionSpace("edges","LagrangeP1",bounds) );
  edges.metadata().set("type",static_cast<int>(Entity::FACES));
  FieldT<int>& edge_nodes = edges.create_field<int>("nodes",2);
  FieldT<int>& edge_glb_idx = edges.create_field<int>("glb_idx",1);
  FieldT<int>& edge_master_glb_idx = edges.create_field<int>("master_glb_idx",1);
  FieldT<int>& edge_proc = edges.create_field<int>("proc",1);
  
  nodes.metadata().set("nb_owned",rgg.ngptot());
  quads.metadata().set("nb_owned",nquads);
  triags.metadata().set("nb_owned",ntriags);
  edges.metadata().set("nb_owned",0);
  
  nodes.metadata().set("max_glb_idx",rgg.ngptot()+2*rgg.nlat()+1);
  quads.metadata().set("max_glb_idx",nelems+1);
  triags.metadata().set("max_glb_idx",nelems+1);
  edges.metadata().set("max_glb_idx",nelems+1);
  
  return mesh;
}

} // namespace meshgen
} // namespace atlas

