#include <stdexcept>
#include <cmath>
#include <set>

#include "Mesh.hpp"
#include "FunctionSpace.hpp"
#include "Field.hpp"
#include "BuildDualMesh.hpp"
#include "Parameters.hpp"

namespace ecmwf {

void build_centroids( FunctionSpace& func_space, FieldT<double>& coords)
{
  int nb_elems = func_space.bounds()[1];
  FieldT<int>& elem_nodes = func_space.field<int>( "nodes" );
  int nb_nodes_per_elem = elem_nodes.bounds()[0];
  FieldT<int>& elem_glb_idx = func_space.field<int>( "glb_idx" );
  FieldT<double>& elem_centroids = func_space.create_field<double>( "centroids", 2 );
  for (int e=0; e<nb_elems; ++e)
  {
    elem_centroids(XX,e) = 0.;
    elem_centroids(YY,e) = 0.;
    for (int n=0; n<nb_nodes_per_elem; ++n)
    {
      elem_centroids(XX,e) += coords(XX, C_IDX(elem_nodes(n,e)) );
      elem_centroids(YY,e) += coords(YY, C_IDX(elem_nodes(n,e)) );
    }
    elem_centroids(XX,e) /= static_cast<double>(nb_nodes_per_elem);
    elem_centroids(YY,e) /= static_cast<double>(nb_nodes_per_elem);
  }
}

void add_dual_volume_contribution(
    FunctionSpace& elements,
    FunctionSpace& edges,
    FieldT<double>& dual_volumes )
{
  int nb_elems = elements.bounds()[1];
  FieldT<double>& elem_centroids = elements.field<double>("centroids");
  FieldT<int>& elem_to_edges = elements.field<int>("to_edge");
  FieldT<double>& edge_centroids = edges.field<double>("centroids");
  FieldT<int>& edge_nodes = edges.field<int>("nodes");
  FieldT<int>& edge_to_elem = edges.field<int>("to_elem");
  FieldT<double>& node_coords = dual_volumes.function_space().field<double>("coordinates");
  int nb_edges_per_elem = elem_to_edges.bounds()[0];
  for (int elem=0; elem<nb_elems; ++elem)
  {
    double x0 = elem_centroids(XX,elem);
    double y0 = elem_centroids(YY,elem);
    for (int jedge=0; jedge<nb_edges_per_elem; ++jedge)
    {
      int edge = C_IDX(elem_to_edges(jedge,elem));
      double x1 = edge_centroids(XX,edge);
      double y1 = edge_centroids(YY,edge);
      for( int j=0; j<2; ++j )
      {
        int node = C_IDX(edge_nodes(j,edge));
        double x2 = node_coords(XX,node);
        double y2 = node_coords(YY,node);
        double triag_area = std::abs( x0*(y1-y2)+x1*(y2-y0)+x2*(y0-y1) )*0.5;
        dual_volumes(0,node) += triag_area;
      }
    }
  }
}

void add_dual_volume_contribution(
    FunctionSpace& edges,
    FieldT<double>& dual_volumes )
{
  FieldT<int>& node_glb_idx = dual_volumes.function_space().field<int>("glb_idx");
  FieldT<double>& edge_centroids = edges.field<double>("centroids");
  FieldT<int>& edge_nodes = edges.field<int>("nodes");
  FieldT<int>& edge_glb_idx = edges.field<int>("glb_idx");
  FieldT<int>& edge_to_elem = edges.field<int>("to_elem");
  FieldT<double>& node_coords = dual_volumes.function_space().field<double>("coordinates");
  int nb_edges = edges.bounds()[1];
  std::map<int,std::vector<int> > node_to_bdry_edge;
  for(int edge=0; edge<nb_edges; ++edge)
  {
    if (C_IDX(edge_to_elem(0,edge)) >= 0 && C_IDX(edge_to_elem(3,edge)) < 0)
    {
      node_to_bdry_edge[ C_IDX(edge_nodes(0,edge)) ].push_back(edge);
      node_to_bdry_edge[ C_IDX(edge_nodes(1,edge)) ].push_back(edge);
    }
  }

  double pi = acos(-1.);
  double tol = 1.e-6;
  double ymax = node_coords.function_space().metadata<double>("ymax");
  double ymin = node_coords.function_space().metadata<double>("ymin");
  double dx   = node_coords.function_space().metadata<double>("xmax")
                - node_coords.function_space().metadata<double>("xmin");
  std::set<int> visited_nodes;

  std::map<int,std::vector<int> >::iterator it;
  for( it=node_to_bdry_edge.begin(); it!=node_to_bdry_edge.end(); ++it)
  {
    int node = (*it).first;
    std::vector<int>& bdry_edges = (*it).second;
    double x0 = node_coords(XX,node);
    double y0 = node_coords(YY,node);
    double x1,y1, x2,y2;
    for (int jedge=0; jedge<bdry_edges.size(); ++jedge)
    {
      int edge = bdry_edges[jedge];
      x1 = edge_centroids(XX,edge);
      y1 = edge_centroids(YY,edge);
      x2 = x1;
      y2 = 0.;
      if ( std::abs(y1-ymax)<tol )
        y2 = 0.5*pi;
      else if ( std::abs(y1-ymin)<tol )
        y2 = -0.5*pi;

      if( y2!=0 )
      {
        //std::cout << "adding contribution for node " << node_glb_idx(node) << " and edge " << edge_glb_idx(edge) << std::endl;
        double triag_area = std::abs( x0*(y1-y2)+x1*(y2-y0)+x2*(y0-y1) )*0.5;
        dual_volumes(0,node) += triag_area;
      }
    }
  }
}


void build_dual_normals( Mesh& mesh )
{
  std::vector< FieldT<double>* > elem_centroids( mesh.nb_function_spaces() );
  for (int func_space_idx=0; func_space_idx<mesh.nb_function_spaces(); ++func_space_idx)
  {
    FunctionSpace& func_space = mesh.function_space(func_space_idx);
    if( func_space.has_field("centroids") )
      elem_centroids[func_space_idx] = &func_space.field<double>("centroids");
  }

  FunctionSpace&  nodes_2d = mesh.function_space("nodes_2d");
  FieldT<double>& node_coords = nodes_2d.field<double>("coordinates");
  double ymax = nodes_2d.metadata<double>("ymax");
  double ymin = nodes_2d.metadata<double>("ymin");
  double pi = acos(-1.);
  double tol = 1.e-6;

  double xl, yl, xr, yr;
  FunctionSpace&  edges = mesh.function_space("edges");
  FieldT<int>&    edge_to_elem = edges.field<int>("to_elem");
  FieldT<int>&    edge_nodes = edges.field<int>("nodes");
  FieldT<double>& edge_centroids = edges.field<double>("centroids");
  FieldT<double>& dual_normals = edges.create_field<double>("dual_normals",2);
  int nb_edges = edges.bounds()[1];

  std::map<int,std::vector<int> > node_to_bdry_edge;
  for(int edge=0; edge<nb_edges; ++edge)
  {
    if (C_IDX(edge_to_elem(0,edge)) >= 0 && C_IDX(edge_to_elem(3,edge)) < 0)
    {
      node_to_bdry_edge[ C_IDX(edge_nodes(0,edge)) ].push_back(edge);
      node_to_bdry_edge[ C_IDX(edge_nodes(1,edge)) ].push_back(edge);
    }
  }

  for (int edge=0; edge<nb_edges; ++edge)
  {
    if( edge_to_elem(0,edge) < 0 )
    {
      // this is a pole edge
      // only compute for one node
      for (int n=0; n<2; ++n)
      {
        int node = edge_nodes(0,edge);
        std::vector<int>& bdry_edges = node_to_bdry_edge[node];
        double x[2];
        int cnt=0;
        for (int jedge=0; jedge<bdry_edges.size(); ++jedge)
        {
          int bdry_edge = bdry_edges[jedge];
          if ( std::abs(edge_centroids(YY,bdry_edge)-ymax)<tol )
          {
            edge_centroids(YY,edge) = 0.5*pi;
            x[cnt] = edge_centroids(XX,bdry_edge);
            ++cnt;
          }
          else if ( std::abs(edge_centroids(YY,bdry_edge)-ymin)<tol )
          {
            edge_centroids(YY,edge) = -0.5*pi;
            x[cnt] = edge_centroids(XX,bdry_edge);
            ++cnt;
          }
        }
        if (cnt == 2 )
        {
          dual_normals(XX,edge) = 0;
          dual_normals(YY,edge) = -x[1]+x[0];
          //std::cout << "pole dual_normal = " << dual_normals(YY,edge) << std::endl;
          break;
        }
      }
    }
    else
    {
      int left_func_space_idx  = C_IDX(edge_to_elem(0,edge));
      int left_elem            = C_IDX(edge_to_elem(1,edge));
      int right_func_space_idx = C_IDX(edge_to_elem(2,edge));
      int right_elem           = C_IDX(edge_to_elem(3,edge));
      xl = (*elem_centroids[left_func_space_idx])(XX,left_elem);
      yl = (*elem_centroids[left_func_space_idx])(YY,left_elem);
      if( right_elem < 0 )
      {
        xr = edge_centroids(XX,edge);
        yr = edge_centroids(YY,edge);;
        if ( std::abs(yr-ymax)<tol )
          yr = 0.5*pi;
        else if( std::abs(yr-ymin)<tol )
          yr = -0.5*pi;
      }
      else
      {
        xr = (*elem_centroids[right_func_space_idx])(XX,right_elem);
        yr = (*elem_centroids[right_func_space_idx])(YY,right_elem);
      }
      dual_normals(XX,edge) =  yl-yr;
      dual_normals(YY,edge) = -xl+xr;
    }
  }
}

void build_dual_mesh( Mesh& mesh )
{
  FunctionSpace& nodes_2d   = mesh.function_space( "nodes_2d" );
  FieldT<double>& coords    = nodes_2d.field<double>( "coordinates" );
  FieldT<int>& glb_idx      = nodes_2d.field<int>( "glb_idx" );
  FieldT<int>& master_glb_idx      = nodes_2d.field<int>( "master_glb_idx" );
  FieldT<int>& proc      = nodes_2d.field<int>( "proc" );
  FieldT<double>& dual_volumes = nodes_2d.create_field<double>( "dual_volumes", 1 );
  int nb_nodes = nodes_2d.bounds()[1];

  FunctionSpace& quads       = mesh.function_space( "quads" );
  FunctionSpace& triags      = mesh.function_space( "triags" );
  FunctionSpace& edges       = mesh.function_space( "edges" );
  FieldT<int>& edge_proc  = edges.field<int>( "proc" );
  FieldT<int>& edge_glb_idx  = edges.field<int>( "glb_idx" );
  FieldT<int>& edge_master_glb_idx  = edges.field<int>( "master_glb_idx" );

  build_centroids(quads,  coords);
  build_centroids(triags, coords);
  build_centroids(edges,  coords);

  for (int node=0; node<nb_nodes; ++node)
    dual_volumes(0,node) = 0.;

  add_dual_volume_contribution( quads, edges, dual_volumes );
  add_dual_volume_contribution( triags, edges, dual_volumes );

  add_dual_volume_contribution( edges, dual_volumes );

  build_dual_normals( mesh );



//  std::cout << "proc" << std::endl;
//  for (int node=0; node<nb_nodes; ++node)
//    std::cout << glb_idx(node) << "  :  " << proc(0,node) << std::endl;

//  std::cout << "master_glb_idx" << std::endl;
//  for (int node=0; node<nb_nodes; ++node)
//    std::cout << glb_idx(node) << "  :  " << master_glb_idx(0,node) << std::endl;


  nodes_2d.parallelise();
  edges.parallelise();


//  std::cout << "dual_volumes" << std::endl;
//  for (int node=0; node<nb_nodes; ++node)
//    std::cout << glb_idx(node) << "  :  " << dual_volumes(0,node) << std::endl;

  dual_volumes.halo_exchange();

//  std::cout << "dual_volumes" << std::endl;
//  for (int node=0; node<nb_nodes; ++node)
//    std::cout << glb_idx(node) << "  :  " << dual_volumes(0,node) << std::endl;

  FieldT<double>& dual_normals = edges.field<double>("dual_normals");
  int nb_edges = edges.bounds()[1];
//  std::cout << "dual_normals" << std::endl;
//  for (int edge=0; edge<nb_edges; ++edge)
//    std::cout << edge_glb_idx(edge) << "  :  " << dual_normals(XX,edge) << " " << dual_normals(YY,edge) << std::endl;


  dual_normals.halo_exchange();

//  std::cout << "dual_normals" << std::endl;
//  for (int edge=0; edge<nb_edges; ++edge)
//    std::cout << edge_glb_idx(edge) << "  :  " << dual_normals(XX,edge) << " " << dual_normals(YY,edge) << std::endl;


}


// ------------------------------------------------------------------
// C wrapper interfaces to C++ routines

void ecmwf__build_dual_mesh ( Mesh* mesh) {
  build_dual_mesh(*mesh);
}

// ------------------------------------------------------------------


} // namespace ecmwf

