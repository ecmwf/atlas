/*
 * (C) Copyright 1996-2014 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include <stdexcept>
#include <cmath>
#include <set>
#include <limits>
#include <iostream>
#include <algorithm>    // std::sort

#include "atlas/mpl/Checksum.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/mesh/FunctionSpace.h"
#include "atlas/mesh/Field.h"
#include "atlas/actions/BuildDualMesh.h"
#include "atlas/mesh/Parameters.h"
#include "atlas/mesh/Util.h"
#include "atlas/util/ArrayView.h"
#include "atlas/util/IndexView.h"

namespace atlas {
namespace actions {


namespace {

void global_bounding_box( FunctionSpace& nodes, double min[2], double max[2] )
{
  ArrayView<double,2> latlon( nodes.field("coordinates") );
  const int nb_nodes = nodes.shape(0);
  min[XX] =  std::numeric_limits<double>::max();
  min[YY] =  std::numeric_limits<double>::max();
  max[XX] = -std::numeric_limits<double>::max();
  max[YY] = -std::numeric_limits<double>::max();
  for (int node=0; node<nb_nodes; ++node)
  {
    min[XX] = std::min( min[XX], latlon(node,XX) );
    min[YY] = std::min( min[YY], latlon(node,YY) );
    max[XX] = std::max( max[XX], latlon(node,XX) );
    max[YY] = std::max( max[YY], latlon(node,YY) );
  }
  MPL_CHECK_RESULT( MPI_Allreduce( MPI_IN_PLACE, &min[XX], 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD ) );
  MPL_CHECK_RESULT( MPI_Allreduce( MPI_IN_PLACE, &min[YY], 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD ) );
  MPL_CHECK_RESULT( MPI_Allreduce( MPI_IN_PLACE, &max[XX], 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD ) );
  MPL_CHECK_RESULT( MPI_Allreduce( MPI_IN_PLACE, &max[YY], 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD ) );
}

struct Node
{
  Node() {}
  Node(int gid, int idx)
  {
    g = gid;
    i = idx;
  }
  int g,i;
  bool operator < (const Node& other) const
  {
    return ( g < other.g );
  }
};

inline double sqr(double a) { return a*a; }

}

void build_centroids( FunctionSpace& func_space, ArrayView<double,2>& coords);
void add_median_dual_volume_contribution(
    FunctionSpace& elements,
    FunctionSpace& edges,
    FunctionSpace& nodes,
    ArrayView<double,1>& dual_volumes );
void add_median_dual_volume_contribution(
    FunctionSpace& edges,
    FunctionSpace& nodes,
    ArrayView<double,1>& dual_volumes );
void add_centroid_dual_volume_contribution(
    Mesh& mesh,
    ArrayView<double,1>& dual_volumes );
void build_dual_normals( Mesh& mesh );
void build_skewness(Mesh& mesh );
void make_dual_normals_outward( Mesh& mesh );


void build_median_dual_mesh( Mesh& mesh )
{
  FunctionSpace& nodes   = mesh.function_space( "nodes" );
  ArrayView<double,2> coords        ( nodes.field( "coordinates"    ) );
  ArrayView<double,1> dual_volumes  ( nodes.create_field<double>( "dual_volumes", 1 ) );

  FunctionSpace& quads       = mesh.function_space( "quads" );
  FunctionSpace& triags      = mesh.function_space( "triags" );
  FunctionSpace& edges       = mesh.function_space( "edges" );

  build_centroids(quads,  coords);
  build_centroids(triags, coords);
  build_centroids(edges,  coords);

  add_median_dual_volume_contribution( quads,  edges, nodes, dual_volumes );
  add_median_dual_volume_contribution( triags, edges, nodes, dual_volumes );
  add_median_dual_volume_contribution( edges,  nodes, dual_volumes );

  build_dual_normals( mesh );

  ArrayView<double,1> skewness      ( edges.create_field<double>("skewness",1) );
  ArrayView<double,1> alpha         ( edges.create_field<double>("alpha",1) );
  skewness = 0.;
  alpha = 0.5;


  nodes.parallelise();
  nodes.halo_exchange().execute(dual_volumes);

  ArrayView<double,2> dual_normals  ( edges.field( "dual_normals" ) );
  edges.parallelise();
  edges.halo_exchange().execute(dual_normals);
  make_dual_normals_outward(mesh);

}

void build_centroid_dual_mesh( Mesh& mesh )
{
  FunctionSpace& nodes   = mesh.function_space( "nodes" );
  ArrayView<double,2> coords        ( nodes.field( "coordinates"    ) );
  ArrayView<double,1> dual_volumes  ( nodes.create_field<double>( "dual_volumes", 1 ) );

  FunctionSpace& quads       = mesh.function_space( "quads" );
  FunctionSpace& triags      = mesh.function_space( "triags" );
  FunctionSpace& edges       = mesh.function_space( "edges" );

  build_centroids(quads,  coords);
  build_centroids(triags, coords);
  build_centroids(edges,  coords);

  add_centroid_dual_volume_contribution( mesh, dual_volumes );

  build_dual_normals( mesh );

  build_skewness( mesh );

  nodes.parallelise();
  nodes.halo_exchange().execute(dual_volumes);

  ArrayView<double,2> dual_normals  ( edges.field( "dual_normals" ) );
  edges.parallelise();
  edges.halo_exchange().execute(dual_normals);
}



void build_centroids( FunctionSpace& func_space, ArrayView<double,2>& coords)
{
  if( !func_space.has_field("centroids") )
  {
    int nb_elems = func_space.shape(0);
    IndexView<int,2> elem_nodes( func_space.field( "nodes" ) );
    int nb_nodes_per_elem = elem_nodes.shape(1);
    ArrayView<int,1> elem_glb_idx( func_space.field( "glb_idx" ) );
    ArrayView<double,2> elem_centroids( func_space.create_field<double>( "centroids", 2 ) );
    for (int e=0; e<nb_elems; ++e)
    {
      elem_centroids(e,XX) = 0.;
      elem_centroids(e,YY) = 0.;
      for (int n=0; n<nb_nodes_per_elem; ++n)
      {
        elem_centroids(e,XX) += coords( elem_nodes(e,n), XX );
        elem_centroids(e,YY) += coords( elem_nodes(e,n), YY );
      }
      elem_centroids(e,XX) /= static_cast<double>(nb_nodes_per_elem);
      elem_centroids(e,YY) /= static_cast<double>(nb_nodes_per_elem);
    }
  }
}

void add_median_dual_volume_contribution(
    FunctionSpace& elements,
    FunctionSpace& edges,
    FunctionSpace& nodes,
    ArrayView<double,1>& dual_volumes )
{
  int nb_elems = elements.shape(0);
  ArrayView<double,2> elem_centroids ( elements.field("centroids") );
  IndexView<int,   2> elem_to_edges  ( elements.field("to_edge") );
  ArrayView<double,2> edge_centroids ( edges.field("centroids") );
  IndexView<int,   2> edge_nodes     ( edges.field("nodes") );
  ArrayView<double,2> node_coords    ( nodes.field("coordinates") );
  ArrayView<int,   1> elem_glb_idx   ( elements.field("glb_idx") );
  ArrayView<int,   1> edge_glb_idx   ( edges.field("glb_idx") );
  int nb_edges_per_elem = elem_to_edges.shape(1);


  // special ordering for bit-identical results
  std::vector<Node> ordering(nb_elems);
  for (int elem=0; elem<nb_elems; ++elem)
  {
    ordering[elem] = Node( LatLonPoint(elem_centroids[elem]).uid(), elem );
  }
  std::sort( ordering.data(), ordering.data()+nb_elems );


  for (int jelem=0; jelem<nb_elems; ++jelem)
  {
    int elem = ordering[jelem].i;
    double x0 = elem_centroids(elem,XX);
    double y0 = elem_centroids(elem,YY);
    for (int jedge=0; jedge<nb_edges_per_elem; ++jedge)
    {
      int edge = elem_to_edges(elem,jedge);
      double x1 = edge_centroids(edge,XX);
      double y1 = edge_centroids(edge,YY);
      for( int j=0; j<2; ++j )
      {
        int node = edge_nodes(edge,j);
        double x2 = node_coords(node,XX);
        double y2 = node_coords(node,YY);
        double triag_area = std::abs( x0*(y1-y2)+x1*(y2-y0)+x2*(y0-y1) )*0.5;
        dual_volumes(node) += triag_area;
      }
    }
  }
}

void add_median_dual_volume_contribution(
    FunctionSpace& edges,
    FunctionSpace& nodes,
    ArrayView<double,1>& dual_volumes )
{
  ArrayView<int,   1> node_glb_idx  ( nodes.field("glb_idx"    ) );
  ArrayView<double,2> edge_centroids( edges.field("centroids"  ) );
  IndexView<int,   2> edge_nodes    ( edges.field("nodes"      ) );
  ArrayView<int,   1> edge_glb_idx  ( edges.field("glb_idx"    ) );
  IndexView<int,   2> edge_to_elem  ( edges.field("to_elem"    ) );
  ArrayView<double,2> node_coords   ( nodes.field("coordinates") );
  int nb_edges = edges.shape(0);
  std::map<int,std::vector<int> > node_to_bdry_edge;
  for(int edge=0; edge<nb_edges; ++edge)
  {
    if ( edge_to_elem(edge,0) >= 0 && edge_to_elem(edge,3) < 0)
    {
      node_to_bdry_edge[ edge_nodes(edge,0) ].push_back(edge);
      node_to_bdry_edge[ edge_nodes(edge,1) ].push_back(edge);
    }
  }

  double tol = 1.e-6;
  double min[2], max[2];
  global_bounding_box( nodes, min, max );

  std::set<int> visited_nodes;

  std::map<int,std::vector<int> >::iterator it;
  for( it=node_to_bdry_edge.begin(); it!=node_to_bdry_edge.end(); ++it)
  {
    int node = (*it).first;
    std::vector<int>& bdry_edges = (*it).second;
    double x0 = node_coords(node,XX);
    double y0 = node_coords(node,YY);
    double x1,y1, x2,y2;
    for (int jedge=0; jedge<bdry_edges.size(); ++jedge)
    {
      int edge = bdry_edges[jedge];
      x1 = edge_centroids(edge,XX);
      y1 = edge_centroids(edge,YY);
      x2 = x1;
      y2 = 0.;
      if ( std::abs(y1-max[YY])<tol )
        y2 = M_PI_2;
      else if ( std::abs(y1-min[YY])<tol )
        y2 = -M_PI_2;

      if( y2!=0 )
      {
        //std::cout << "edge " << edge_glb_idx(edge) << " adding contribution for node " << node_glb_idx(node) << std::endl;
        double quad_area = std::abs( (x1-x0)*(y2-y0) );
        dual_volumes(node) += quad_area;
      }
    }
  }
}

void add_centroid_dual_volume_contribution(
    Mesh& mesh,
    ArrayView<double,1>& dual_volumes )
{
  FunctionSpace& nodes = mesh.function_space("nodes");
  FunctionSpace& edges = mesh.function_space("edges");
  ArrayView<int,   1> node_glb_idx  ( nodes.field("glb_idx"    ) );
  ArrayView<double,2> edge_centroids( edges.field("centroids"  ) );
  IndexView<int,   2> edge_nodes    ( edges.field("nodes"      ) );
  ArrayView<int,   1> edge_glb_idx  ( edges.field("glb_idx"    ) );
  IndexView<int,   2> edge_to_elem  ( edges.field("to_elem"    ) );
  ArrayView<double,2> node_coords   ( nodes.field("coordinates") );
  std::vector< ArrayView<double,2> > elem_centroids(mesh.nb_function_spaces());
  for( int f=0; f<mesh.nb_function_spaces(); ++f )
  {
    FunctionSpace& elements = mesh.function_space(f);
    if( elements.metadata().get<int>("type") == Entity::ELEMS )
    {
      elem_centroids[f] = ArrayView<double,2>( elements.field("centroids") );
    }
  }
  double tol = 1.e-6;
  double min[2], max[2];
  global_bounding_box( nodes, min, max );

  int nb_edges = edges.shape(0);

  // special ordering for bit-identical results
  std::vector<Node> ordering(nb_edges);
  for (int edge=0; edge<nb_edges; ++edge)
  {
    ordering[edge] = Node( LatLonPoint(edge_centroids[edge]).uid(), edge );
  }
  std::sort( ordering.data(), ordering.data()+nb_edges );


  for(int jedge=0; jedge<nb_edges; ++jedge)
  {
    int edge = ordering[jedge].i;
    if ( edge_to_elem(edge,0) >= 0 && edge_to_elem(edge,2) >= 0 )
    {
      double x0 = elem_centroids[edge_to_elem(edge,0)](edge_to_elem(edge,1),XX);
      double y0 = elem_centroids[edge_to_elem(edge,0)](edge_to_elem(edge,1),YY);
      double x1 = elem_centroids[edge_to_elem(edge,2)](edge_to_elem(edge,3),XX);
      double y1 = elem_centroids[edge_to_elem(edge,2)](edge_to_elem(edge,3),YY);
      for( int jnode=0; jnode<2; ++jnode )
      {
        int node = edge_nodes(edge,jnode);
        double x2 = node_coords( node, XX );
        double y2 = node_coords( node, YY );
        double triag_area = std::abs( x0*(y1-y2)+x1*(y2-y0)+x2*(y0-y1) )*0.5;
        dual_volumes(node) += triag_area;
      }
    }
    else if ( edge_to_elem(edge,0) >= 0 && edge_to_elem(edge,2) < 0  )
    {
      // This is a boundary edge
      double x0 = elem_centroids[edge_to_elem(edge,0)](edge_to_elem(edge,1),XX);
      double y0 = elem_centroids[edge_to_elem(edge,0)](edge_to_elem(edge,1),YY);
      double x1 = x0;
      double y1 = 0;
      double y_edge = edge_centroids(edge,YY);
      if ( std::abs(y_edge-max[YY])<tol )
        y1 = M_PI_2;
      else if ( std::abs(y_edge-min[YY])<tol )
        y1 = -M_PI_2;

      if( y1 != 0. )
      {
        for( int jnode=0; jnode<2; ++jnode )
        {
          int node = edge_nodes(edge,jnode);
          double x2 = node_coords( node, XX );
          double y2 = node_coords( node, YY );
          double triag_area = std::abs( x0*(y1-y2)+x1*(y2-y0)+x2*(y0-y1) )*0.5;
          dual_volumes(node) += triag_area;
          double x3 = x2;
          double y3 = y1;
          triag_area = std::abs( x3*(y1-y2)+x1*(y2-y3)+x2*(y3-y1) )*0.5;
          dual_volumes(node) += triag_area;
        }
      }
    }
  }
}


void build_dual_normals( Mesh& mesh )
{
  std::vector< ArrayView<double,2> > elem_centroids( mesh.nb_function_spaces() );
  for (int func_space_idx=0; func_space_idx<mesh.nb_function_spaces(); ++func_space_idx)
  {
    FunctionSpace& func_space = mesh.function_space(func_space_idx);
    if( func_space.has_field("centroids") )
      elem_centroids[func_space_idx] = ArrayView<double,2>( func_space.field<double>("centroids") );
  }

  FunctionSpace&  nodes = mesh.function_space("nodes");
  ArrayView<double,2> node_coords( nodes.field("coordinates") );
  double min[2], max[2];
  global_bounding_box( nodes, min, max );
  double tol = 1.e-6;

  double xl, yl, xr, yr, dx, dy;
  FunctionSpace&  edges = mesh.function_space("edges");
  IndexView<int,   2> edge_to_elem  ( edges.field("to_elem"  ) );
  IndexView<int,   2> edge_nodes    ( edges.field("nodes"    ) );
  ArrayView<double,2> edge_centroids( edges.field("centroids") );
  ArrayView<double,2> dual_normals  ( edges.create_field<double>("dual_normals",2) );
  int nb_edges = edges.shape(0);

  std::map<int,std::vector<int> > node_to_bdry_edge;
  for(int edge=0; edge<nb_edges; ++edge)
  {
    if ( edge_to_elem(edge,0) >= 0 && edge_to_elem(edge,3) < 0)
    {
      node_to_bdry_edge[ edge_nodes(edge,0) ].push_back(edge);
      node_to_bdry_edge[ edge_nodes(edge,1) ].push_back(edge);
    }
  }

  for (int edge=0; edge<nb_edges; ++edge)
  {
    if( edge_to_elem(edge,0) < 0 )
    {
      // this is a pole edge
      // only compute for one node
      for (int n=0; n<2; ++n)
      {
        int node = edge_nodes(edge,n);
        std::vector<int>& bdry_edges = node_to_bdry_edge[node];
        double x[2];
        int cnt=0;
        for (int jedge=0; jedge<bdry_edges.size(); ++jedge)
        {
          int bdry_edge = bdry_edges[jedge];
          if ( std::abs(edge_centroids(bdry_edge,YY)-max[YY])<tol )
          {
            edge_centroids(edge,YY) = M_PI_2;
            x[cnt] = edge_centroids(bdry_edge,XX);
            ++cnt;
          }
          else if ( std::abs(edge_centroids(bdry_edge,YY)-min[YY])<tol )
          {
            edge_centroids(edge,YY) = -M_PI_2;
            x[cnt] = edge_centroids(bdry_edge,XX);
            ++cnt;
          }
        }
        if (cnt == 2 )
        {
          dual_normals(edge,XX) = 0;
          if (node_coords(node,YY) < 0.)
            dual_normals(edge,YY) = -std::abs(x[1]-x[0]);
          else if (node_coords(node,YY) > 0.)
            dual_normals(edge,YY) = std::abs(x[1]-x[0]);

          //std::cout << "pole dual_normal = " << dual_normals(YY,edge) << std::endl;
          break;
        }
      }
    }
    else
    {
      int left_func_space_idx  = edge_to_elem(edge,0);
      int left_elem            = edge_to_elem(edge,1);
      int right_func_space_idx = edge_to_elem(edge,2);
      int right_elem           = edge_to_elem(edge,3);
      xl = elem_centroids[left_func_space_idx](left_elem,XX);
      yl = elem_centroids[left_func_space_idx](left_elem,YY);
      if( right_elem < 0 )
      {
        xr = edge_centroids(edge,XX);
        yr = edge_centroids(edge,YY);;
        if ( std::abs(yr-max[YY])<tol )
          yr = M_PI_2;
        else if( std::abs(yr-min[YY])<tol )
          yr = -M_PI_2;
      }
      else
      {
        xr = elem_centroids[right_func_space_idx](right_elem,XX);
        yr = elem_centroids[right_func_space_idx](right_elem,YY);
      }

      dual_normals(edge,XX) =  yl-yr;
      dual_normals(edge,YY) = -xl+xr;
    }
  }
}

void make_dual_normals_outward( Mesh& mesh )
{

  FunctionSpace&  nodes = mesh.function_space("nodes");
  ArrayView<double,2> node_coords( nodes.field("coordinates") );

  FunctionSpace&  edges = mesh.function_space("edges");
  IndexView<int,   2> edge_to_elem  ( edges.field("to_elem"  ) );
  IndexView<int,   2> edge_nodes    ( edges.field("nodes"    ) );
  ArrayView<double,2> dual_normals  ( edges.field("dual_normals") );
  int nb_edges = edges.shape(0);

  for (int edge=0; edge<nb_edges; ++edge)
  {
    if( edge_to_elem(edge,0) < 0 )
    {

    }
    else
    {
      // Make normal point from node 1 to node 2
      double dx = node_coords( edge_nodes(edge,1), XX ) - node_coords( edge_nodes(edge,0), XX );
      double dy = node_coords( edge_nodes(edge,1), YY ) - node_coords( edge_nodes(edge,0), YY );
      if( dx*dual_normals(edge,XX) + dy*dual_normals(edge,YY) < 0 )
      {
        dual_normals(edge,XX) = - dual_normals(edge,XX);
        dual_normals(edge,YY) = - dual_normals(edge,YY);
      }
    }
  }
}



void build_skewness( Mesh& mesh )
{
  std::vector< ArrayView<double,2> > elem_centroids( mesh.nb_function_spaces() );
  for (int func_space_idx=0; func_space_idx<mesh.nb_function_spaces(); ++func_space_idx)
  {
    FunctionSpace& func_space = mesh.function_space(func_space_idx);
    if( func_space.has_field("centroids") )
      elem_centroids[func_space_idx] = ArrayView<double,2>( func_space.field<double>("centroids") );
  }

  FunctionSpace&  nodes = mesh.function_space("nodes");
  ArrayView<double,2> node_coords( nodes.field<double>("coordinates") );
  double min[2], max[2];
  global_bounding_box( nodes, min, max );
  double pi = acos(-1.);
  double tol = 1.e-6;

  double x1, y1, x2, y2, xc1, yc1, xc2, yc2, xi, yi;
  FunctionSpace&  edges = mesh.function_space("edges");
  IndexView<int   ,2> edge_to_elem  ( edges.field("to_elem"  ) );
  IndexView<int   ,2> edge_nodes    ( edges.field("nodes"    ) );
  ArrayView<double,2> edge_centroids( edges.field("centroids") );
  ArrayView<double,1> skewness      ( edges.create_field<double>("skewness",1) );
  ArrayView<double,1> alpha         ( edges.create_field<double>("alpha",1) );
  int nb_edges = edges.shape(0);

  // special ordering for bit-identical results
  std::vector<Node> ordering(nb_edges);
  for (int edge=0; edge<nb_edges; ++edge)
  {
    ordering[edge] = Node( LatLonPoint(edge_centroids[edge]).uid(), edge );
  }
  std::sort( ordering.data(), ordering.data()+nb_edges );


  for(int jedge=0; jedge<nb_edges; ++jedge)
  {
    int edge = ordering[jedge].i;
    if( edge_to_elem(edge,0) < 0 )
    {
      // this is a pole edge
      // only compute for one node
      skewness(edge) = 0.;
    }
    else
    {
      int left_func_space_idx  = edge_to_elem(edge,0);
      int left_elem            = edge_to_elem(edge,1);
      int right_func_space_idx = edge_to_elem(edge,2);
      int right_elem           = edge_to_elem(edge,3);
      xc1 = elem_centroids[left_func_space_idx](left_elem,XX);
      yc1 = elem_centroids[left_func_space_idx](left_elem,YY);
      if( right_elem < 0 )
      {
        xc2 = edge_centroids(edge,XX);
        yc2 = edge_centroids(edge,YY);
        if ( std::abs(yc2-max[YY])<tol )
          yc2 = M_PI_2;
        else if( std::abs(yc2-min[YY])<tol )
          yc2 = -M_PI_2;
      }
      else
      {
        xc2 = elem_centroids[right_func_space_idx](right_elem,XX);
        yc2 = elem_centroids[right_func_space_idx](right_elem,YY);
      }

      x1 = node_coords(edge_nodes(edge,0),XX);
      y1 = node_coords(edge_nodes(edge,0),YY);
      x2 = node_coords(edge_nodes(edge,1),XX);
      y2 = node_coords(edge_nodes(edge,1),YY);

      xi = ( x1*(xc2*(-y2 + yc1) + xc1*(y2 - yc2))
             + x2*(xc2*(y1 - yc1) + xc1*(-y1 + yc2)) ) /
           (-((xc1 - xc2)*(y1 - y2)) + (x1 - x2)*(yc1 - yc2));
      yi = ( xc2*(y1 - y2)*yc1 + x1*y2*yc1
             - xc1*y1*yc2 - x1*y2*yc2 + xc1*y2*yc2
             + x2*y1*(-yc1 + yc2))/
             (-((xc1 - xc2)*(y1 - y2)) + (x1 - x2)*(yc1 - yc2));
      double r1 = 0;
      double r2 = sqrt( sqr(x2-x1) + sqr(y2-y1) );
      double rs = sqrt( sqr(xi-x1) + sqr(yi-y1) );
      skewness(edge) = (r1-2.*rs+r2)/(r2-r1);
      alpha(edge) = 0.5*(skewness(edge)+1.);
    }
  }
}


// ------------------------------------------------------------------
// C wrapper interfaces to C++ routines

void atlas__build_median_dual_mesh ( Mesh* mesh) {
  build_median_dual_mesh(*mesh);
}
void atlas__build_centroid_dual_mesh ( Mesh* mesh) {
  build_centroid_dual_mesh(*mesh);
}

// ------------------------------------------------------------------

} // namespace actions
} // namespace atlas

