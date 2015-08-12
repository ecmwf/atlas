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

#include "atlas/atlas_config.h"
#include "atlas/runtime/ErrorHandling.h"
#include "atlas/mpl/Checksum.h"
#include "atlas/Mesh.h"
#include "atlas/FunctionSpace.h"
#include "atlas/Field.h"
#include "atlas/actions/BuildDualMesh.h"
#include "atlas/Parameters.h"
#include "atlas/util/ArrayView.h"
#include "atlas/util/IndexView.h"
#include "atlas/util/Unique.h"
#include "atlas/grids/ReducedGrid.h"

namespace atlas {
namespace actions {

namespace {

void global_bounding_box( FunctionSpace& nodes, double min[2], double max[2] )
{
  ArrayView<double,2> lonlat( nodes.field("lonlat") );
  const int nb_nodes = nodes.shape(0);
  min[LON] =  std::numeric_limits<double>::max();
  min[LAT] =  std::numeric_limits<double>::max();
  max[LON] = -std::numeric_limits<double>::max();
  max[LAT] = -std::numeric_limits<double>::max();
  for (int node=0; node<nb_nodes; ++node)
  {
    min[LON] = std::min( min[LON], lonlat(node,LON) );
    min[LAT] = std::min( min[LAT], lonlat(node,LAT) );
    max[LON] = std::max( max[LON], lonlat(node,LON) );
    max[LAT] = std::max( max[LAT], lonlat(node,LAT) );
  }
  ECKIT_MPI_CHECK_RESULT( MPI_Allreduce( MPI_IN_PLACE, &min[LON], 1, MPI_DOUBLE, MPI_MIN, eckit::mpi::comm() ) );
  ECKIT_MPI_CHECK_RESULT( MPI_Allreduce( MPI_IN_PLACE, &min[LAT], 1, MPI_DOUBLE, MPI_MIN, eckit::mpi::comm() ) );
  ECKIT_MPI_CHECK_RESULT( MPI_Allreduce( MPI_IN_PLACE, &max[LON], 1, MPI_DOUBLE, MPI_MAX, eckit::mpi::comm() ) );
  ECKIT_MPI_CHECK_RESULT( MPI_Allreduce( MPI_IN_PLACE, &max[LAT], 1, MPI_DOUBLE, MPI_MAX, eckit::mpi::comm() ) );
}

struct Node
{
  Node() {}
  Node(gidx_t gid, int idx)
  {
    g = gid;
    i = idx;
  }
  gidx_t g;
  int i;

  bool operator < (const Node& other) const
  {
    return ( g < other.g );
  }
};

inline double sqr(double a) { return a*a; }

}

void build_centroids( FunctionSpace& func_space, ArrayView<double,2>& lonlat);
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
  ArrayView<double,2> lonlat        ( nodes.field( "lonlat"    ) );
  ArrayView<double,1> dual_volumes  ( nodes.create_field<double>( "dual_volumes", 1 ) );

  FunctionSpace& quads       = mesh.function_space( "quads" );
  FunctionSpace& triags      = mesh.function_space( "triags" );
  FunctionSpace& edges       = mesh.function_space( "edges" );

  build_centroids(quads,  lonlat);
  build_centroids(triags, lonlat);
  build_centroids(edges,  lonlat);

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
  ArrayView<double,2> lonlat        ( nodes.field( "lonlat"    ) );
  ArrayView<double,1> dual_volumes  ( nodes.create_field<double>( "dual_volumes", 1 ) );

  FunctionSpace& quads       = mesh.function_space( "quads" );
  FunctionSpace& triags      = mesh.function_space( "triags" );
  FunctionSpace& edges       = mesh.function_space( "edges" );

  build_centroids(quads,  lonlat);
  build_centroids(triags, lonlat);
  build_centroids(edges,  lonlat);

  add_centroid_dual_volume_contribution( mesh, dual_volumes );

  build_dual_normals( mesh );

  build_skewness( mesh );

  nodes.parallelise();
  nodes.halo_exchange().execute(dual_volumes);

  ArrayView<double,2> dual_normals  ( edges.field( "dual_normals" ) );
  edges.parallelise();
  edges.halo_exchange().execute(dual_normals);
}



void build_centroids( FunctionSpace& func_space, ArrayView<double,2>& lonlat)
{
  if( !func_space.has_field("centroids") )
  {
    int nb_elems = func_space.shape(0);
    IndexView<int,2> elem_nodes( func_space.field( "nodes" ) );
    int nb_nodes_per_elem = elem_nodes.shape(1);
    ArrayView<gidx_t,1> elem_glb_idx( func_space.field( "glb_idx" ) );
    ArrayView<double,2> elem_centroids( func_space.create_field<double>( "centroids", 2 ) );
    for (int e=0; e<nb_elems; ++e)
    {
      elem_centroids(e,LON) = 0.;
      elem_centroids(e,LAT) = 0.;
      for (int n=0; n<nb_nodes_per_elem; ++n)
      {
        elem_centroids(e,LON) += lonlat( elem_nodes(e,n), LON );
        elem_centroids(e,LAT) += lonlat( elem_nodes(e,n), LAT );
      }
      elem_centroids(e,LON) /= static_cast<double>(nb_nodes_per_elem);
      elem_centroids(e,LAT) /= static_cast<double>(nb_nodes_per_elem);
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
  ArrayView<double,2> node_lonlat    ( nodes.field("lonlat") );
  ArrayView<gidx_t,1> elem_glb_idx   ( elements.field("glb_idx") );
  ArrayView<gidx_t,1> edge_glb_idx   ( edges.field("glb_idx") );
  int nb_edges_per_elem = elem_to_edges.shape(1);


  // special ordering for bit-identical results
  std::vector<Node> ordering(nb_elems);
  for (int elem=0; elem<nb_elems; ++elem)
  {
    ordering[elem] = Node( util::unique_lonlat(elem_centroids[elem]), elem );
  }
  std::sort( ordering.data(), ordering.data()+nb_elems );


  for (int jelem=0; jelem<nb_elems; ++jelem)
  {
    int elem = ordering[jelem].i;
    double x0 = elem_centroids(elem,LON);
    double y0 = elem_centroids(elem,LAT);
    for (int jedge=0; jedge<nb_edges_per_elem; ++jedge)
    {
      int edge = elem_to_edges(elem,jedge);
      double x1 = edge_centroids(edge,LON);
      double y1 = edge_centroids(edge,LAT);
      for( int j=0; j<2; ++j )
      {
        int node = edge_nodes(edge,j);
        double x2 = node_lonlat(node,LON);
        double y2 = node_lonlat(node,LAT);
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
  ArrayView<gidx_t,1> node_glb_idx  ( nodes.field("glb_idx"    ) );
  ArrayView<double,2> edge_centroids( edges.field("centroids"  ) );
  IndexView<int,   2> edge_nodes    ( edges.field("nodes"      ) );
  ArrayView<gidx_t,1> edge_glb_idx  ( edges.field("glb_idx"    ) );
  IndexView<int,   2> edge_to_elem  ( edges.field("to_elem"    ) );
  ArrayView<double,2> node_lonlat   ( nodes.field("lonlat") );
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
    double x0 = node_lonlat(node,LON);
    double y0 = node_lonlat(node,LAT);
    double x1, y1, y2;
    for (size_t jedge = 0; jedge < bdry_edges.size(); ++jedge)
    {
      int edge = bdry_edges[jedge];
      x1 = edge_centroids(edge,LON);
      y1 = edge_centroids(edge,LAT);

      y2 = 0.;
      if ( std::abs(y1-max[LAT])<tol )
        y2 = 90.;
      else if ( std::abs(y1-min[LAT])<tol )
        y2 = -90.;

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
  ArrayView<gidx_t,1> node_glb_idx  ( nodes.field("glb_idx"    ) );
  ArrayView<double,2> edge_centroids( edges.field("centroids"  ) );
  IndexView<int,   2> edge_nodes    ( edges.field("nodes"      ) );
  ArrayView<gidx_t,1> edge_glb_idx  ( edges.field("glb_idx"    ) );
  IndexView<int,   2> edge_to_elem  ( edges.field("to_elem"    ) );
  ArrayView<double,2> node_lonlat   ( nodes.field("lonlat") );
  std::vector< ArrayView<double,2> > elem_centroids(mesh.nb_function_spaces());
  for(size_t f = 0; f < mesh.nb_function_spaces(); ++f)
  {
    FunctionSpace& elements = mesh.function_space(f);
    if( elements.metadata().get<long>("type") == Entity::ELEMS )
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
    ordering[edge] = Node( util::unique_lonlat(edge_centroids[edge]), edge );
  }
  std::sort( ordering.data(), ordering.data()+nb_edges );


  for(int jedge=0; jedge<nb_edges; ++jedge)
  {
    int edge = ordering[jedge].i;
    if ( edge_to_elem(edge,0) >= 0 && edge_to_elem(edge,2) >= 0 )
    {
      double x0 = elem_centroids[edge_to_elem(edge,0)](edge_to_elem(edge,1),LON);
      double y0 = elem_centroids[edge_to_elem(edge,0)](edge_to_elem(edge,1),LAT);
      double x1 = elem_centroids[edge_to_elem(edge,2)](edge_to_elem(edge,3),LON);
      double y1 = elem_centroids[edge_to_elem(edge,2)](edge_to_elem(edge,3),LAT);
      for( int jnode=0; jnode<2; ++jnode )
      {
        int node = edge_nodes(edge,jnode);
        double x2 = node_lonlat( node, LON );
        double y2 = node_lonlat( node, LAT );
        double triag_area = std::abs( x0*(y1-y2)+x1*(y2-y0)+x2*(y0-y1) )*0.5;
        dual_volumes(node) += triag_area;
      }
    }
    else if ( edge_to_elem(edge,0) >= 0 && edge_to_elem(edge,2) < 0  )
    {
      // This is a boundary edge
      double x0 = elem_centroids[edge_to_elem(edge,0)](edge_to_elem(edge,1),LON);
      double y0 = elem_centroids[edge_to_elem(edge,0)](edge_to_elem(edge,1),LAT);
      double x1 = x0;
      double y1 = 0;
      double y_edge = edge_centroids(edge,LAT);
      if ( std::abs(y_edge-max[LAT])<tol )
        y1 = 90.;
      else if ( std::abs(y_edge-min[LAT])<tol )
        y1 = -90.;

      if( y1 != 0. )
      {
        for( int jnode=0; jnode<2; ++jnode )
        {
          int node = edge_nodes(edge,jnode);
          double x2 = node_lonlat( node, LON );
          double y2 = node_lonlat( node, LAT );
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
  for (size_t func_space_idx = 0; func_space_idx < mesh.nb_function_spaces(); ++func_space_idx)
  {
    FunctionSpace& func_space = mesh.function_space(func_space_idx);
    if( func_space.has_field("centroids") )
      elem_centroids[func_space_idx] = ArrayView<double,2>( func_space.field("centroids") );
  }

  FunctionSpace&  nodes = mesh.function_space("nodes");
  ArrayView<double,2> node_lonlat( nodes.field("lonlat") );
  double min[2], max[2];
  global_bounding_box( nodes, min, max );
  double tol = 1.e-6;

  double xl, yl, xr, yr;
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
        for (size_t jedge = 0; jedge < bdry_edges.size(); ++jedge)
        {
          int bdry_edge = bdry_edges[jedge];
          if ( std::abs(edge_centroids(bdry_edge,LAT)-max[LAT])<tol )
          {
            edge_centroids(edge,LAT) = 90.;
            x[cnt] = edge_centroids(bdry_edge,LON);
            ++cnt;
          }
          else if ( std::abs(edge_centroids(bdry_edge,LAT)-min[LAT])<tol )
          {
            edge_centroids(edge,LAT) = -90.;
            x[cnt] = edge_centroids(bdry_edge,LON);
            ++cnt;
          }
        }
        if (cnt == 2 )
        {
          dual_normals(edge,LON) = 0;
          if (node_lonlat(node,LAT) < 0.)
            dual_normals(edge,LAT) = -std::abs(x[1]-x[0]);
          else if (node_lonlat(node,LAT) > 0.)
            dual_normals(edge,LAT) = std::abs(x[1]-x[0]);

          //std::cout << "pole dual_normal = " << dual_normals(LAT,edge) << std::endl;
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
      xl = elem_centroids[left_func_space_idx](left_elem,LON);
      yl = elem_centroids[left_func_space_idx](left_elem,LAT);
      if( right_elem < 0 )
      {
        xr = edge_centroids(edge,LON);
        yr = edge_centroids(edge,LAT);;
        if ( std::abs(yr-max[LAT])<tol )
          yr = 90.;
        else if( std::abs(yr-min[LAT])<tol )
          yr = -90.;
      }
      else
      {
        xr = elem_centroids[right_func_space_idx](right_elem,LON);
        yr = elem_centroids[right_func_space_idx](right_elem,LAT);
      }

      dual_normals(edge,LON) =  yl-yr;
      dual_normals(edge,LAT) = -xl+xr;
    }
  }
}

void make_dual_normals_outward( Mesh& mesh )
{

  FunctionSpace&  nodes = mesh.function_space("nodes");
  ArrayView<double,2> node_lonlat( nodes.field("lonlat") );

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
      double dx = node_lonlat( edge_nodes(edge,1), LON ) - node_lonlat( edge_nodes(edge,0), LON );
      double dy = node_lonlat( edge_nodes(edge,1), LAT ) - node_lonlat( edge_nodes(edge,0), LAT );
      if( dx*dual_normals(edge,LON) + dy*dual_normals(edge,LAT) < 0 )
      {
        dual_normals(edge,LON) = - dual_normals(edge,LON);
        dual_normals(edge,LAT) = - dual_normals(edge,LAT);
      }
    }
  }
}



void build_skewness( Mesh& mesh )
{
  std::vector< ArrayView<double,2> > elem_centroids( mesh.nb_function_spaces() );
  for (size_t func_space_idx = 0; func_space_idx < mesh.nb_function_spaces(); ++func_space_idx)
  {
    FunctionSpace& func_space = mesh.function_space(func_space_idx);
    if( func_space.has_field("centroids") )
      elem_centroids[func_space_idx] = ArrayView<double,2>( func_space.field("centroids") );
  }

  FunctionSpace&  nodes = mesh.function_space("nodes");
  ArrayView<double,2> node_lonlat( nodes.field("lonlat") );
  double min[2], max[2];
  global_bounding_box( nodes, min, max );
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
    ordering[edge] = Node( util::unique_lonlat(edge_centroids[edge]), edge );
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
      xc1 = elem_centroids[left_func_space_idx](left_elem,LON);
      yc1 = elem_centroids[left_func_space_idx](left_elem,LAT);
      if( right_elem < 0 )
      {
        xc2 = edge_centroids(edge,LON);
        yc2 = edge_centroids(edge,LAT);
        if ( std::abs(yc2-max[LAT])<tol )
          yc2 = 90.;
        else if( std::abs(yc2-min[LAT])<tol )
          yc2 = -90.;
      }
      else
      {
        xc2 = elem_centroids[right_func_space_idx](right_elem,LON);
        yc2 = elem_centroids[right_func_space_idx](right_elem,LAT);
      }

      x1 = node_lonlat(edge_nodes(edge,0),LON);
      y1 = node_lonlat(edge_nodes(edge,0),LAT);
      x2 = node_lonlat(edge_nodes(edge,1),LON);
      y2 = node_lonlat(edge_nodes(edge,1),LAT);

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




void build_brick_dual_mesh( Mesh& mesh )
{
  if( const grids::ReducedGrid* g = dynamic_cast<const grids::ReducedGrid*>(&mesh.grid()) )
  {
    if( eckit::mpi::size() != 1 )
      throw eckit::UserError("Cannot build_brick_dual_mesh with more than 1 task",Here());

    FunctionSpace& nodes   = mesh.function_space( "nodes" );
    ArrayView<double,2> lonlat        ( nodes.field( "lonlat"    ) );
    ArrayView<double,1> dual_volumes  ( nodes.create_field<double>( "dual_volumes", 1 ) );
    ArrayView<gidx_t,1> gidx  ( nodes.field( "glb_idx" ) );

    int c=0;
    int n=0;
    for(size_t jlat = 0; jlat < g->nlat(); ++jlat)
    {
      double lat = g->lat(jlat);
      double latN = (jlat==0) ? 90. : 0.5*(lat+g->lat(jlat-1));
      double latS = (jlat==g->nlat()-1) ? -90. : 0.5*(lat+g->lat(jlat+1));
      double dlat = (latN-latS);
      double dlon = 360./static_cast<double>(g->nlon(jlat));

      for(size_t jlon = 0; jlon < g->nlon(jlat); ++jlon)
      {
        while( gidx(c) != n+1 ) c++;
        ASSERT( lonlat(c,LON) == g->lon(jlat,jlon) );
        ASSERT( lonlat(c,LAT) == lat );
        dual_volumes(c) = dlon*dlat;
        ++n;
      }

    }

    nodes.parallelise();
    nodes.halo_exchange().execute(dual_volumes);
  }
  else
  {
    throw eckit::BadCast("Cannot build_brick_dual_mesh with mesh provided grid type",Here());
  }
}


// ------------------------------------------------------------------
// C wrapper interfaces to C++ routines

void atlas__build_median_dual_mesh ( Mesh* mesh) {
  ATLAS_ERROR_HANDLING( build_median_dual_mesh(*mesh) );
}
void atlas__build_centroid_dual_mesh ( Mesh* mesh) {
  ATLAS_ERROR_HANDLING( build_centroid_dual_mesh(*mesh) );
}

// ------------------------------------------------------------------

} // namespace actions
} // namespace atlas

