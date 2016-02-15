/*
 * (C) Copyright 1996-2015 ECMWF.
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
#include "atlas/mesh/Nodes.h"
#include "atlas/util/ArrayView.h"
#include "atlas/util/IndexView.h"
#include "atlas/util/Unique.h"
#include "atlas/grids/ReducedGrid.h"
#include "atlas/functionspace/Nodes.h"
#include "atlas/functionspace/Edges.h"
#include "atlas/mesh/HybridElements.h"

using atlas::functionspace::Nodes;
using atlas::functionspace::Halo;

namespace atlas {
namespace actions {

namespace {

void global_bounding_box( const mesh::Nodes& nodes, double min[2], double max[2] )
{
  ArrayView<double,2> lonlat( nodes.lonlat() );
  const int nb_nodes = nodes.size();
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
  eckit::mpi::all_reduce(min,2,eckit::mpi::min());
  eckit::mpi::all_reduce(max,2,eckit::mpi::max());
}

struct Node
{
  Node() {}
  Node(gidx_t gid, idx_t idx)
  {
    g = gid;
    i = idx;
  }
  gidx_t g;
  idx_t i;

  bool operator < (const Node& other) const
  {
    return ( g < other.g );
  }
};

inline double sqr(double a) { return a*a; }

}

Array* build_centroids_lonlat( const mesh::HybridElements& , const Field& lonlat);

void add_centroid_dual_volume_contribution(
    Mesh& mesh,
    ArrayView<double,1>& dual_volumes );
void add_median_dual_volume_contribution_cells(
    const mesh::HybridElements& cells,
    const mesh::HybridElements& edges,
    const mesh::Nodes& nodes,
    Array& array_dual_volumes );
void add_median_dual_volume_contribution_poles(
    const mesh::HybridElements& edges,
    const mesh::Nodes& nodes,
    Array& array_dual_volumes );
void build_dual_normals( Mesh& mesh );
void build_dual_normals_new( Mesh& mesh );
void build_skewness(Mesh& mesh );
void make_dual_normals_outward( Mesh& mesh );
void make_dual_normals_outward_new( Mesh& mesh );

void build_median_dual_mesh_new( Mesh& mesh )
{
  mesh::Nodes& nodes   = mesh.nodes();
  mesh::HybridElements& edges = mesh.edges();
  Field& dual_volumes = nodes.add( Field::create<double>( "dual_volumes", make_shape(nodes.size()) ) );

  if( ! mesh.cells().has_field("centroids_lonlat") )
    mesh.cells().add( Field::create("centroids_lonlat",build_centroids_lonlat(mesh.cells(),mesh.nodes().lonlat())) );

  if( ! mesh.edges().has_field("centroids_lonlat") )
    mesh.edges().add( Field::create("centroids_lonlat",build_centroids_lonlat(mesh.edges(),mesh.nodes().lonlat())) );

  add_median_dual_volume_contribution_cells(mesh.cells(),mesh.edges(),mesh.nodes(),nodes.field("dual_volumes").array());
  add_median_dual_volume_contribution_poles(mesh.edges(),mesh.nodes(),nodes.field("dual_volumes").array());

  build_dual_normals_new( mesh );

  ArrayView<double,1> skewness ( mesh.edges().add( Field::create<double>("skewness",make_shape(mesh.edges().size()))) );
  ArrayView<double,1> alpha    ( mesh.edges().add( Field::create<double>("alpha",make_shape(mesh.edges().size()))) );
  skewness = 0.;
  alpha = 0.5;

  eckit::SharedPtr<functionspace::Nodes> nodes_fs( new functionspace::Nodes(mesh,Halo(mesh)) );
  nodes_fs->haloExchange(nodes.field( "dual_volumes" ));

  eckit::SharedPtr<functionspace::Edges> edges_fs( new functionspace::Edges(mesh,Halo(mesh)) );
  edges_fs->haloExchange(edges.field( "dual_normals" ));

  make_dual_normals_outward_new(mesh);
}

#if !DEPRECATE_OLD_FUNCTIONSPACE
void build_median_dual_mesh_convert_to_old( Mesh& mesh )
{
  deprecated::FunctionSpace&  edges = mesh.function_space("edges");
  ArrayView<double,2> dual_normals  ( edges.create_field<double>("dual_normals",2) );
  const ArrayView<double,2> new_dual_normals( mesh.edges().field( "dual_normals") );
  for( size_t jedge=0; jedge<mesh.edges().size(); ++jedge )
  {
    dual_normals(jedge,LON) = new_dual_normals(jedge,LON);
    dual_normals(jedge,LAT) = new_dual_normals(jedge,LAT);
  }
}
#endif


void build_median_dual_mesh( Mesh& mesh )
{
  build_median_dual_mesh_new( mesh );
#if !DEPRECATE_OLD_FUNCTIONSPACE
  build_median_dual_mesh_convert_to_old( mesh );
#endif
}


#if !DEPRECATE_OLD_FUNCTIONSPACE

void build_centroids( deprecated::FunctionSpace& func_space, ArrayView<double,2>& lonlat);

void build_centroid_dual_mesh( Mesh& mesh )
{
  mesh::Nodes& nodes   = mesh.nodes();
  ArrayView<double,2> lonlat        ( nodes.lonlat() );
  ArrayView<double,1> dual_volumes  ( nodes.add( Field::create<double>( "dual_volumes", make_shape(nodes.size(),1) ) ) );

  deprecated::FunctionSpace& quads       = mesh.function_space( "quads" );
  deprecated::FunctionSpace& triags      = mesh.function_space( "triags" );
  deprecated::FunctionSpace& edges       = mesh.function_space( "edges" );

  build_centroids(quads,  lonlat);
  build_centroids(triags, lonlat);
  build_centroids(edges,  lonlat);

  add_centroid_dual_volume_contribution( mesh, dual_volumes );

  build_dual_normals( mesh );

  build_skewness( mesh );

  eckit::SharedPtr<functionspace::Nodes> nodes_fs( new functionspace::Nodes(mesh,Halo(mesh)) );
  nodes_fs->haloExchange(nodes.field("dual_volumes"));

  ArrayView<double,2> dual_normals  ( edges.field( "dual_normals" ) );
  edges.parallelise();
  edges.halo_exchange().execute(dual_normals);
}



void build_centroids( deprecated::FunctionSpace& func_space, ArrayView<double,2>& lonlat)
{
  if( !func_space.has_field("centroids") )
  {
    int nb_elems = func_space.shape(0);
    IndexView<int,2> elem_nodes( func_space.field( "nodes" ) );
    int nb_nodes_per_elem = elem_nodes.shape(1);
    ArrayView<double,2> elem_centroids( func_space.create_field<double>( "centroids", 2 ) );
    for (size_t e=0; e<nb_elems; ++e)
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
#else
void build_centroid_dual_mesh( Mesh& mesh )
{
  NOTIMP;
}
#endif

Array* build_centroids_lonlat( const mesh::HybridElements& elements, const Field& field_lonlat )
{
  const ArrayView<double,2> lonlat( field_lonlat );
  Array* array_centroids = Array::create<double>(make_shape(elements.size(),2));
  ArrayView<double,2> centroids( *array_centroids );
  size_t nb_elems = elements.size();
  const mesh::HybridElements::Connectivity& elem_nodes = elements.node_connectivity();
  for (size_t e=0; e<nb_elems; ++e)
  {
    centroids(e,LON) = 0.;
    centroids(e,LAT) = 0.;
    const size_t nb_nodes_per_elem = elem_nodes.cols(e);
    const double average_coefficient = 1./static_cast<double>(nb_nodes_per_elem);
    for (size_t n=0; n<nb_nodes_per_elem; ++n)
    {
      centroids(e,LON) += lonlat( elem_nodes(e,n), LON );
      centroids(e,LAT) += lonlat( elem_nodes(e,n), LAT );
    }
    centroids(e,LON) *= average_coefficient;
    centroids(e,LAT) *= average_coefficient;
  }
  return array_centroids;
}

void add_median_dual_volume_contribution_cells(
    const mesh::HybridElements& cells,
    const mesh::HybridElements& edges,
    const mesh::Nodes& nodes,
    Array& array_dual_volumes )
{
  ArrayView<double,1> dual_volumes( array_dual_volumes );

  const ArrayView<double,2> lonlat ( nodes.lonlat() );
  const ArrayView<double,2> cell_centroids ( cells.field("centroids_lonlat") );
  const ArrayView<double,2> edge_centroids ( edges.field("centroids_lonlat") );
  const mesh::HybridElements::Connectivity& cell_edge_connectivity = cells.edge_connectivity();
  const mesh::HybridElements::Connectivity& edge_node_connectivity = edges.node_connectivity();

  // special ordering for bit-identical results
  size_t nb_cells = cells.size();
  std::vector<Node> ordering(nb_cells);
  for (int jcell=0; jcell<nb_cells; ++jcell)
    ordering[jcell] = Node( util::unique_lonlat(cell_centroids[jcell]), jcell );
  std::sort( ordering.data(), ordering.data()+nb_cells );

  for (size_t jcell=0; jcell<nb_cells; ++jcell)
  {
    size_t icell = ordering[jcell].i;
    double x0 = cell_centroids(icell,LON);
    double y0 = cell_centroids(icell,LAT);

    for (size_t jedge=0; jedge<cell_edge_connectivity.cols(icell); ++jedge)
    {
      size_t iedge = cell_edge_connectivity(icell,jedge);
      double x1 = edge_centroids(iedge,LON);
      double y1 = edge_centroids(iedge,LAT);
      for( size_t jnode=0; jnode<2; ++jnode )
      {
        size_t inode = edge_node_connectivity(iedge,jnode);
        double x2 = lonlat(inode,LON);
        double y2 = lonlat(inode,LAT);
        double triag_area = std::abs( x0*(y1-y2)+x1*(y2-y0)+x2*(y0-y1) )*0.5;
        dual_volumes(inode) += triag_area;
      }
    }
  }
}

void add_median_dual_volume_contribution_poles(
    const mesh::HybridElements& edges,
    const mesh::Nodes& nodes,
    Array& array_dual_volumes )
{
  ArrayView<double,1> dual_volumes( array_dual_volumes );
  const ArrayView<double,2> lonlat ( nodes.lonlat() );
  const ArrayView<double,2> edge_centroids ( edges.field("centroids_lonlat") );
  const mesh::HybridElements::Connectivity& edge_node_connectivity = edges.node_connectivity();
  const mesh::HybridElements::Connectivity& edge_cell_connectivity = edges.cell_connectivity();

  const size_t nb_edges = edges.size();
  std::map<idx_t,std::vector<idx_t> > node_to_bdry_edge;
  for(size_t jedge=0; jedge<nb_edges; ++jedge)
  {
    if (    edge_cell_connectivity(jedge,0)!=edge_cell_connectivity.missing_value()
         && edge_cell_connectivity(jedge,1)==edge_cell_connectivity.missing_value() )
    {
      node_to_bdry_edge[ edge_node_connectivity(jedge,0) ].push_back(jedge);
      node_to_bdry_edge[ edge_node_connectivity(jedge,1) ].push_back(jedge);
    }
  }

  const double tol = 1.e-6;
  double min[2], max[2];
  global_bounding_box( nodes, min, max );

  std::map<idx_t,std::vector<idx_t> >::iterator it;
  for( it=node_to_bdry_edge.begin(); it!=node_to_bdry_edge.end(); ++it)
  {
    const size_t jnode = (*it).first;
    std::vector<idx_t>& bdry_edges = (*it).second;
    const double x0 = lonlat(jnode,LON);
    const double y0 = lonlat(jnode,LAT);
    double x1, y1, y2;
    for (size_t jedge = 0; jedge < bdry_edges.size(); ++jedge)
    {
      const size_t iedge = bdry_edges[jedge];
      x1 = edge_centroids(iedge,LON);
      y1 = edge_centroids(iedge,LAT);

      y2 = 0.;
      if ( std::abs(y1-max[LAT])<tol )
        y2 = 90.;
      else if ( std::abs(y1-min[LAT])<tol )
        y2 = -90.;

      if( y2!=0 )
      {
        const double quad_area = std::abs( (x1-x0)*(y2-y0) );
        dual_volumes(jnode) += quad_area;
      }
    }
  }

}

#if !DEPRECATE_OLD_FUNCTIONSPACE
void add_centroid_dual_volume_contribution(
    Mesh& mesh,
    ArrayView<double,1>& dual_volumes )
{
  mesh::Nodes& nodes = mesh.nodes();
  deprecated::FunctionSpace& edges = mesh.function_space("edges");
  ArrayView<gidx_t,1> node_glb_idx  ( nodes.field("glb_idx"    ) );
  ArrayView<double,2> edge_centroids( edges.field("centroids"  ) );
  IndexView<int,   2> edge_nodes    ( edges.field("nodes"      ) );
  ArrayView<gidx_t,1> edge_glb_idx  ( edges.field("glb_idx"    ) );
  IndexView<int,   2> edge_to_elem  ( edges.field("to_elem"    ) );
  ArrayView<double,2> node_lonlat   ( nodes.lonlat() );
  std::vector< ArrayView<double,2> > elem_centroids(mesh.nb_function_spaces());
  for(size_t f = 0; f < mesh.nb_function_spaces(); ++f)
  {
    deprecated::FunctionSpace& elements = mesh.function_space(f);
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
#endif

void build_dual_normals_new( Mesh& mesh )
{
  ArrayView<double,2> elem_centroids( mesh.cells().field("centroids_lonlat") );

  mesh::Nodes&  nodes = mesh.nodes();
  mesh::HybridElements&  edges = mesh.edges();
  const size_t nb_edges = edges.size();

  ArrayView<double,2> node_lonlat( nodes.lonlat() );
  double min[2], max[2];
  global_bounding_box( nodes, min, max );
  double tol = 1.e-6;

  double xl, yl, xr, yr;
  ArrayView<double,2> edge_centroids( edges.field("centroids_lonlat") );
  ArrayView<double,2> dual_normals  ( edges.add( Field::create<double>("dual_normals",make_shape(nb_edges,2)) ) );

  const mesh::HybridElements::Connectivity& edge_node_connectivity = edges.node_connectivity();
  const mesh::HybridElements::Connectivity& edge_cell_connectivity = edges.cell_connectivity();

  std::map<idx_t,std::vector<idx_t> > node_to_bdry_edge;
  for(size_t jedge=0; jedge<nb_edges; ++jedge)
  {
    if (    edge_cell_connectivity(jedge,0)!=edge_cell_connectivity.missing_value()
         && edge_cell_connectivity(jedge,1)==edge_cell_connectivity.missing_value() )
    {
      node_to_bdry_edge[ edge_node_connectivity(jedge,0) ].push_back(jedge);
      node_to_bdry_edge[ edge_node_connectivity(jedge,1) ].push_back(jedge);
    }
  }

  for (int edge=0; edge<nb_edges; ++edge)
  {
    if( edge_cell_connectivity(edge,0) == edge_cell_connectivity.missing_value() )
    {
      // this is a pole edge
      // only compute for one node
      for (size_t n=0; n<2; ++n)
      {
        idx_t node = edge_node_connectivity(edge,n);
        std::vector<idx_t>& bdry_edges = node_to_bdry_edge[node];
        double x[2];
        int cnt=0;
        for (size_t jedge = 0; jedge < bdry_edges.size(); ++jedge)
        {
          idx_t bdry_edge = bdry_edges[jedge];
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
      idx_t left_elem          = edge_cell_connectivity(edge,0);
      idx_t right_elem         = edge_cell_connectivity(edge,1);
      xl = elem_centroids(left_elem,LON);
      yl = elem_centroids(left_elem,LAT);
      if( right_elem == edge_cell_connectivity.missing_value() )
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
        xr = elem_centroids(right_elem,LON);
        yr = elem_centroids(right_elem,LAT);
      }

      dual_normals(edge,LON) =  yl-yr;
      dual_normals(edge,LAT) = -xl+xr;
    }
  }
}

#if !DEPRECATE_OLD_FUNCTIONSPACE

void build_dual_normals( Mesh& mesh )
{

  std::vector< ArrayView<double,2> > elem_centroids( mesh.nb_function_spaces() );
  for (size_t func_space_idx = 0; func_space_idx < mesh.nb_function_spaces(); ++func_space_idx)
  {
    deprecated::FunctionSpace& func_space = mesh.function_space(func_space_idx);
    if( func_space.has_field("centroids") )
      elem_centroids[func_space_idx] = ArrayView<double,2>( func_space.field("centroids") );
  }

  mesh::Nodes&  nodes = mesh.nodes();
  ArrayView<double,2> node_lonlat( nodes.lonlat() );
  double min[2], max[2];
  global_bounding_box( nodes, min, max );
  double tol = 1.e-6;

  double xl, yl, xr, yr;
  deprecated::FunctionSpace&  edges = mesh.function_space("edges");
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
#endif

void make_dual_normals_outward_new( Mesh& mesh )
{
  mesh::Nodes&  nodes = mesh.nodes();
  ArrayView<double,2> node_lonlat( nodes.lonlat() );

  mesh::HybridElements& edges = mesh.edges();
  const mesh::HybridElements::Connectivity& edge_cell_connectivity = edges.cell_connectivity();
  const mesh::HybridElements::Connectivity& edge_node_connectivity = edges.node_connectivity();
  ArrayView<double,2> dual_normals  ( edges.field("dual_normals") );
  const size_t nb_edges = edges.size();

  for (size_t edge=0; edge<nb_edges; ++edge)
  {
    if( edge_cell_connectivity(edge,0) != edge_cell_connectivity.missing_value() )
    {
      // Make normal point from node 1 to node 2
      const size_t ip1 = edge_node_connectivity(edge,0);
      const size_t ip2 = edge_node_connectivity(edge,1);
      double dx = node_lonlat( ip2, LON ) - node_lonlat( ip1, LON );
      double dy = node_lonlat( ip2, LAT ) - node_lonlat( ip1, LAT );
      if( dx*dual_normals(edge,LON) + dy*dual_normals(edge,LAT) < 0 )
      {
        dual_normals(edge,LON) = - dual_normals(edge,LON);
        dual_normals(edge,LAT) = - dual_normals(edge,LAT);
      }
    }
  }
}

#if !DEPRECATE_OLD_FUNCTIONSPACE
void make_dual_normals_outward( Mesh& mesh )
{

  mesh::Nodes&  nodes = mesh.nodes();
  ArrayView<double,2> node_lonlat( nodes.lonlat() );

  deprecated::FunctionSpace&  edges = mesh.function_space("edges");
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
#endif

#if !DEPRECATE_OLD_FUNCTIONSPACE
void build_skewness( Mesh& mesh )
{
  std::vector< ArrayView<double,2> > elem_centroids( mesh.nb_function_spaces() );
  for (size_t func_space_idx = 0; func_space_idx < mesh.nb_function_spaces(); ++func_space_idx)
  {
    deprecated::FunctionSpace& func_space = mesh.function_space(func_space_idx);
    if( func_space.has_field("centroids") )
      elem_centroids[func_space_idx] = ArrayView<double,2>( func_space.field("centroids") );
  }

  mesh::Nodes&  nodes = mesh.nodes();
  ArrayView<double,2> node_lonlat( nodes.lonlat() );
  double min[2], max[2];
  global_bounding_box( nodes, min, max );
  double tol = 1.e-6;

  double x1, y1, x2, y2, xc1, yc1, xc2, yc2, xi, yi;
  deprecated::FunctionSpace&  edges = mesh.function_space("edges");
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
#endif



void build_brick_dual_mesh( Mesh& mesh )
{
  if( const grids::ReducedGrid* g = dynamic_cast<const grids::ReducedGrid*>(&mesh.grid()) )
  {
    if( eckit::mpi::size() != 1 )
      throw eckit::UserError("Cannot build_brick_dual_mesh with more than 1 task",Here());

    mesh::Nodes& nodes   = mesh.nodes();
    ArrayView<double,2> lonlat        ( nodes.lonlat() );
    ArrayView<double,1> dual_volumes  ( nodes.add( Field::create<double>("dual_volumes",make_shape(nodes.size(),1) ) ) );
    ArrayView<gidx_t,1> gidx  ( nodes.global_index() );

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

    eckit::SharedPtr<functionspace::Nodes> nodes_fs( new functionspace::Nodes(mesh,Halo(mesh)) );
    nodes_fs->haloExchange(nodes.field("dual_volumes"));
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

#if !DEPRECATE_OLD_FUNCTIONSPACE
void atlas__build_centroid_dual_mesh ( Mesh* mesh) {
  ATLAS_ERROR_HANDLING( build_centroid_dual_mesh(*mesh) );
}
#endif
// ------------------------------------------------------------------

} // namespace actions
} // namespace atlas

