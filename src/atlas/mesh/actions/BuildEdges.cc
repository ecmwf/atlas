/*
 * (C) Copyright 1996-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include <set>
#include <cmath>
#include <limits>
#include <iostream>
#include <stdexcept>
#include "atlas/atlas_config.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/mesh/actions/BuildEdges.h"
#include "atlas/mesh/Nodes.h"
#include "atlas/mesh/ElementType.h"
#include "atlas/mesh/HybridElements.h"
#include "atlas/mesh/Elements.h"
#include "atlas/field/Field.h"
#include "atlas/functionspace/FunctionSpace.h"
#include "atlas/internals/AccumulateFaces.h"
#include "atlas/internals/Bitflags.h"
#include "atlas/internals/Unique.h"
#include "atlas/internals/LonLatMicroDeg.h"
#include "atlas/internals/Functions.h"
#include "atlas/internals/Parameters.h"
#include "atlas/array/ArrayView.h"
#include "atlas/array/Array.h"
#include "atlas/array/IndexView.h"
#include "atlas/util/runtime/ErrorHandling.h"

using atlas::internals::accumulate_facets;
using atlas::internals::Topology;
using atlas::internals::UniqueLonLat;
using atlas::internals::microdeg;

namespace atlas {
namespace mesh {
namespace actions {

#if !DEPRECATE_OLD_FUNCTIONSPACE
void build_edges_convert_to_old( Mesh& mesh );
#endif

//----------------------------------------------------------------------------------------------------------------------

namespace {
struct Sort
{
  Sort() {}
  Sort(gidx_t gid, int idx)
  {
    g = gid;
    i = idx;
  }
  gidx_t g;
  int i;
  bool operator < (const Sort& other) const
  {
    return ( g < other.g );
  }
};
}

#if !DEPRECATE_OLD_FUNCTIONSPACE
void build_element_to_edge_connectivity_convert_to_old( Mesh& mesh )
{
  std::vector< array::IndexView<int,2> > elem_to_edge( mesh.nb_function_spaces() );
  size_t etype=0;
  for( size_t func_space_idx=0; func_space_idx<mesh.nb_function_spaces(); ++func_space_idx)
  {
    deprecated::FunctionSpace& func_space = mesh.function_space(func_space_idx);
    if( func_space.metadata().get<long>("type") == Entity::ELEMS )
    {
      int nb_edges_per_elem;
      if (func_space.name() == "quads")  nb_edges_per_elem = 4;
      if (func_space.name() == "triags") nb_edges_per_elem = 3;
      elem_to_edge[etype++] =
          array::IndexView<int,2>(func_space.create_field<int>("to_edge",nb_edges_per_elem));
    }
  }
  const MultiBlockConnectivity& cell_edge_connectivity = mesh.cells().edge_connectivity();
  for( size_t t=0; t<mesh.cells().nb_types(); ++t)
  {
    deprecated::FunctionSpace& func_space = mesh.function_space(t);
    size_t nb_edges_per_elem = mesh.cells().element_type(t).nb_edges();

    ASSERT(elem_to_edge[t].shape(0) == mesh.cells().elements(t).size());
    ASSERT(elem_to_edge[t].shape(1) == nb_edges_per_elem);

    for( size_t jelem=0; jelem<mesh.cells().elements(t).size(); ++jelem )
    {
      for( size_t jedge=0; jedge<nb_edges_per_elem; ++jedge )
      {
        elem_to_edge[t](jelem,jedge) = mesh.cells().edge_connectivity()(t,jelem,jedge);
      }
    }
  }

}
#endif

void build_element_to_edge_connectivity_new( Mesh& mesh )
{
  mesh::HybridElements::Connectivity& cell_edge_connectivity = mesh.cells().edge_connectivity();
  cell_edge_connectivity.clear();

  // Allocate cell_edge_connectivity
  for( size_t t=0; t<mesh.cells().nb_types(); ++t)
  {
    size_t nb_elements = mesh.cells().elements(t).size();
    size_t nb_edges_per_elem = mesh.cells().element_type(t).nb_edges();
    std::vector<idx_t> init(mesh.cells().elements(t).size()*nb_edges_per_elem,cell_edge_connectivity.missing_value());
    cell_edge_connectivity.add(nb_elements,nb_edges_per_elem,init.data());
  }

  size_t nb_edges = mesh.edges().size();
  mesh::HybridElements::Connectivity& edge_cell_connectivity = mesh.edges().cell_connectivity();
  mesh::HybridElements::Connectivity& edge_node_connectivity = mesh.edges().node_connectivity();

  bool has_pole_edges(false);
  array::ArrayView<int,1> is_pole_edge;
  if( mesh.edges().has_field("is_pole_edge") )
  {
    has_pole_edges = true;
    is_pole_edge = array::ArrayView<int,1>( mesh.edges().field("is_pole_edge") );
  }

  // Sort edges for bit-reproducibility
  std::vector<Sort> edge_sort(nb_edges);
  {
    UniqueLonLat compute_uid( mesh.nodes() );

    for( size_t jedge=0; jedge<nb_edges; ++jedge )
      edge_sort[jedge] = Sort( compute_uid(edge_node_connectivity.row(jedge)), jedge );

    std::sort( edge_sort.data(), edge_sort.data()+nb_edges );
  }

  // Fill in cell_edge_connectivity
  std::vector<size_t> edge_cnt( mesh.cells().size() );
  for( size_t jedge=0; jedge<nb_edges; ++jedge)
  {
    int iedge = edge_sort[jedge].i;
    for( size_t j=0; j<2; ++j)
    {
      idx_t elem = edge_cell_connectivity(iedge,j);

      if ( elem != edge_cell_connectivity.missing_value() )
      {
        cell_edge_connectivity.set(elem,edge_cnt[elem]++, iedge);
      }
      else
      {
        if( !( has_pole_edges && is_pole_edge(iedge) ) )
        {
          if( j==0 )
            throw eckit::SeriousBug("edge has no element connected",Here());
        }
      }
    }
  }


  // Verify that all edges have been found
  for( size_t jcell=0; jcell<mesh.cells().size(); ++jcell )
  {
    for( size_t jcol=0; jcol<cell_edge_connectivity.cols(jcell); ++jcol )
    {
      if( cell_edge_connectivity(jcell,jcol) == cell_edge_connectivity.missing_value() )
      {
        const array::ArrayView<gidx_t,1> gidx (mesh.nodes().global_index() );
        std::stringstream msg; msg << "Could not find edge " << jcol << " for " << mesh.cells().name(jcell) << " elem " << jcell << " with nodes ( ";
        for( size_t jnode=0; jnode<mesh.cells().node_connectivity().cols(jcell); ++jnode )
        {
          msg << gidx( mesh.cells().node_connectivity()(jcell,jnode) ) <<" ";
        }
        msg << ")";
        throw eckit::SeriousBug(msg.str(),Here());
      }
    }
  }
}


void build_element_to_edge_connectivity( Mesh& mesh )
{
  build_element_to_edge_connectivity_new(mesh);
#if !DEPRECATE_OLD_FUNCTIONSPACE
  build_element_to_edge_connectivity_convert_to_old(mesh);
#endif
}


void build_node_to_edge_connectivity_new( Mesh& mesh )
{
  mesh::Nodes& nodes = mesh.nodes();
  const size_t nb_edges = mesh.edges().size();

  mesh::HybridElements::Connectivity& edge_node_connectivity = mesh.edges().node_connectivity();

  std::vector<size_t> to_edge_size (nodes.size(),0);
  for(size_t jedge = 0; jedge < nb_edges; ++jedge)
  {
    for( int j=0; j<2; ++j)
    {
      ++to_edge_size[ edge_node_connectivity(jedge,j) ];
    }
  }

  mesh::Nodes::Connectivity& node_to_edge = nodes.edge_connectivity();
  node_to_edge.add( nodes.size(), to_edge_size.data() );
  for( size_t jnode=0; jnode<nodes.size(); ++jnode )
    to_edge_size[jnode] = 0;

  UniqueLonLat compute_uid( nodes );
  std::vector<Sort> edge_sort(nb_edges);
  for( size_t jedge=0; jedge<nb_edges; ++jedge )
    edge_sort[jedge] = Sort( compute_uid(edge_node_connectivity.row(jedge)), jedge );
  std::stable_sort( edge_sort.data(), edge_sort.data()+nb_edges );

  for( size_t jedge=0; jedge<nb_edges; ++jedge)
  {
    size_t iedge = edge_sort[jedge].i;
    for( size_t j=0; j<2; ++j)
    {
      idx_t node = edge_node_connectivity(iedge,j);
      node_to_edge.set(node, to_edge_size[node]++, iedge );
    }
  }
}

#if !DEPRECATE_OLD_FUNCTIONSPACE
void build_node_to_edge_connectivity_convert_to_old( Mesh& mesh )
{
  mesh::Nodes& nodes = mesh.nodes();
  const mesh::Nodes::Connectivity& node_edge_connectivity = mesh.nodes().edge_connectivity();

  // Get max_edge_cnt
  array::ArrayView<int,1> to_edge_size ( nodes.add( field::Field::create<int>( "to_edge_size", array::make_shape(nodes.size(),1) ) ) );
  int max_edge_cnt(0);
  for( size_t jnode=0; jnode<nodes.size(); ++jnode )
  {
    to_edge_size(jnode) = node_edge_connectivity.cols(jnode);
    max_edge_cnt = std::max(max_edge_cnt,to_edge_size(jnode));
  }

  ECKIT_MPI_CHECK_RESULT( MPI_Allreduce( MPI_IN_PLACE, &max_edge_cnt, 1, MPI_INT, MPI_MAX, eckit::mpi::comm() ) );

  array::IndexView<int,2> node_to_edge ( nodes.add( field::Field::create<int>("to_edge",array::make_shape(nodes.size(),max_edge_cnt))) );

  for( size_t jnode=0; jnode<nodes.size(); ++jnode )
  {
    for( size_t jedge=0; jedge<node_edge_connectivity.cols(jnode); ++jedge )
    {
      node_to_edge(jnode,jedge) = node_edge_connectivity(jnode,jedge);
    }
  }
}
#endif

void build_node_to_edge_connectivity( Mesh& mesh )
{
  build_node_to_edge_connectivity_new( mesh );
#if !DEPRECATE_OLD_FUNCTIONSPACE
  build_node_to_edge_connectivity_convert_to_old( mesh );
#endif
}


void accumulate_pole_edges( mesh::Nodes& nodes, std::vector<idx_t>& pole_edge_nodes, size_t& nb_pole_edges )
{
  enum { NORTH=0, SOUTH=1 };

  array::ArrayView<double,2> lonlat    ( nodes.lonlat() );
  array::ArrayView<int,   1> flags     ( nodes.field( "flags"       ) );
  array::ArrayView<int,   1> part      ( nodes.partition() );
  const size_t nb_nodes = nodes.size();

  double min[2], max[2];
  min[internals::LON] =  std::numeric_limits<double>::max();
  min[internals::LAT] =  std::numeric_limits<double>::max();
  max[internals::LON] = -std::numeric_limits<double>::max();
  max[internals::LAT] = -std::numeric_limits<double>::max();
  for (size_t node=0; node<nb_nodes; ++node)
  {
    min[internals::LON] = std::min( min[internals::LON], lonlat(node,internals::LON) );
    min[internals::LAT] = std::min( min[internals::LAT], lonlat(node,internals::LAT) );
    max[internals::LON] = std::max( max[internals::LON], lonlat(node,internals::LON) );
    max[internals::LAT] = std::max( max[internals::LAT], lonlat(node,internals::LAT) );
  }

  eckit::mpi::all_reduce( min, 2, eckit::mpi::min() );
  eckit::mpi::all_reduce( max, 2, eckit::mpi::max() );

  double tol = 1e-6;

  // Collect all nodes closest to poles
  std::vector< std::set<int> > pole_nodes(2);
  for (size_t node=0; node<nb_nodes; ++node)
  {
      if ( std::abs(lonlat(node,internals::LAT)-max[internals::LAT])<tol )
      {
        pole_nodes[NORTH].insert(node);
      }
      else if ( std::abs(lonlat(node,internals::LAT)-min[internals::LAT])<tol )
      {
        pole_nodes[SOUTH].insert(node);
      }
  }

  // Sanity check
  {
    for( size_t NS = 0; NS<2; ++NS )
    {
      int npart=-1;
      for( std::set<int>::iterator it=pole_nodes[NS].begin(); it!=pole_nodes[NS].end(); ++it)
      {
        int node = *it;
        if( npart == -1 ) npart = part(node);
        else if ( part(node) != npart )
        {
          // Not implemented yet, when pole-lattitude is split.
          std::stringstream msg;
          msg << "Split pole-latitude is not supported yet...  node "<<node<<"[p"<<part(node)<<"] should belong to part " << npart;
          throw eckit::NotImplemented(msg.str(),Here());
        }
      }
    }
  }

  // Create connections over the poles and store in pole_edge_nodes
  nb_pole_edges = 0;
  for( size_t NS = 0; NS<2; ++NS )
  {
    for( std::set<int>::iterator it=pole_nodes[NS].begin(); it!=pole_nodes[NS].end(); ++it)
    {
      int node = *it;
      if( !Topology::check(flags(node),Topology::PERIODIC|Topology::GHOST) )
      {
        int x2 = microdeg( lonlat(node,internals::LON) + 180. );
        for( std::set<int>::iterator itr=pole_nodes[NS].begin(); itr!=pole_nodes[NS].end(); ++itr)
        {
          int other_node = *itr;
          if( microdeg( lonlat(other_node,internals::LON) ) == x2 )
          {
            if( !Topology::check(flags(other_node),Topology::PERIODIC) )
            {
              pole_edge_nodes.push_back(node);
              pole_edge_nodes.push_back(other_node);
              ++nb_pole_edges;
            }
          }
        }
      }
    }
  }

}


struct ComputeUniquePoleEdgeIndex
{
  ComputeUniquePoleEdgeIndex( const mesh::Nodes& nodes )
  {
    lonlat = array::ArrayView<double,2> ( nodes.lonlat() );
  }

  gidx_t operator()( const array::IndexView<int,1>& edge_nodes ) const
  {
    double centroid[2];
    centroid[internals::LON] = 0.;
    centroid[internals::LAT] = 0.;
    for( size_t jnode=0; jnode<2; ++jnode )
    {
      centroid[internals::LON] += lonlat( edge_nodes(jnode), internals::LON );
      centroid[internals::LAT] += lonlat( edge_nodes(jnode), internals::LAT );
    }
    centroid[internals::LON] /= 2.;
    centroid[internals::LAT] /= 2.;
    if( centroid[internals::LAT] > 0 )
      centroid[internals::LAT] =  90.;
    else
      centroid[internals::LAT] = -90.;
    /// FIXME make this into `internals::unique_lonlat(centroid)` but this causes weird parallel behavior
    return internals::detail::unique32( internals::microdeg(centroid[internals::LON]), internals::microdeg(centroid[internals::LON]) );
  }

  array::ArrayView<double,2> lonlat;
};


void build_edges_new( Mesh& mesh )
{
  mesh::Nodes& nodes   = mesh.nodes();
  array::ArrayView<gidx_t,1> glb_idx ( nodes.global_index() );
  array::ArrayView<int,1>    part    ( nodes.partition() );
  array::ArrayView<double,2> lonlat  ( nodes.lonlat() );

  size_t nb_nodes = nodes.size();

  // storage for edge-to-node-connectivity shape=(nb_edges,2)
  std::vector< idx_t > edge_nodes_data;
  std::vector< idx_t > edge_to_elem_data;
  size_t nb_edges;
  size_t nb_inner_edges;
  idx_t missing_value;

  accumulate_facets(mesh.cells(),mesh.nodes(),edge_nodes_data,edge_to_elem_data,nb_edges,nb_inner_edges,missing_value);
  // Build edges
  mesh.edges().add( new mesh::temporary::Line(), nb_edges, edge_nodes_data.data() );
  mesh::HybridElements::Connectivity& edge_nodes = mesh.edges().node_connectivity();
  mesh::HybridElements::Connectivity& cell_nodes = mesh.cells().node_connectivity();

  UniqueLonLat compute_uid( nodes );

  array::IndexView<idx_t,1>  edge_ridx     ( mesh.edges().remote_index() );
  array::ArrayView<int,1>    edge_part     ( mesh.edges().partition() );
  array::ArrayView<gidx_t,1> edge_glb_idx  ( mesh.edges().global_index() );

  for( size_t edge=0; edge<nb_edges; ++edge )
  {
    const int ip1 = edge_nodes(edge,0);
    const int ip2 = edge_nodes(edge,1);
    if( compute_uid(ip1) > compute_uid(ip2) )
    {
      idx_t swapped[2] = {ip2,ip1};
      edge_nodes.set(edge,swapped);
    }

    ASSERT( size_t(edge_nodes(edge,0)) < nb_nodes );
    ASSERT( size_t(edge_nodes(edge,1)) < nb_nodes );
    edge_glb_idx(edge)   = compute_uid(edge_nodes.row(edge));
    edge_part(edge)      = std::min( part(edge_nodes(edge,0)), part(edge_nodes(edge,1) ) );
    edge_ridx(edge)      = edge;

    const idx_t e1 = edge_to_elem_data[2*edge+0];
    const idx_t e2 = edge_to_elem_data[2*edge+1];

    if( e2 != cell_nodes.missing_value() )
    {
      // Swap order to ensure bit-reproducibility
      if( compute_uid(cell_nodes.row(e1)) > compute_uid(cell_nodes.row(e2)) )
      {
        edge_to_elem_data[edge*2+0] = e2;
        edge_to_elem_data[edge*2+1] = e1;
      }
    }
  }

  mesh.edges().cell_connectivity().add( nb_edges, 2, edge_to_elem_data.data() );
}

#if !DEPRECATE_OLD_FUNCTIONSPACE
void build_edges_convert_to_old( Mesh& mesh )
{
  deprecated::FunctionSpace& quads       = mesh.function_space( "quads" );
  deprecated::FunctionSpace& triags      = mesh.function_space( "triags" );
  const size_t nb_edges = mesh.edges().size();
  if( ! mesh.has_function_space("edges") )
  {
    mesh.create_function_space( "edges", "shapefunc", array::make_shape(nb_edges,deprecated::FunctionSpace::UNDEF_VARS) );
  }
  deprecated::FunctionSpace& edges = mesh.function_space("edges");
  edges.metadata().set<long>("type",Entity::FACES);
  edges.resize(array::make_shape(nb_edges,deprecated::FunctionSpace::UNDEF_VARS));

  if( ! edges.has_field("nodes")      )  edges.create_field<int>("nodes",     2);
  if( ! edges.has_field("glb_idx")    )  edges.create_field<gidx_t>("glb_idx",   1);
  if( ! edges.has_field("partition")  )  edges.create_field<int>("partition", 1);
  if( ! edges.has_field("to_elem")    )  edges.create_field<int>("to_elem",   4);
  if( ! edges.has_field("remote_idx") )  edges.create_field<int>("remote_idx",1);

  array::IndexView<int,2> edge_nodes   ( edges.field( "nodes"      ) );
  array::ArrayView<gidx_t,1> edge_glb_idx ( edges.field( "glb_idx"    ) );
  array::ArrayView<int,1> edge_part    ( edges.field( "partition"  ) );
  array::IndexView<int,1> edge_ridx    ( edges.field( "remote_idx" ) );
  array::IndexView<int,3> edge_to_elem ( edges.field( "to_elem"    ).data<int>(), array::make_shape(nb_edges,2,2) );

  const array::ArrayView<gidx_t,1> edge_glb_idx_new ( mesh.edges().field( "glb_idx"    ) );
  const array::ArrayView<int,1> edge_part_new    ( mesh.edges().field( "partition"  ) );
  const array::IndexView<int,1> edge_ridx_new    ( mesh.edges().field( "remote_idx" ) );
  const mesh::HybridElements::Connectivity&  edge_cell_connectivity = mesh.edges().cell_connectivity();
  const mesh::HybridElements::Connectivity&  edge_node_connectivity = mesh.edges().node_connectivity();
  for( size_t jedge=0; jedge<nb_edges; ++jedge )
  {
    edge_glb_idx(jedge) = edge_glb_idx_new(jedge);
    edge_part(jedge) = edge_part_new(jedge);
    edge_ridx(jedge) = edge_ridx_new(jedge);
    for( size_t jnode=0; jnode<2; ++jnode )
    {
      edge_nodes(jedge,jnode) = edge_node_connectivity(jedge,jnode);
    }
    for( size_t jelem=0; jelem<2; ++jelem )
    {
      if( edge_cell_connectivity(jedge,jelem) != edge_cell_connectivity.missing_value() )
      {
        int f = mesh.cells().type_idx(edge_cell_connectivity(jedge,jelem));
        int e = edge_cell_connectivity(jedge,jelem) - mesh.cells().elements(f).begin();
        edge_to_elem(jedge,jelem,0) = f;
        edge_to_elem(jedge,jelem,1) = e;
      }
      else
      {
        edge_to_elem(jedge,jelem,0) = -1;
        edge_to_elem(jedge,jelem,1) = -1;
      }
    }
  }
}
#endif


void build_edges( Mesh& mesh )
{
  build_edges_new(mesh);

#if !DEPRECATE_OLD_FUNCTIONSPACE
  build_edges_convert_to_old(mesh);
#endif

  build_element_to_edge_connectivity(mesh);

}

void build_pole_edges_new( Mesh& mesh )
{
  mesh::Nodes& nodes   = mesh.nodes();
  mesh::HybridElements& edges = mesh.edges();

  size_t nb_cell_edges = edges.size();

  size_t nb_pole_edges;
  std::vector<idx_t> pole_edge_nodes;
  accumulate_pole_edges( nodes, pole_edge_nodes, nb_pole_edges );

  edges.add(new mesh::temporary::Line(), nb_pole_edges, pole_edge_nodes.data() );

  if( ! edges.has_field("is_pole_edge") )
    edges.add(field::Field::create<int>("is_pole_edge",array::make_shape(edges.size())));

  array::ArrayView<int,1> node_part       ( nodes.partition() );

  array::ArrayView<gidx_t,1> edge_glb_idx ( edges.global_index() );
  array::ArrayView<int,1> edge_part       ( edges.partition() );
  array::IndexView<int,1> edge_ridx       ( edges.remote_index() );
  array::ArrayView<int,1> is_pole_edge    ( edges.field( "is_pole_edge" ) );

  mesh::HybridElements::Connectivity& edge_nodes   = edges.node_connectivity();
  IrregularConnectivity& edge_to_elem = edges.cell_connectivity();
  edge_to_elem.add(nb_pole_edges,2);

  for(size_t edge=0; edge<nb_cell_edges; ++edge)
  {
    is_pole_edge(edge) = 0;
  }

  size_t cnt = 0;
  ComputeUniquePoleEdgeIndex compute_uid( nodes );
  for(size_t edge = nb_cell_edges; edge < nb_cell_edges + nb_pole_edges; ++edge)
  {
    idx_t ip1 = pole_edge_nodes[cnt++];
    idx_t ip2 = pole_edge_nodes[cnt++];
    idx_t enodes[] = {ip1,ip2};
    edge_nodes.set( edge, enodes );
    edge_glb_idx(edge)   = compute_uid( edge_nodes.row(edge) );
    edge_part(edge)      = std::min( node_part(edge_nodes(edge,0)), node_part(edge_nodes(edge,1) ) );
    edge_ridx(edge)      = edge;
    is_pole_edge(edge) = 1;
  }
}

#if !DEPRECATE_OLD_FUNCTIONSPACE
void build_pole_edges_convert_to_old( Mesh& mesh )
{
  mesh::Nodes& nodes   = mesh.nodes();
  deprecated::FunctionSpace& edges = mesh.function_space("edges");

  size_t nb_edges = edges.shape(0);
  edges.resize( array::make_shape(mesh.edges().size(), deprecated::FunctionSpace::UNDEF_VARS) );

  if( ! edges.has_field("is_pole_edge") )  edges.create_field<int>("is_pole_edge",1);

  array::IndexView<int,2>    edge_nodes   ( edges.field( "nodes"      ) );
  array::ArrayView<gidx_t,1> edge_glb_idx ( edges.field( "glb_idx"    ) );
  array::ArrayView<int,1>    edge_part    ( edges.field( "partition"  ) );
  array::IndexView<int,1>    edge_ridx    ( edges.field( "remote_idx" ) );
  array::ArrayView<int,1>    is_pole_edge ( edges.field( "is_pole_edge" ) );
  array::IndexView<int,3>    edge_to_elem ( edges.field( "to_elem"    ).data<int>(), array::make_shape(mesh.edges().size(),2,2) );

  mesh::HybridElements::Connectivity& new_edge_nodes   = mesh.edges().node_connectivity();
  mesh::HybridElements::Connectivity& new_edge_to_elem = mesh.edges().cell_connectivity();
  array::ArrayView<gidx_t,1> new_edge_glb_idx ( mesh.edges().global_index() );
  array::ArrayView<int,1>    new_edge_part    ( mesh.edges().partition() );
  array::IndexView<int,1>    new_edge_ridx    ( mesh.edges().remote_index() );
  array::ArrayView<int,1>    new_is_pole_edge ( mesh.edges().field( "is_pole_edge" ) );


  for(size_t edge=0; edge<nb_edges; ++edge)
  {
    is_pole_edge(edge) = 0;
  }

  for(size_t edge = nb_edges; edge < mesh.edges().size(); ++edge)
  {
    for( size_t i=0; i<2; ++i )
    {
      edge_nodes(edge,i)     = new_edge_nodes(edge,i);
      edge_to_elem(edge,i,0) = new_edge_to_elem(edge,i);
      edge_to_elem(edge,i,1) = new_edge_to_elem(edge,i);
    }
    edge_glb_idx(edge) = new_edge_glb_idx(edge);
    edge_part(edge)    = new_edge_part(edge);
    edge_ridx(edge)    = new_edge_ridx(edge);
    is_pole_edge(edge) = new_is_pole_edge(edge);
  }
}
#endif

void build_pole_edges( Mesh& mesh )
{
  build_pole_edges_new( mesh );

#if !DEPRECATE_OLD_FUNCTIONSPACE
  build_pole_edges_convert_to_old( mesh );
#endif
}


//----------------------------------------------------------------------------------------------------------------------
// C wrapper interfaces to C++ routines

void atlas__build_edges ( Mesh* mesh) {
  ATLAS_ERROR_HANDLING( build_edges(*mesh) );
}
void atlas__build_pole_edges ( Mesh* mesh) {
  ATLAS_ERROR_HANDLING( build_pole_edges(*mesh) );
}
void atlas__build_node_to_edge_connectivity ( Mesh* mesh) {
  ATLAS_ERROR_HANDLING( build_node_to_edge_connectivity(*mesh) );
}

//----------------------------------------------------------------------------------------------------------------------

} // namespace actions
} // namespace mesh
} // namespace atlas

