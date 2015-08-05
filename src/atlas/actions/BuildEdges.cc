/*
 * (C) Copyright 1996-2014 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include <iostream>
#include <stdexcept>
#include <cmath>
#include <limits>
#include <set>

#include "atlas/atlas_config.h"
#include "atlas/runtime/ErrorHandling.h"
#include "atlas/Mesh.h"
#include "atlas/FunctionSpace.h"
#include "atlas/Field.h"
#include "atlas/actions/BuildEdges.h"
#include "atlas/Parameters.h"
#include "atlas/Nodes.h"

#include "atlas/util/ArrayView.h"
#include "atlas/util/Array.h"
#include "atlas/util/IndexView.h"
#include "atlas/util/AccumulateFaces.h"
#include "atlas/util/Bitflags.h"
#include "atlas/util/Unique.h"
#include "atlas/util/LonLatMicroDeg.h"
#include "atlas/util/Functions.h"

using atlas::util::Face;
using atlas::util::accumulate_faces;
using atlas::util::Topology;
using atlas::util::UniqueLonLat;
using atlas::util::microdeg;

namespace atlas {
namespace actions {

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

void build_element_to_edge_connectivity( Mesh& mesh )
{

  std::vector< IndexView<int,2> > elem_to_edge( mesh.nb_function_spaces() );
  std::vector< std::vector<int> > edge_cnt( mesh.nb_function_spaces() );

  for( size_t func_space_idx=0; func_space_idx<mesh.nb_function_spaces(); ++func_space_idx)
  {
    FunctionSpace& func_space = mesh.function_space(func_space_idx);
    if( func_space.metadata().get<long>("type") == Entity::ELEMS )
    {
      int nb_edges_per_elem;
      if (func_space.name() == "quads")  nb_edges_per_elem = 4;
      if (func_space.name() == "triags") nb_edges_per_elem = 3;
      elem_to_edge[func_space_idx] =
          IndexView<int,2>(func_space.create_field<int>("to_edge",nb_edges_per_elem));
      elem_to_edge[func_space_idx] = -1;
      edge_cnt[func_space_idx].resize( func_space.shape(0), 0);
    }
  }

  Nodes& nodes = mesh.nodes();
  FunctionSpace& edges = mesh.function_space("edges");

  size_t nb_edges = edges.shape(0);

  IndexView<int,3> edge_to_elem ( edges.field( "to_elem" ).data<int>(), make_shape(nb_edges,2,2) );
  IndexView<int,2> edge_nodes   ( edges.field( "nodes" ) );

  bool has_pole_edges(false);
  ArrayView<int,1> is_pole_edge;
  if( edges.has_field("is_pole_edge") )
  {
    has_pole_edges = true;
    is_pole_edge = ArrayView<int,1>( edges.field("is_pole_edge") );
  }

  UniqueLonLat uid( nodes );

  std::vector<Sort> edge_sort(nb_edges);
  for( size_t edge=0; edge<nb_edges; ++edge )
    edge_sort[edge] = Sort( uid(edge_nodes[edge]), edge );

  std::sort( edge_sort.data(), edge_sort.data()+nb_edges );


//  FunctionSpace& edges = mesh.function_space("edges");
//  int nb_edges = edges.shape(0);
//  IndexView<int,3> edge_to_elem ( edges.field( "to_elem" ).data<int>(), Extents(nb_edges,2,2) );
//  ArrayView<gidx_t,1> edge_gidx    ( edges.field( "glb_idx" ) );
//  bool has_pole_edges(false);
//  ArrayView<int,1> is_pole_edge;
//  if( edges.has_field("is_pole_edge") )
//  {
//    has_pole_edges = true;
//    is_pole_edge = ArrayView<int,1>( edges.field("is_pole_edge") );
//  }


//  std::vector<Sort> edge_sort(nb_edges);
//  for( int edge=0; edge<nb_edges; ++edge )
//    edge_sort[edge] = Sort(edge_gidx(edge),edge);
//  std::sort( edge_sort.data(), edge_sort.data()+nb_edges );


  for( size_t jedge=0; jedge<nb_edges; ++jedge)
  {
    int edge = edge_sort[jedge].i;
    for( size_t j=0; j<2; ++j)
    {
      int func_space_idx = edge_to_elem(edge,j,0);
      int elem           = edge_to_elem(edge,j,1);

      if ( elem >= 0 )
      {
        elem_to_edge[func_space_idx](elem,edge_cnt[func_space_idx][elem]++) = edge;
      }
      else
      {
        if( !( has_pole_edges && is_pole_edge(edge) ) )
        {
          if( func_space_idx >= 0)
            throw eckit::SeriousBug("func_space_idx not negative",Here());
          if( j==0 )
            throw eckit::SeriousBug("edge has no element connected",Here());
        }
      }
    }
  }


	// Verify that all edges have been found
	ASSERT( nb_edges > 0 );
	for( size_t func_space_idx=0; func_space_idx<mesh.nb_function_spaces(); ++func_space_idx)
  {
    FunctionSpace& func_space = mesh.function_space(func_space_idx);
	if( func_space.metadata().get<long>("type") == Entity::ELEMS )
    {
      size_t nb_edges_per_elem;
      if (func_space.name() == "quads")  nb_edges_per_elem = 4;
      if (func_space.name() == "triags") nb_edges_per_elem = 3;
			for( size_t jelem=0; jelem< func_space.shape(0); ++jelem )
			{
				for( size_t jedge=0; jedge<nb_edges_per_elem; ++jedge )
				{
					if( elem_to_edge[func_space_idx](jelem,jedge) < 0 )
					{
						const IndexView<int,2> elem_nodes ( func_space.field("nodes") );
						const ArrayView<gidx_t,1> gidx (nodes.field("glb_idx"));

						std::stringstream msg; msg << "Could not find edge " << jedge << " for " << func_space.name() << " elem " << jelem << " with nodes ( ";
						for( size_t jnode=0; jnode<elem_nodes.shape(1); ++jnode )
						{
							msg << gidx(elem_nodes(jelem,jnode)) <<" ";
						}
						msg << ")";
						throw eckit::SeriousBug(msg.str(),Here());
					}
				}
			}
    }
  }

}

void build_node_to_edge_connectivity( Mesh& mesh )
{
  Nodes& nodes = mesh.nodes();
  FunctionSpace& edges = mesh.function_space("edges");
  int nb_nodes = nodes.shape(0);
  int nb_edges = edges.shape(0);

  IndexView<int,2> edge_nodes   ( edges.field( "nodes" ) );

  // Get max_edge_cnt
  ArrayView<int,1> to_edge_size ( nodes.create_field<int>("to_edge_size",1) );
  to_edge_size = 0.;
  for( int jedge=0; jedge<nb_edges; ++jedge)
  {
    for( int j=0; j<2; ++j)
    {
      ++to_edge_size( edge_nodes(jedge,j) );
    }
  }

  int max_edge_cnt(0);
  for( int jnode=0; jnode<nb_nodes; ++jnode )
  {
    max_edge_cnt = std::max(max_edge_cnt,to_edge_size(jnode));
  }
  ECKIT_MPI_CHECK_RESULT( MPI_Allreduce( MPI_IN_PLACE, &max_edge_cnt, 1, MPI_INT, MPI_MAX, eckit::mpi::comm() ) );

  IndexView<int,2> node_to_edge ( nodes.create_field<int>("to_edge",max_edge_cnt) );

  UniqueLonLat uid( nodes );
  std::vector<Sort> edge_sort(nb_edges);
  for( size_t edge=0; edge<nb_edges; ++edge )
    edge_sort[edge] = Sort( uid(edge_nodes[edge]), edge );
  std::stable_sort( edge_sort.data(), edge_sort.data()+nb_edges );

//  ArrayView<gidx_t,1> edge_gidx    ( edges.field( "glb_idx" ) );
//  std::vector<Sort> edge_sort(nb_edges);
//  for( int edge=0; edge<nb_edges; ++edge )
//    edge_sort[edge] = Sort(edge_gidx(edge),edge);
//  std::sort( edge_sort.data(), edge_sort.data()+nb_edges );

  to_edge_size = 0.;
  for( size_t jedge=0; jedge<nb_edges; ++jedge)
  {
    int edge = edge_sort[jedge].i;
    for( size_t j=0; j<2; ++j)
    {
      int node = edge_nodes(edge,j);
      node_to_edge( node, to_edge_size(node)++ ) = edge;
    }
  }
}


void accumulate_pole_edges( Mesh& mesh, std::vector<int>& pole_edge_nodes, int& nb_pole_edges )
{
  Nodes& nodes   = mesh.nodes();
  ArrayView<double,2> lonlat    ( nodes.lonlat() );
  ArrayView<gidx_t,1> glb_idx   ( nodes.field( "glb_idx"     ) );
  ArrayView<int,   1> part      ( nodes.field( "partition"   ) );
  ArrayView<int,   1> flags     ( nodes.field( "flags"       ) );
  IndexView<int,   1> ridx      ( nodes.field( "remote_idx"  ) );
  size_t nb_nodes = nodes.shape(0);

  double min[2], max[2];
  min[LON] =  std::numeric_limits<double>::max();
  min[LAT] =  std::numeric_limits<double>::max();
  max[LON] = -std::numeric_limits<double>::max();
  max[LAT] = -std::numeric_limits<double>::max();
  for (size_t node=0; node<nb_nodes; ++node)
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

  double tol = 1e-6;
  std::vector<int> north_pole_edges;
  std::vector<int> south_pole_edges;

  std::vector< std::set<int> > pole_nodes(2);

  enum { NORTH=0, SOUTH=1 };

  for (size_t node=0; node<nb_nodes; ++node)
  {
    //std::cout << "node " << node << "   " << std::abs(lonlat(LAT,node)-ymax) << std::endl;

//    // Only add edges that connect non-ghost nodes
//    if( ridx(node) == node && part(node) == eckit::mpi::rank() )
//    {
      if ( std::abs(lonlat(node,LAT)-max[LAT])<tol )
      {
        pole_nodes[NORTH].insert(node);
      }
      else if ( std::abs(lonlat(node,LAT)-min[LAT])<tol )
      {
        pole_nodes[SOUTH].insert(node);
      }
//    }
  }

  // Sanity check
  {
    for( int NS = 0; NS<2; ++NS )
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

  nb_pole_edges = 0;
  for( size_t NS = 0; NS<2; ++NS )
  {
    for( std::set<int>::iterator it=pole_nodes[NS].begin(); it!=pole_nodes[NS].end(); ++it)
    {
      int node = *it;
      if( !Topology::check(flags(node),Topology::PERIODIC|Topology::GHOST) )
      {
        int x1 = microdeg( lonlat(node,LON) );
        int x2 = microdeg( lonlat(node,LON) + 180. );
        for( std::set<int>::iterator itr=pole_nodes[NS].begin(); itr!=pole_nodes[NS].end(); ++itr)
        {
          int other_node = *itr;
          if( microdeg( lonlat(other_node,LON) ) == x2 )
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

// This no longer works as only edges that connect non-ghost nodes are added
//        if ( std::abs( std::abs(lonlat(recip_node,LON) - lonlat(node,LON)) - 180.) > tol )
//        {
//          //std::cout << eckit::mpi::rank() << "  :  distance = " << lonlat(recip_node,LON) - lonlat(node,LON) << std::endl;
//          if( eckit::mpi::rank() == part(node) )
//          {
//            throw eckit::SeriousBug("Not implemented yet, when pole-lattitude is split, "
//                                    "or non-even number of longitudes at pole",Here());
//          }
//          else
//          {
//            // pole is in halo of other partition, and is not completely full to connect edge
//            std::stringstream msg;
//            msg << "Pole is in halo of partition " << eckit::mpi::rank()
//                << " and cannot create edge to connect to other side of pole";
//            throw eckit::SeriousBug(msg.str(),Here());
//          }
//        }
//        else
//        {
//          pole_edge_nodes.push_back(node);
//          pole_edge_nodes.push_back(recip_node);
//          ++nb_pole_edges;
//        }
}


struct ComputeUniquePoleEdgeIndex
{
  ComputeUniquePoleEdgeIndex( const Nodes& nodes )
  {
    lonlat = ArrayView<double,2> ( nodes.lonlat() );
  }

  gidx_t operator()( const IndexView<int,1>& edge_nodes ) const
  {
    double centroid[2];
    centroid[LON] = 0.;
    centroid[LAT] = 0.;
    for( size_t jnode=0; jnode<2; ++jnode )
    {
      centroid[LON] += lonlat( edge_nodes(jnode), LON );
      centroid[LAT] += lonlat( edge_nodes(jnode), LAT );
    }
    centroid[LON] /= 2.;
    centroid[LAT] /= 2.;
    if( centroid[LAT] > 0 )
      centroid[LAT] =  90.;
    else
      centroid[LAT] = -90.;
    /// FIXME make this into `util::unique_lonlat(centroid)` but this causes weird parallel behavior
    return util::detail::unique32( util::microdeg(centroid[LON]), util::microdeg(centroid[LON]) );
  }

  ArrayView<double,2> lonlat;
};

void build_edges( Mesh& mesh )
{
  Nodes& nodes   = mesh.nodes();
  ArrayView<gidx_t,1> glb_idx(        nodes.field( "glb_idx" ) );
  ArrayView<int,1> part   (        nodes.field( "partition" ) );
  ArrayView<double,2> lonlat (     nodes.lonlat() );
  size_t nb_nodes = nodes.shape(0);

  FunctionSpace& quads       = mesh.function_space( "quads" );
  FunctionSpace& triags      = mesh.function_space( "triags" );

  std::vector< std::vector<int> > node_to_face(nb_nodes);
  std::vector< int > face_nodes_data; face_nodes_data.reserve(4*nb_nodes);
  std::vector< Face > face_to_elem;
  face_to_elem.reserve(4*nb_nodes);
  int nb_faces = 0;
  int nb_inner_faces = 0;

  accumulate_faces(quads, node_to_face,face_nodes_data,face_to_elem,nb_faces,nb_inner_faces);
  accumulate_faces(triags,node_to_face,face_nodes_data,face_to_elem,nb_faces,nb_inner_faces);


  size_t extents[] = {nb_faces,2};
  ArrayView<int,2> face_nodes(face_nodes_data.data(),extents);

  // Build edges
  int nb_edges = nb_faces;
  if( ! mesh.has_function_space("edges") )
  {
    mesh.create_function_space( "edges", "shapefunc", make_shape(nb_edges,FunctionSpace::UNDEF_VARS) );
  }
  FunctionSpace& edges = mesh.function_space("edges");
  edges.metadata().set<long>("type",Entity::FACES);
  edges.resize(make_shape(nb_edges,FunctionSpace::UNDEF_VARS));

  if( ! edges.has_field("nodes")      )  edges.create_field<int>("nodes",     2);
  if( ! edges.has_field("glb_idx")    )  edges.create_field<gidx_t>("glb_idx",   1);
  if( ! edges.has_field("partition")  )  edges.create_field<int>("partition", 1);
  if( ! edges.has_field("to_elem")    )  edges.create_field<int>("to_elem",   4);
  if( ! edges.has_field("remote_idx") )  edges.create_field<int>("remote_idx",1);

  IndexView<int,2> edge_nodes   ( edges.field( "nodes"      ) );
  ArrayView<gidx_t,1> edge_glb_idx ( edges.field( "glb_idx"    ) );
  ArrayView<int,1> edge_part    ( edges.field( "partition"  ) );
  IndexView<int,1> edge_ridx    ( edges.field( "remote_idx" ) );
  IndexView<int,3> edge_to_elem ( edges.field( "to_elem"    ).data<int>(), make_shape(nb_edges,2,2) );

  std::vector< IndexView<int,2> > elem_nodes( mesh.nb_function_spaces() );

  for( size_t func_space_idx=0; func_space_idx<mesh.nb_function_spaces(); ++func_space_idx)
  {
    FunctionSpace& func_space = mesh.function_space(func_space_idx);
	if( func_space.metadata().get<long>("type") == Entity::ELEMS )
    {
      elem_nodes[func_space_idx] = IndexView<int,2>(func_space.field("nodes"));
    }
  }

  UniqueLonLat compute_uid( nodes );

  int cnt=0;
  for( size_t edge=0; edge<nb_edges; ++edge )
  {
    const int ip1 = face_nodes(edge,0);
    const int ip2 = face_nodes(edge,1);
    edge_nodes(edge,0) = ip1;
    edge_nodes(edge,1) = ip2;
    //if( glb_idx(ip1) > glb_idx(ip2) )
    if( compute_uid(ip1) > compute_uid(ip2) )
    {
      edge_nodes(edge,0) = ip2;
      edge_nodes(edge,1) = ip1;
    }

    ASSERT( edge_nodes(edge,0) < nb_nodes );
    ASSERT( edge_nodes(edge,1) < nb_nodes );
    edge_glb_idx(edge)   = compute_uid(edge_nodes[edge]);
    edge_part(edge)      = std::min( part(edge_nodes(edge,0)), part(edge_nodes(edge,1) ) );
    edge_ridx(edge)      = edge;

    const int f1 = face_to_elem[edge][0].f;
    const int f2 = face_to_elem[edge][1].f;
    const int e1 = face_to_elem[edge][0].e;
    const int e2 = face_to_elem[edge][1].e;

    edge_to_elem(edge,0,0) = f1;
    edge_to_elem(edge,0,1) = e1;
    edge_to_elem(edge,1,0) = f2;
    edge_to_elem(edge,1,1) = e2;

    if( f2 >= 0 )
    {
      if( compute_uid(elem_nodes[f1][e1]) > compute_uid(elem_nodes[f2][e2]) )
      {
        edge_to_elem(edge,0,0) = f2;
        edge_to_elem(edge,0,1) = e2;
        edge_to_elem(edge,1,0) = f1;
        edge_to_elem(edge,1,1) = e1;
      }
    }
  }
  build_element_to_edge_connectivity(mesh);
}

void build_pole_edges( Mesh& mesh )
{
  Nodes& nodes   = mesh.nodes();

  ArrayView<int,1> part   (        nodes.field( "partition" ) );
  size_t nb_edges = 0;

  if( ! mesh.has_function_space("edges") )
	mesh.create_function_space( "edges","shapefunc", make_shape(nb_edges,FunctionSpace::UNDEF_VARS) );

  FunctionSpace& edges = mesh.function_space("edges");
  edges.metadata().set<long>("type",Entity::FACES);

  nb_edges = edges.shape(0);

  int nb_pole_edges;
  std::vector<int> pole_edge_nodes;
  accumulate_pole_edges( mesh, pole_edge_nodes, nb_pole_edges );
  edges.resize( make_shape(nb_edges+nb_pole_edges, FunctionSpace::UNDEF_VARS) );


  if( ! edges.has_field("nodes")      )    edges.create_field<int>("nodes",     2);
  if( ! edges.has_field("glb_idx")    )    edges.create_field<gidx_t>("glb_idx",   1);
  if( ! edges.has_field("partition")  )    edges.create_field<int>("partition", 1);
  if( ! edges.has_field("to_elem")    )    edges.create_field<int>("to_elem",   4);
  if( ! edges.has_field("remote_idx") )    edges.create_field<int>("remote_idx",1);
  if( ! edges.has_field("is_pole_edge") )  edges.create_field<int>("is_pole_edge",1);

  IndexView<int,2> edge_nodes   ( edges.field( "nodes"      ) );
  ArrayView<gidx_t,1> edge_glb_idx ( edges.field( "glb_idx"    ) );
  ArrayView<int,1> edge_part    ( edges.field( "partition"  ) );
  IndexView<int,1> edge_ridx    ( edges.field( "remote_idx" ) );
  ArrayView<int,1> is_pole_edge ( edges.field( "is_pole_edge" ) );
  IndexView<int,3> edge_to_elem ( edges.field( "to_elem"    ).data<int>(), make_shape(nb_edges+nb_pole_edges,2,2) );

  for(size_t edge=0; edge<nb_edges; ++edge)
  {
    is_pole_edge(edge) = 0;
  }

  size_t cnt = 0;
  ComputeUniquePoleEdgeIndex compute_uid( nodes );
  for(int edge=nb_edges; edge<nb_edges+nb_pole_edges; ++edge)
  {
    edge_nodes(edge,0)   = pole_edge_nodes[cnt++];
    edge_nodes(edge,1)   = pole_edge_nodes[cnt++];
    edge_glb_idx(edge)   = compute_uid( edge_nodes[edge] );
    edge_part(edge)      = std::min( part(edge_nodes(edge,0)), part(edge_nodes(edge,1) ) );
    edge_ridx(edge)      = edge;
    edge_to_elem(edge,0,0) = -1;
    edge_to_elem(edge,0,1) = -1;
    edge_to_elem(edge,1,0) = -1;
    edge_to_elem(edge,1,1) = -1;
    is_pole_edge(edge) = 1;
  }
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
} // namespace atlas

