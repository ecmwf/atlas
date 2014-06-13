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
#include "atlas/mesh/Mesh.hpp"
#include "atlas/mesh/FunctionSpace.hpp"
#include "atlas/mesh/Field.hpp"
#include "atlas/actions/BuildEdges.hpp"
#include "atlas/mesh/Parameters.hpp"
#include "atlas/mesh/ArrayView.hpp"
#include "atlas/mesh/Array.hpp"
#include "atlas/mesh/IndexView.hpp"
#include "atlas/mesh/Util.hpp"

namespace atlas {
namespace actions {

void build_element_to_edge_connectivity( Mesh& mesh, IndexView<int,3>& edge_to_elem )
{
  std::vector< IndexView<int,2> > elem_to_edge( mesh.nb_function_spaces() );
  std::vector< std::vector<int> > edge_cnt( mesh.nb_function_spaces() );

  for( int func_space_idx=0; func_space_idx<mesh.nb_function_spaces(); ++func_space_idx)
  {
    FunctionSpace& func_space = mesh.function_space(func_space_idx);
    if( func_space.metadata<int>("type") == Entity::ELEMS )
    {
      int nb_edges_per_elem;
      if (func_space.name() == "quads")  nb_edges_per_elem = 4;
      if (func_space.name() == "triags") nb_edges_per_elem = 3;
      elem_to_edge[func_space_idx] =
          IndexView<int,2>(func_space.create_field<int>("to_edge",nb_edges_per_elem));
      elem_to_edge[func_space_idx] = -1;
      edge_cnt[func_space_idx].resize( func_space.extents()[0], 0);
    }
  }

  ArrayView<int,1> edge_glb_idx( mesh.function_space("edges").field("glb_idx") );

  int nb_edges = edge_to_elem.extents()[0];
  for( int edge=0; edge<nb_edges; ++edge)
  {
    for( int j=0; j<2; ++j)
    {
      int func_space_idx = edge_to_elem(edge,j,0);
      int elem           = edge_to_elem(edge,j,1);
      if ( elem >= 0 )
      {
        elem_to_edge[func_space_idx](elem,edge_cnt[func_space_idx][elem]++) = edge;
      }
      else
      {
        if( func_space_idx >= 0)
          throw eckit::SeriousBug("func_space_idx not negative",Here());
        if( j==0 )
          throw eckit::SeriousBug("edge has no element connected",Here());
      }
    }
  }
}


void accumulate_pole_edges( Mesh& mesh, std::vector<int>& pole_edge_nodes, int& nb_pole_edges )
{
  FunctionSpace& nodes   = mesh.function_space( "nodes" );
  ArrayView<double,2> coords    ( nodes.field( "coordinates" ) );
  ArrayView<int,   1> glb_idx   ( nodes.field( "glb_idx"     ) );
  ArrayView<int,   1> part      ( nodes.field( "partition"   ) );
  int nb_nodes = nodes.extents()[0];

  double min[2], max[2];
  min[XX] =  std::numeric_limits<double>::max();
  min[YY] =  std::numeric_limits<double>::max();
  max[XX] = -std::numeric_limits<double>::max();
  max[YY] = -std::numeric_limits<double>::max();
  for (int node=0; node<nb_nodes; ++node)
  {
    min[XX] = std::min( min[XX], coords(node,XX) );
    min[YY] = std::min( min[YY], coords(node,YY) );
    max[XX] = std::max( max[XX], coords(node,XX) );
    max[YY] = std::max( max[YY], coords(node,YY) );
  }
  MPL_CHECK_RESULT( MPI_Allreduce( MPI_IN_PLACE, &min[XX], 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD ) );
  MPL_CHECK_RESULT( MPI_Allreduce( MPI_IN_PLACE, &min[YY], 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD ) );
  MPL_CHECK_RESULT( MPI_Allreduce( MPI_IN_PLACE, &max[XX], 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD ) );
  MPL_CHECK_RESULT( MPI_Allreduce( MPI_IN_PLACE, &max[YY], 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD ) );

  double tol = 1e-6;
  std::vector<int> north_pole_edges;
  std::vector<int> south_pole_edges;

  std::vector< std::set<int> > pole_nodes(2);

  enum { NORTH=0, SOUTH=1 };

  for (int node=0; node<nb_nodes; ++node)
  {
    //std::cout << "node " << node << "   " << std::abs(coords(YY,node)-ymax) << std::endl;
    if ( std::abs(coords(node,YY)-max[YY])<tol )
    {
      pole_nodes[NORTH].insert(node);
    }
    else if ( std::abs(coords(node,YY)-min[YY])<tol )
    {
      pole_nodes[SOUTH].insert(node);
    }
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
          throw eckit::NotImplemented(Here());
        }
      }
    }
  }

  nb_pole_edges = 0;
  for( int NS = 0; NS<2; ++NS )
  {
    for( std::set<int>::iterator it=pole_nodes[NS].begin(); it!=pole_nodes[NS].end(); ++it)
    {
      int node = *it;
      double x1 = coords(node,XX);
      double x2 = coords(node,XX) + M_PI;
      double dist = 2.*M_PI;
      if( x1>=min[XX]-tol && x1<=(max[XX]-min[XX])*0.5+tol )
      {
        int recip_node = -1;
        for( std::set<int>::iterator itr=pole_nodes[NS].begin(); itr!=pole_nodes[NS].end(); ++itr)
        {
          int other_node = *itr;
          if( std::abs(coords(other_node,XX)-x2 )<dist )
          {
            dist = std::abs(coords(other_node,XX)-x2);
            recip_node = other_node;
          }
        }
        if ( std::abs( std::abs(coords(recip_node,XX) - coords(node,XX)) - M_PI) > tol )
        {
          //std::cout << MPL::rank() << "  :  distance = " << coords(recip_node,XX) - coords(node,XX) << std::endl;
          if( MPL::rank() == part(node) )
          {
            throw eckit::SeriousBug("Not implemented yet, when pole-lattitude is split, "
                                    "or non-even number of longitudes at pole",Here());
          }
          else
          {
            // pole is in halo of other partition, and is not completely full to connect edge
            std::stringstream msg;
            msg << "Pole is in halo of partition " << MPL::rank()
                << " and cannot create edge to connect to other side of pole";
            throw eckit::SeriousBug(msg.str(),Here());
          }
        }
        else
        {
          pole_edge_nodes.push_back(node);
          pole_edge_nodes.push_back(recip_node);
          ++nb_pole_edges;
        }
      }
    }
  }
}


struct ComputeUniquePoleEdgeIndex
{
  ComputeUniquePoleEdgeIndex( const FunctionSpace& nodes )
  {
    coords = ArrayView<double,2> ( nodes.field("coordinates") );
  }

  int operator()( const IndexView<int,1>& edge_nodes ) const
  {
    double centroid[2];
    centroid[XX] = 0.;
    centroid[YY] = 0.;
    for( int jnode=0; jnode<2; ++jnode )
    {
      centroid[XX] += coords( edge_nodes(jnode), XX );
      centroid[YY] += coords( edge_nodes(jnode), YY );
    }
    centroid[XX] /= 2.;
    centroid[YY] /= 2.;
    if( centroid[YY] > 0 )
      centroid[YY] =  M_PI_2;
    else
      centroid[YY] = -M_PI_2;
    return LatLonPoint( centroid[XX], centroid[YY] ).uid();
  }

  ArrayView<double,2> coords;
};

void build_edges( Mesh& mesh )
{
  FunctionSpace& nodes   = mesh.function_space( "nodes" );
  ArrayView<int,1> glb_idx(        nodes.field( "glb_idx" ) );
  ArrayView<int,1> part   (        nodes.field( "partition" ) );
  int nb_nodes = nodes.extents()[0];

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


  int extents[] = {nb_faces,2};
  ArrayView<int,2> face_nodes(face_nodes_data.data(),extents);

  // Build edges
  int nb_edges = nb_faces;
  if( ! mesh.has_function_space("edges") )
  {
    mesh.add_function_space( new FunctionSpace("edges","shapefunc", Extents(nb_edges,Field::UNDEF_VARS)) );
  }
  FunctionSpace& edges = mesh.function_space("edges");
  edges.metadata().set("type",static_cast<int>(Entity::FACES));
  edges.resize(Extents(nb_edges,Field::UNDEF_VARS));

  if( ! edges.has_field("nodes")      )  edges.create_field<int>("nodes",     2);
  if( ! edges.has_field("glb_idx")    )  edges.create_field<int>("glb_idx",   1);
  if( ! edges.has_field("partition")  )  edges.create_field<int>("partition", 1);
  if( ! edges.has_field("to_elem")    )  edges.create_field<int>("to_elem",   4);
  if( ! edges.has_field("remote_idx") )  edges.create_field<int>("remote_idx",1);

  IndexView<int,2> edge_nodes   ( edges.field( "nodes"      ) );
  ArrayView<int,1> edge_glb_idx ( edges.field( "glb_idx"    ) );
  ArrayView<int,1> edge_part    ( edges.field( "partition"  ) );
  IndexView<int,1> edge_ridx    ( edges.field( "remote_idx" ) );
  IndexView<int,3> edge_to_elem ( edges.field( "to_elem"    ).data<int>(), Extents(nb_edges,2,2) );

  ComputeUniqueElementIndex uid( nodes );

  int cnt=0;
  for( int edge=0; edge<nb_edges; ++edge )
  {
    edge_nodes(edge,0)   = face_nodes(edge,0);
    edge_nodes(edge,1)   = face_nodes(edge,1);
    ASSERT( edge_nodes(edge,0) < nb_nodes );
    ASSERT( edge_nodes(edge,1) < nb_nodes );
    edge_glb_idx(edge)   = uid(edge_nodes[edge]);
    edge_part(edge)      = std::min( glb_idx(edge_nodes(edge,0)), glb_idx(edge_nodes(edge,1) ) );
    edge_ridx(edge)      = edge;
    edge_to_elem(edge,0,0) = face_to_elem[edge][0].f;
    edge_to_elem(edge,0,1) = face_to_elem[edge][0].e;
    edge_to_elem(edge,1,0) = face_to_elem[edge][1].f;
    edge_to_elem(edge,1,1) = face_to_elem[edge][1].e;
  }

  // Element to edge connectivity
  build_element_to_edge_connectivity(mesh,edge_to_elem);

}

void build_pole_edges( Mesh& mesh )
{
  FunctionSpace& nodes   = mesh.function_space( "nodes" );
  ArrayView<int,1> part   (        nodes.field( "partition" ) );
  int nb_edges = 0;
  if( ! mesh.has_function_space("edges") )
    mesh.add_function_space( new FunctionSpace("edges","shapefunc", Extents(nb_edges,Field::UNDEF_VARS)) );
  FunctionSpace& edges = mesh.function_space("edges");
  edges.metadata().set("type",static_cast<int>(Entity::FACES));

  nb_edges = edges.extents()[0];

  int nb_pole_edges;
  std::vector<int> pole_edge_nodes;
  accumulate_pole_edges( mesh, pole_edge_nodes, nb_pole_edges );
  edges.resize( Extents(nb_edges+nb_pole_edges, Field::UNDEF_VARS) );


  if( ! edges.has_field("nodes")      )  edges.create_field<int>("nodes",     2);
  if( ! edges.has_field("glb_idx")    )  edges.create_field<int>("glb_idx",   1);
  if( ! edges.has_field("partition")  )  edges.create_field<int>("partition", 1);
  if( ! edges.has_field("to_elem")    )  edges.create_field<int>("to_elem",   4);
  if( ! edges.has_field("remote_idx") )  edges.create_field<int>("remote_idx",1);

  IndexView<int,2> edge_nodes   ( edges.field( "nodes"      ) );
  ArrayView<int,1> edge_glb_idx ( edges.field( "glb_idx"    ) );
  ArrayView<int,1> edge_part    ( edges.field( "partition"  ) );
  IndexView<int,1> edge_ridx    ( edges.field( "remote_idx" ) );
  IndexView<int,3> edge_to_elem ( edges.field( "to_elem"    ).data<int>(), Extents(nb_edges+nb_pole_edges,2,2) );

  int cnt = 0;
  ComputeUniquePoleEdgeIndex uid( nodes );
  for(int edge=nb_edges; edge<nb_edges+nb_pole_edges; ++edge)
  {
    edge_nodes(edge,0)   = pole_edge_nodes[cnt++];
    edge_nodes(edge,1)   = pole_edge_nodes[cnt++];
    edge_glb_idx(edge)   = uid( edge_nodes[edge] );
    edge_part(edge)      = std::min( part(edge_nodes(edge,0)), part(edge_nodes(edge,1) ) );
    edge_ridx(edge)      = edge;
    edge_to_elem(edge,0,0) = -1;
    edge_to_elem(edge,0,1) = -1;
    edge_to_elem(edge,1,0) = -1;
    edge_to_elem(edge,1,1) = -1;
  }
}


// ------------------------------------------------------------------
// C wrapper interfaces to C++ routines

void atlas__build_edges ( Mesh* mesh) {
  build_edges(*mesh);
}

// ------------------------------------------------------------------

} // namespace actions
} // namespace atlas

