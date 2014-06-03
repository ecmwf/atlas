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
#include "atlas/mesh/Mesh.hpp"
#include "atlas/mesh/FunctionSpace.hpp"
#include "atlas/mesh/Field.hpp"
#include "atlas/actions/BuildHalo.hpp"
#include "atlas/mesh/Parameters.hpp"
#include "atlas/mesh/ArrayView.hpp"
#include "atlas/mesh/IndexView.hpp"
#include "atlas/mesh/Array.hpp"

namespace atlas {

struct Face
{
  ElementRef& operator[](int i) { return elems[i]; }
  bool is_bdry() const { return (elems[1].f < 0); }
  ElementRef elems[2];
};

void scan_function_space(
    FunctionSpace& func_space,
    std::vector< std::vector<int> >& node_to_face,
    std::vector<int>& tmp_edges,
    std::vector< Face >& connectivity_edge_to_elem,
    int& nb_faces,
    int& nb_inner_faces );


void build_halo( Mesh& mesh )
{
  FunctionSpace& nodes         = mesh.function_space( "nodes" );
  ArrayView<double,2> coords(         nodes.field( "coordinates"    ) );
  ArrayView<int   ,1> glb_idx(        nodes.field( "glb_idx"        ) );
  ArrayView<int   ,1> master_glb_idx( nodes.field( "master_glb_idx" ) );
  ArrayView<int   ,1> node_proc(      nodes.field( "proc"           ) );
  int nb_nodes = nodes.extents()[0];

  std::vector< std::vector< ElementRef > > node_to_elem(nb_nodes);
  
  std::vector< ArrayView<int,2> > elem_nodes( mesh.nb_function_spaces() );

  for( int func_space_idx=0; func_space_idx<mesh.nb_function_spaces(); ++func_space_idx)
  {
    FunctionSpace& elements = mesh.function_space(func_space_idx);
    if( elements.metadata<int>("type") == Entity::ELEMS )
    {
      elem_nodes[func_space_idx] = ArrayView<int,2>( elements.field("nodes") );
      int nb_elems = elem_nodes[func_space_idx].extents()[0];
      int nb_nodes_per_elem = elem_nodes[func_space_idx].extents()[1];
      for (int elem=0; elem<nb_elems; ++elem)
      {
        for (int n=0; n<nb_nodes_per_elem; ++n)
        {
          int node = elem_nodes[func_space_idx](elem,n);
          node_to_elem[node].push_back( ElementRef(elements.index(),elem) );
        }
      }
    }
  }
  
  /*
  1) Find nodes at boundary of partition
  2) Communicate glb_index of these boundary nodes to other partitions
  3) Find received glb_index in glb_node_to_local_node list
  4) Find elements in node_to_elem list that belong to me
  5) Make list of all nodes that complete the elements
  6) Communicate elements and nodes back
  7) Adapt mesh
  */


  /*
  1) Find boundary of partition:
  - find unique edges
  - if edge is bdry_edge, then the nodes are bdry nodes
  */
  std::set<int> bdry_nodes_set;
  
  FunctionSpace& quads       = mesh.function_space( "quads" );
  FunctionSpace& triags      = mesh.function_space( "triags" );

  std::vector< std::vector<int> > node_to_face(nb_nodes);
  std::vector< int > faces_nodes_data; faces_nodes_data.reserve(4*nb_nodes);
  std::vector< Face > face_to_elem;
  face_to_elem.reserve(4*nb_nodes);
  int nb_faces = 0;
  int nb_inner_faces = 0;

  scan_function_space(quads, node_to_face,faces_nodes_data,face_to_elem,nb_faces,nb_inner_faces);
  scan_function_space(triags,node_to_face,faces_nodes_data,face_to_elem,nb_faces,nb_inner_faces);
  
  int extents[] = {nb_faces,4};
  ArrayView<int,2> face_nodes(faces_nodes_data.data(),extents);

  for( int jface=0; jface<nb_faces; ++jface )
  {
    if( face_to_elem[jface].is_bdry() )
    {
      for( int jnode=0; jnode<4; ++jnode)
      {
        if( face_nodes(jface,jnode) >= 0 )
        {
          bdry_nodes_set.insert(face_nodes(jface,jnode));
        }
      }
    }
  }
  std::vector<int> bdry_nodes( bdry_nodes.begin(), bdry_nodes.end());

  /*
  2) Communicate glb_index of these boundary nodes to other partitions
  3) Find received glb_index in glb_node_to_local_node list
  */

  std::map<int,int> node_glb_to_loc;
  for( int jnode=0; jnode<nb_nodes; ++jnode )
  {
    node_glb_to_loc[ glb_idx(jnode) ] = jnode;
  }

  int nb_bdry_nodes = bdry_nodes.size();
  std::vector<int> bdry_nodes_glb_idx( nb_bdry_nodes );
  for( int jnode=0; jnode<nb_bdry_nodes; ++jnode )
  {
    bdry_nodes_glb_idx[jnode] = glb_idx[bdry_nodes[jnode]];
  }

  std::vector< int > par_nb_bdry_nodes( MPL::size() );
  MPL_CHECK_RESULT( MPI_Allgather(&nb_bdry_nodes,            1, MPI_INT,
                                   par_nb_bdry_nodes.data(), 1, MPI_INT, MPI_COMM_WORLD ) );

  int recvcnt=0;
  std::vector<int> recvcounts( MPL::size() );
  std::vector<int> recvdispls( MPL::size() );
  recvdispls[0] = 0;
  recvcounts[0] = par_nb_bdry_nodes[0];
  for( int jproc=1; jproc<MPL::size(); ++jproc )
  {
    recvcounts[jproc] = par_nb_bdry_nodes[jproc];
    recvdispls[jproc] = recvdispls[jproc-1] + recvcounts[jproc-1];
    recvcnt += recvcounts[jproc];
  }
  std::vector<int> recvbuf(recvcnt);
  MPL_CHECK_RESULT( MPI_Allgatherv( bdry_nodes_glb_idx.data(), nb_bdry_nodes, MPI_INT,
                    recvbuf.data(), recvcounts.data(), recvdispls.data(),
                    MPI_INT, MPI_COMM_WORLD) );

  std::vector< std::vector<int> > recv_bdry_nodes_glb_idx( MPL::size() );

  for (int jproc=0; jproc<MPL::size(); ++jproc)
  {
    int* recv_bdry_nodes_glb_idx = recvbuf.data()+recvdispls[jproc];
    int recv_nb_bdry_nodes = par_nb_bdry_nodes[jproc];

    // Find elements that have these nodes
    // In order to do this, check the node_to_elem list

    std::set<ElementRef> send_bdry_elements_set;
    if( jproc != MPL::rank() )
    {
      for( int jrecv=0; jrecv<recv_nb_bdry_nodes; ++jrecv )
      {
        int recv_glb_idx = recv_bdry_nodes_glb_idx[jrecv];
        if( node_glb_to_loc.count(recv_glb_idx) )
        {
          int loc = node_glb_to_loc[recv_glb_idx];
          send_bdry_elements_set.insert( node_to_elem[loc].begin(), node_to_elem[loc].end() );
        }
      }
    }
    std::vector<ElementRef> send_bdry_elements(send_bdry_elements_set.begin(), send_bdry_elements_set.end() );
    int nb_bdry_elems=send_bdry_elements.size();

    // Collect all nodes needed to complete the element
    std::set<int> send_bdry_nodes_glb_idx_set;
    if( jproc != MPL::rank() )
    {
      for( int jelem=0; jelem<nb_bdry_elems; ++jelem )
      {
        int f = send_bdry_elements[jelem].f;
        int e = send_bdry_elements[jelem].e;
        int nb_elem_nodes = elem_nodes[f].extents()[1];
        for( int n=0; n<nb_elem_nodes; ++n )
        {
          send_bdry_nodes_glb_idx_set.insert( glb_idx[ elem_nodes[f](e,n) ] );
        }
      }
      // Remove nodes we already have received, as we won't need to send them
      for( int jrecv=0; jrecv<recv_nb_bdry_nodes; ++jrecv)
      {
        send_bdry_nodes_glb_idx_set.erase(recv_bdry_nodes_glb_idx[jrecv]);
      }
    }
    std::vector<int> send_bdry_nodes_glb_idx(send_bdry_nodes_glb_idx_set.begin(),send_bdry_nodes_glb_idx_set.end());
    int nb_send_bdry_nodes = send_bdry_nodes_glb_idx.size();

    std::vector<int> send_bdry_nodes_loc_idx( nb_send_bdry_nodes );
    for( int jnode=0; jnode<nb_send_bdry_nodes; ++jnode)
    {
      send_bdry_nodes_loc_idx[jnode] = node_glb_to_loc[send_bdry_nodes_glb_idx[jnode]];
    }
  }

}

void scan_function_space(
    FunctionSpace& func_space,
    std::vector< std::vector<int> >& node_to_face,
    std::vector<int>& tmp_edges,
    std::vector< Face >& connectivity_edge_to_elem,
    int& nb_faces,
    int& nb_inner_faces )
{
  IndexView<int,2> elem_nodes( func_space.field( "nodes" ) );
  int nb_elems = func_space.extents()[0];
  int nb_nodes_in_face = 2;

  std::vector< std::vector<int> > face_node_numbering;
  int nb_faces_in_elem;
  if (func_space.name() == "quads")
  {
    nb_faces_in_elem = 4;
    face_node_numbering.resize(nb_faces_in_elem, std::vector<int>(nb_nodes_in_face) );
    face_node_numbering[0][0] = 0;
    face_node_numbering[0][1] = 1;
    face_node_numbering[1][0] = 1;
    face_node_numbering[1][1] = 2;
    face_node_numbering[2][0] = 2;
    face_node_numbering[2][1] = 3;
    face_node_numbering[3][0] = 3;
    face_node_numbering[3][1] = 0;
  }
  else if (func_space.name() == "triags")
  {
    nb_faces_in_elem = 3;
    face_node_numbering.resize(nb_faces_in_elem, std::vector<int>(nb_nodes_in_face) );
    face_node_numbering[0][0] = 0;
    face_node_numbering[0][1] = 1;
    face_node_numbering[1][0] = 1;
    face_node_numbering[1][1] = 2;
    face_node_numbering[2][0] = 2;
    face_node_numbering[2][1] = 0;
  }
  else
  {
    throw std::runtime_error(func_space.name()+" is not \"quads\" or \"triags\"");
  }

  for (int e=0; e<nb_elems; ++e)
  {
    for (int f=0; f<nb_faces_in_elem; ++f)
    {
      bool found_face = false;

      std::vector<int> face_nodes(nb_nodes_in_face);
      for (int jnode=0; jnode<nb_nodes_in_face; ++jnode)
        face_nodes[jnode] = elem_nodes(e,face_node_numbering[f][jnode]);

      int node = face_nodes[0];
      for( int jface=0; jface< node_to_face[node].size(); ++jface )
      {
        int face = node_to_face[node][jface];
        int nb_matched_nodes = 0;
        if (nb_nodes_in_face>1) // 2D or 3D
        {
          for( int jnode=0; jnode<nb_nodes_in_face; ++jnode)
          {
            int other_node = face_nodes[jnode];
            for( int iface=0; iface<node_to_face[other_node].size(); ++iface )
            {
              if( node_to_face[face_nodes[jnode]][iface] == face )
              {
                ++nb_matched_nodes;
                break;
              }
            }
          }
          if (nb_matched_nodes == nb_nodes_in_face)
          {
            connectivity_edge_to_elem[face][1].f = func_space.index();
            connectivity_edge_to_elem[face][1].e = e;
            ++nb_inner_faces;
            found_face = true;
            break;
          }
        }
      }

      if (found_face == false)
      {
        connectivity_edge_to_elem.push_back( Face() );
        connectivity_edge_to_elem[nb_faces][0].f = func_space.index();
        connectivity_edge_to_elem[nb_faces][0].e = e;
        // if 3rd or 4th element stays negative, it is a bdry face
        connectivity_edge_to_elem[nb_faces][0].f = -1;
        connectivity_edge_to_elem[nb_faces][0].e = -1;
        for (int n=0; n<nb_nodes_in_face; ++n)
        {
          node_to_face[face_nodes[n]].push_back(nb_faces);
          tmp_edges.push_back(face_nodes[n]);
        }
        ++nb_faces;
      }
    }
  }
}

// ------------------------------------------------------------------
// C wrapper interfaces to C++ routines

void atlas__build_halo ( Mesh* mesh) {
  build_halo(*mesh);
}

// ------------------------------------------------------------------


} // namespace atlas

