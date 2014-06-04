/*
 * (C) Copyright 1996-2014 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */


/// @warning Still doesn't know about periodic BC to enlarge Halo

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

void accumulate_faces(
    FunctionSpace& func_space,
    std::vector< std::vector<int> >& node_to_face,
    std::vector<int>& face_nodes_data,
    std::vector< Face >& connectivity_edge_to_elem,
    int& nb_faces,
    int& nb_inner_faces );


void build_halo( Mesh& mesh )
{
  FunctionSpace& nodes         = mesh.function_space( "nodes" );
  ArrayView<double,2> latlon(         nodes.field( "coordinates"    ) );
  ArrayView<int   ,1> glb_idx(        nodes.field( "glb_idx"        ) );
  ArrayView<int   ,1> mglb_idx( nodes.field( "master_glb_idx" ) );
  ArrayView<int   ,1> node_proc(      nodes.field( "proc"           ) );
  int nb_nodes = nodes.extents()[0];

  std::vector< std::vector< ElementRef > > node_to_elem(nb_nodes);
  
  std::vector< ArrayView<int,2> > elem_nodes( mesh.nb_function_spaces() );
  std::vector< ArrayView<int,1> > elem_proc ( mesh.nb_function_spaces() );
  std::vector< ArrayView<int,1> > elem_glb_idx ( mesh.nb_function_spaces() );

  for( int func_space_idx=0; func_space_idx<mesh.nb_function_spaces(); ++func_space_idx)
  {
    FunctionSpace& elements = mesh.function_space(func_space_idx);
    if( elements.metadata<int>("type") == Entity::ELEMS )
    {
      elem_nodes[func_space_idx] = ArrayView<int,2>( elements.field("nodes") );
      elem_proc [func_space_idx] = ArrayView<int,1>( elements.field("proc") );
      elem_glb_idx [func_space_idx] = ArrayView<int,1>( elements.field("glb_idx") );
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

  accumulate_faces(quads, node_to_face,faces_nodes_data,face_to_elem,nb_faces,nb_inner_faces);
  accumulate_faces(triags,node_to_face,faces_nodes_data,face_to_elem,nb_faces,nb_inner_faces);
  
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

  std::vector<int> recvcounts( MPL::size() );
  std::vector<int> recvdispls( MPL::size() );
  recvdispls[0] = 0;
  recvcounts[0] = par_nb_bdry_nodes[0];
  int recvcnt = recvcounts[0];
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


  // sfn stands for "send_found_nodes"
  std::vector< std::vector<int>    > sfn_proc( MPL::size() );
  std::vector< std::vector<int>    > sfn_mglb_idx( MPL::size() );
  std::vector< std::vector<int>    > sfn_glb_idx( MPL::size() );
  std::vector< std::vector<double> > sfn_latlon ( MPL::size() );
  // sfn stands for "send_found_elems"
  std::vector< std::vector< std::vector<int> > > sfe_glb_idx( mesh.nb_function_spaces(), std::vector< std::vector<int> >( MPL::size() ) );
  std::vector< std::vector< std::vector<int> > > sfe_nodes  ( mesh.nb_function_spaces(), std::vector< std::vector<int> >( MPL::size() ) );

  for (int jproc=0; jproc<MPL::size(); ++jproc)
  {
    int* recv_bdry_nodes_glb_idx = recvbuf.data()+recvdispls[jproc];
    int recv_nb_bdry_nodes = par_nb_bdry_nodes[jproc];

    // Find elements that have these nodes
    // In order to do this, check the node_to_elem list
    // Warning: only add elements that are owned!

    std::vector< std::set<int> > found_bdry_elements_set( mesh.nb_function_spaces() );
    if( jproc != MPL::rank() )
    {
      for( int jrecv=0; jrecv<recv_nb_bdry_nodes; ++jrecv )
      {
        int recv_glb_idx = recv_bdry_nodes_glb_idx[jrecv];
        if( node_glb_to_loc.count(recv_glb_idx) )
        {
          int loc = node_glb_to_loc[recv_glb_idx];
          for( int jelem=0; jelem<node_to_elem[loc].size(); ++jelem )
          {
            int f = node_to_elem[loc][jelem].f;
            int e = node_to_elem[loc][jelem].e;
            if( elem_proc[f](e) == MPL::rank() )
            {
              found_bdry_elements_set[f].insert( e );
            }
          }
        }
      }
    }
    std::vector< std::vector<int> > found_bdry_elements( mesh.nb_function_spaces() );
    std::vector<int> nb_found_bdry_elems( mesh.nb_function_spaces() );
    for( int f=0; f<mesh.nb_function_spaces(); ++f )
    {
      found_bdry_elements[f] = std::vector<int>( found_bdry_elements_set[f].begin(), found_bdry_elements_set[f].end() );
      nb_found_bdry_elems[f] = found_bdry_elements[f].size();
    }

    // Collect all nodes needed to complete the element
    std::set<int> found_bdry_nodes_glb_idx_set;
    if( jproc != MPL::rank() )
    {
      for( int f=0; f<mesh.nb_function_spaces(); ++f )
      {
        for( int jelem=0; jelem<nb_found_bdry_elems[f]; ++jelem )
        {
          int e = found_bdry_elements[f][jelem];
          int nb_elem_nodes = elem_nodes[f].extents()[1];
          for( int n=0; n<nb_elem_nodes; ++n )
          {
            found_bdry_nodes_glb_idx_set.insert( glb_idx[ elem_nodes[f](e,n) ] );
          }
        }
      }
      // Remove nodes we already have received, as we won't need to send them
      for( int jrecv=0; jrecv<recv_nb_bdry_nodes; ++jrecv)
      {
        found_bdry_nodes_glb_idx_set.erase(recv_bdry_nodes_glb_idx[jrecv]);
      }
    }
    int nb_found_bdry_nodes = found_bdry_nodes_glb_idx_set.size();

    sfn_glb_idx[jproc] = std::vector<int>(found_bdry_nodes_glb_idx_set.begin(),found_bdry_nodes_glb_idx_set.end());
    sfn_mglb_idx.resize(nb_found_bdry_nodes);
    sfn_proc.resize(nb_found_bdry_nodes);
    sfn_latlon.resize(2*nb_found_bdry_nodes);
    for( int jnode=0; jnode<nb_found_bdry_nodes; ++jnode)
    {
      int loc = node_glb_to_loc[sfn_glb_idx[jproc][jnode]];
      sfn_mglb_idx[jproc][jnode] = mglb_idx[loc];
      sfn_proc[jproc][jnode] = node_proc[loc];
      sfn_latlon[jproc][jnode*2+XX] = latlon(loc,XX);
      sfn_latlon[jproc][jnode*2+YY] = latlon(loc,YY);
    }

    for( int f=0; f<mesh.nb_function_spaces(); ++f )
    {
      int nb_elem_nodes = elem_nodes[f].extents()[1];
      sfe_glb_idx[f][jproc].resize( nb_found_bdry_elems[f] );
      sfe_nodes[f][jproc].resize( nb_found_bdry_elems[f]*nb_elem_nodes );
      ArrayView<int,2> sfe_nodes_view( sfe_nodes[f][jproc].data(), Extents(nb_found_bdry_elems[f],nb_elem_nodes).data() );
      for( int jelem=0; jelem<nb_found_bdry_elems[f]; ++jelem )
      {
        int e = found_bdry_elements[f][jelem];
        sfe_glb_idx[f][jproc][jelem] = elem_glb_idx[f](e);
        for( int n=0; n<nb_elem_nodes; ++n)
        {
          sfe_nodes_view(jelem,n) = elem_nodes[f](e,n);
        }
      }
    }
  }

  // Now communicate all found fields back

  //    rfn stands for "recv_found_nodes"
  std::vector< std::vector<int> >    rfn_glb_idx(MPL::size());
  std::vector< std::vector<int> >    rfn_mglb_idx(MPL::size());
  std::vector< std::vector<int> >    rfn_proc(MPL::size());
  std::vector< std::vector<double> > rfn_latlon(MPL::size());
  //    rfe stands for "recv_found_elems"
  std::vector< std::vector< std::vector<int> > > rfe_glb_idx( mesh.nb_function_spaces(), std::vector< std::vector<int> >( MPL::size() ) );
  std::vector< std::vector< std::vector<int> > > rfe_nodes  ( mesh.nb_function_spaces(), std::vector< std::vector<int> >( MPL::size() ) );

  MPL::Alltoall(sfn_glb_idx,  rfn_glb_idx);
  MPL::Alltoall(sfn_mglb_idx, rfn_mglb_idx);
  MPL::Alltoall(sfn_proc,     rfn_proc);
  MPL::Alltoall(sfn_latlon,   rfn_latlon);

  for( int f=0; f<mesh.nb_function_spaces(); ++f )
  {
    MPL::Alltoall(sfe_glb_idx[f], sfe_glb_idx[f] );
    MPL::Alltoall(sfe_nodes[f],   sfe_nodes[f]   );
  }

  // We now have everything we need in rfe_ and rfn_ vectors
  // Now adapt the mesh

  // Nodes might be duplicated from different Tasks. We need to identify unique entries
  std::set<int> rfn_glb_idx_unique;
  std::vector< std::vector<int> > rfn_unique_idx(MPL::size());
  for( int jproc=0; jproc<MPL::size(); ++jproc )
  {
    rfn_unique_idx[jproc].reserve(rfn_glb_idx[jproc].size());
  }

  int nb_new_nodes=0;
  for( int jproc=0; jproc<MPL::size(); ++jproc )
  {
    for( int n=0; n<rfn_glb_idx[jproc].size(); ++n )
    {
      bool inserted = rfn_glb_idx_unique.insert( rfn_glb_idx[jproc][n] ).second;
      if( inserted )
      {
        rfn_unique_idx[jproc].push_back(n);
      }
    }
    nb_new_nodes += rfn_unique_idx[jproc].size();
  }

  nodes.resize( Extents( nb_nodes+nb_new_nodes, Field::UNDEF_VARS ) );
  glb_idx   = ArrayView<int,   1>( nodes.field("glb_idx") );
  mglb_idx  = ArrayView<int,   1>( nodes.field("master_glb_idx") );
  node_proc = ArrayView<int,   1>( nodes.field("proc") );
  latlon    = ArrayView<double,2>( nodes.field("coordinates") );

  int new_node=0;
  for( int jproc=0; jproc<MPL::size(); ++jproc )
  {
    for( int n=0; n<rfn_unique_idx[jproc].size(); ++n )
    {
      glb_idx  (nb_nodes+new_node)    = rfn_glb_idx [jproc][rfn_unique_idx[jproc][n]];
      mglb_idx (nb_nodes+new_node)    = rfn_mglb_idx[jproc][rfn_unique_idx[jproc][n]];
      node_proc(nb_nodes+new_node)    = rfn_proc    [jproc][rfn_unique_idx[jproc][n]];
      latlon   (nb_nodes+new_node,XX) = rfn_latlon  [jproc][rfn_unique_idx[jproc][n]*2+XX];
      latlon   (nb_nodes+new_node,XX) = rfn_latlon  [jproc][rfn_unique_idx[jproc][n]*2+YY];
      ++new_node;
    }
  }

  for( int f=0; f<mesh.nb_function_spaces(); ++f )
  {
    FunctionSpace& elements = mesh.function_space(f);
    if( elements.metadata<int>("type") == Entity::ELEMS )
    {
      // Elements might be duplicated from different Tasks. We need to identify unique entries
      std::set<int> rfe_glb_idx_unique;
      std::vector< std::vector<int> > rfe_unique_idx(MPL::size());
      for( int jproc=0; jproc<MPL::size(); ++jproc )
      {
        rfe_unique_idx[jproc].reserve(rfe_glb_idx[f][jproc].size());
      }

      int nb_new_elems=0;
      for( int jproc=0; jproc<MPL::size(); ++jproc )
      {
        for( int e=0; e<rfe_glb_idx[f][jproc].size(); ++e )
        {
          bool inserted = rfe_glb_idx_unique.insert( rfe_glb_idx[f][jproc][e] ).second;
          if( inserted )
          {
            rfe_unique_idx[jproc].push_back(e);
          }
        }
        nb_new_elems += rfe_unique_idx[jproc].size();
      }
      int nb_elems /*old*/ = elements.extents()[0];
      int nb_nodes_per_elem = elem_nodes[f].extents()[1];
      elements.resize( Extents( nb_elems+nb_new_elems, Field::UNDEF_VARS ) );
      elem_glb_idx[f] = ArrayView<int,1>( nodes.field("glb_idx") );
      elem_nodes[f]   = ArrayView<int,2>( nodes.field("nodes")   );
      elem_proc[f]    = ArrayView<int,1>( nodes.field("proc")   );
      int new_elem=0;
      for( int jproc=0; jproc<MPL::size(); ++jproc )
      {
        for( int e=0; e<rfe_unique_idx[jproc].size(); ++e )
        {
          elem_glb_idx[f](nb_elems+new_elem)   = rfe_glb_idx[f][jproc][rfe_unique_idx[jproc][e]];
          elem_proc   [f](nb_elems+new_elem)   = jproc;
          for( int n=0; n<nb_nodes_per_elem; ++n )
            elem_nodes[f](nb_elems+new_elem,n) = rfe_nodes [f][jproc][rfe_unique_idx[jproc][e]*nb_nodes_per_elem+n];
          ++new_elem;
        }
      }
    }
  }
}

void accumulate_faces(
    FunctionSpace& func_space,
    std::vector< std::vector<int> >& node_to_face,
    std::vector<int>& face_nodes_data,
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
        // if 2nd element stays negative, it is a bdry face
        connectivity_edge_to_elem[nb_faces][1].f = -1;
        connectivity_edge_to_elem[nb_faces][1].e = -1;
        for (int n=0; n<nb_nodes_in_face; ++n)
        {
          node_to_face[face_nodes[n]].push_back(nb_faces);
          face_nodes_data.push_back(face_nodes[n]);
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

