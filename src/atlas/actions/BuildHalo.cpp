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

namespace atlas {

void build_halo( Mesh& mesh )
{
  FunctionSpace& nodes         = mesh.function_space( "nodes" );
  ArrayView<double,2> coords(         nodes.field( "coordinates"    ) );
  ArrayView<int   ,1> glb_idx(        nodes.field( "glb_idx"        ) );
  ArrayView<int   ,1> master_glb_idx( nodes.field( "master_glb_idx" ) );
  ArrayView<int   ,1> node_proc(      nodes.field( "proc"           ) );
  int nb_nodes = coords.extents()[0];

  std::vector< std::vector< ElementRef > > node_to_elem(nb_nodes);
  
  for( int func_space_idx=0; func_space_idx<mesh.nb_function_spaces(); ++func_space_idx)
  {
    FunctionSpace& elements = mesh.function_space(func_space_idx);
    if( elements.metadata<int>("type") == Entity::ELEMS )
    {
      ArrayView<int,2> elem_nodes( elements.field("nodes") );
      ArrayView<int,1> elem_proc(  elements.field("proc" ) );
      int nb_elems = elem_nodes.extents()[0];
      int nb_nodes_per_elem = elem_nodes.extents()[1];
      for (int elem=0; elem<nb_elems; ++elem)
      {
        for (int n=0; n<nb_nodes_per_elem; ++n)
        {
          int node = elem_nodes(elem,n);
          node_to_elem[node].push_back( ElementRef(elements.index(),elem) );
        }
      }
    }
  }
  
  /*
  1) Find nodes at boundary of partition
  2) Communicate glb_index of this node to other partitions
  3) Find received glb_index in glb_node_to_local_node list
  4) Find elements in node_to_elem list that belong to me
  5) Make list of all nodes that complete the elements
  6) Communicate elements and nodes back
  7) Adapt mesh
  */
  
  
  
  
  
  
  
  
  
  

#if 0
  std::vector< FieldT<int>* > elem_nodes(mesh.nb_function_spaces());
  std::vector< int > nb_nodes_per_elem(mesh.nb_function_spaces(),0);
  std::vector< int > nb_elems(mesh.nb_function_spaces(),0);
  for (int f=0; f<mesh.nb_function_spaces(); ++f)
  {
    FunctionSpace& func_space = mesh.function_space(f);
    if (func_space.metadata<int>("type") == Entity::ELEMS)
    {
      elem_nodes[f] = &func_space.field<int>("nodes");
      nb_nodes_per_elem[f] = elem_nodes[f]->bounds()[0];
      nb_elems[f] = elem_nodes[f]->bounds()[0];
    }
  }
  std::vector<int> node_parts(nb_nodes);
  for (int node=0; node<nb_nodes; ++node)
  {
    node_parts[node] = node_proc(node);
  }

  int nb_parts = 1;

  // Loop parts
  for (int mypart=0; mypart<nb_parts; ++mypart)
  {
    // Loop elements
    for (int f=0; f<mesh.nb_function_spaces(); ++f)
    {
      for (int elem=0; elem<nb_elems[f]; ++elem)
      {
        // Loop nodes of element
        for (int n=0; n<nb_nodes_per_elem[f]; ++n)
        {
          int node = (*elem_nodes[f])(n,elem);
          if (node_parts[node] == mypart)
          {
            // Loop neighbour elements of node
            for (int jelem=0; node_to_elem[node].size(); ++jelem)
            {
              int neighbour_func_space_idx=node_to_elem[node][jelem].f;
              int neighbour_elem=node_to_elem[node][jelem].e;
              if (neighbour_elem!=elem || neighbour_func_space_idx!=f)
              {
                // Found "neighbour_elem" of element "elem"
                for (int jnode=0; jnode<nb_nodes_per_elem[neighbour_func_space_idx]; ++jnode)
                {
                  int neighbour_node = (*elem_nodes[neighbour_func_space_idx])(jnode,neighbour_elem);
                  // add this node to my part, essentially marking it as ghost
                  node_parts[neighbour_node] = mypart;
                }
              }
            }
          }
        }
      }
    }
  }
  #endif
}

// ------------------------------------------------------------------
// C wrapper interfaces to C++ routines

void atlas__build_halo ( Mesh* mesh) {
  build_halo(*mesh);
}

// ------------------------------------------------------------------


} // namespace atlas

