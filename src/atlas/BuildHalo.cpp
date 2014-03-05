#include <iostream>
#include <stdexcept>
#include <cmath>
#include "Mesh.hpp"
#include "FunctionSpace.hpp"
#include "Field.hpp"
#include "BuildHalo.hpp"
#include "Parameters.hpp"

namespace atlas {

void build_halo( Mesh& mesh )
{
  FunctionSpace& nodes_2d   = mesh.function_space( "nodes_2d" );
  FieldT<double>& coords    = nodes_2d.field<double>( "coordinates" );
  FieldT<int>& glb_idx      = nodes_2d.field<int>( "glb_idx" );
  FieldT<int>& master_glb_idx  = nodes_2d.field<int>( "master_glb_idx" );
  FieldT<int>& node_proc         = nodes_2d.field<int>( "proc" );
  int nb_nodes = nodes_2d.bounds()[1];

  FunctionSpace& edges   = mesh.function_space( "edges" );
  FieldT<int>& edge_to_elem   = edges.field<int>( "to_elem" );
  FieldT<int>& edge_nodes   = edges.field<int>( "nodes" );
  int nb_edges = edges.bounds()[1];

  std::vector< std::vector< ElementRef > > node_to_elem(nb_nodes);

  {
    FunctionSpace& elements = mesh.function_space("quads");
    FieldT<int>& elem_nodes = elements.field<int>("nodes");
    FieldT<int>& elem_proc = elements.field<int>("proc");
    int nb_elems = elements.bounds()[1];
    int nb_nodes_per_elem = elem_nodes.bounds()[0];
    for (int elem=0; elem<nb_elems; ++elem)
    {
      for (int n=0; n<nb_nodes_per_elem; ++n)
      {
        int node = elem_nodes(n,elem);
        node_to_elem[node].push_back( ElementRef(elements.index(),elem) );
      }
    }
  }

  std::vector< FieldT<int>* > elem_nodes(mesh.nb_function_spaces());
  std::vector< int > nb_nodes_per_elem(mesh.nb_function_spaces(),0);
  std::vector< int > nb_elems(mesh.nb_function_spaces(),0);
  for (int f=0; f<mesh.nb_function_spaces(); ++f)
  {
    FunctionSpace& func_space = mesh.function_space(f);
    if (func_space.metadata<int>("type") == ELEMS)
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
}

// ------------------------------------------------------------------
// C wrapper interfaces to C++ routines

void atlas__build_halo ( Mesh* mesh) {
  build_halo(*mesh);
}

// ------------------------------------------------------------------


} // namespace atlas

