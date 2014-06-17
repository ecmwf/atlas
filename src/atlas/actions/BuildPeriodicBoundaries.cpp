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
#include <set>

#include "atlas/mesh/Mesh.hpp"
#include "atlas/mesh/FunctionSpace.hpp"
#include "atlas/mesh/Field.hpp"
#include "atlas/actions/BuildPeriodicBoundaries.hpp"
#include "atlas/mesh/Parameters.hpp"
#include "atlas/util/ArrayView.hpp"
#include "atlas/util/IndexView.hpp"

namespace atlas {

void scan_bdry_elements( FunctionSpace& elements, FunctionSpace& nodes,
                         double min[2], double max[2],
                         std::vector< ElementRef >& east_bdry_elements,
                         std::vector< ElementRef >& west_bdry_elements )
{
  ArrayView<double,2> coords         ( nodes.field("coordinates"   ) );
  ArrayView<int,   1> glb_idx        ( nodes.field("glb_idx"       ) );
  ArrayView<int,   1> master_glb_idx ( nodes.field("master_glb_idx") );
  IndexView<int,   2> elem_nodes     ( elements.field("nodes"      ) );

  int nb_elems = elements.extents()[0];
  int nb_nodes_per_elem = elem_nodes.extents()[1];
  for (int elem=0; elem<nb_elems; ++elem)
  {
    for (int n=0; n<nb_nodes_per_elem; ++n)
    {
      int node = elem_nodes(elem,n);
      if ( coords(node,XX) == max[XX] )
      {
        east_bdry_elements.push_back( ElementRef(elements.index(),elem) );
        break;
      }
      else if ( coords(node,XX) == min[XX] )
      {
        west_bdry_elements.push_back( ElementRef(elements.index(),elem) );
        break;
      }
    }
  }
}


void build_periodic_boundaries( Mesh& mesh )
{
  FunctionSpace& nodes   = mesh.function_space( "nodes" );
  ArrayView<double,2> coords         ( nodes.field("coordinates"   ) );
  ArrayView<int,   1> glb_idx        ( nodes.field("glb_idx"       ) );
  ArrayView<int,   1> master_glb_idx ( nodes.field("master_glb_idx") );
  ArrayView<int,   1> proc           ( nodes.field("proc"          ) );
  
  int nb_nodes = nodes.extents()[0];


  double min[2];
  double max[2];
  max[XX] = -1e10;
  max[YY] = -1e10;
  min[XX] =  1e10;
  min[YY] =  1e10;
  for (int node=0; node<nb_nodes; ++node)
  {
    min[XX] = std::min( min[XX], coords(node,XX) );
    min[YY] = std::min( min[YY], coords(node,YY) );
    max[XX] = std::max( max[XX], coords(node,XX) );
    max[YY] = std::max( max[YY], coords(node,YY) );
  }
  nodes.metadata().set("xmin",min[XX]);
  nodes.metadata().set("xmax",max[XX]);
  nodes.metadata().set("ymin",min[YY]);
  nodes.metadata().set("ymax",max[YY]);

  double dx = max[XX]-min[XX];

  // Find east boundary elements (that connect to a boundary node)
  std::vector< ElementRef > east_bdry_elements;
  std::vector< ElementRef > west_bdry_elements;
  for (int f=0; f<mesh.nb_function_spaces(); ++f)
  {
    if (mesh.function_space(f).metadata<int>("type") == Entity::ELEMS)
      scan_bdry_elements( mesh.function_space(f),  nodes, min, max,
                          east_bdry_elements,
                          west_bdry_elements);

  }
  int nb_east_bdry_elements = east_bdry_elements.size();
  int nb_west_bdry_elements = west_bdry_elements.size();

  FunctionSpace& quads = mesh.function_space("quads");
  FunctionSpace& triags = mesh.function_space("triags");
  std::vector<int> nb_elems(mesh.nb_function_spaces(),0);
  std::vector<int> nb_nodes_per_elem(mesh.nb_function_spaces(),0);
  std::vector< IndexView<int,2> > elem_nodes(mesh.nb_function_spaces());

  int elements_max_glb_idx=0;
  for (int f=0; f<mesh.nb_function_spaces(); ++f)
  {
    if (mesh.function_space(f).metadata<int>("type") == Entity::ELEMS)
    {
      elem_nodes[f] = IndexView<int,2>(mesh.function_space(f).field("nodes"));
      nb_elems[f]   = elem_nodes[f].extents()[0];
      nb_nodes_per_elem[f] = elem_nodes[f].extents()[1];
      elements_max_glb_idx = std::max(elements_max_glb_idx, mesh.function_space(f).metadata<int>("max_glb_idx"));
    }
  }

  for (int elem=0; elem<nb_east_bdry_elements; ++elem)
  {
    int func_space_idx = east_bdry_elements[elem].f;
    ++nb_elems[func_space_idx];
  }


  double tol = 1.e-6;
  std::set<int> new_nodes;
  std::set<int> transform_to_ghost;
  for (int f=0; f<mesh.nb_function_spaces(); ++f)
  {
    if (mesh.function_space(f).metadata<int>("type") == Entity::ELEMS)
    {
      FunctionSpace& elements = mesh.function_space(f);
      std::vector<int> extents = mesh.function_space(f).extents();
      int new_elem = elements.extents()[0];
      extents[0] = nb_elems[f];
      elements.resize(extents);
      ArrayView<int,1> elem_glb_idx       ( elements.field("glb_idx"       ) );
      ArrayView<int,1> elem_proc          ( elements.field("proc"          ) );
      ArrayView<int,1> elem_master_glb_idx( elements.field("master_glb_idx") );
      elem_nodes[f] = IndexView<int,2>( elements.field("nodes") );

      // Add elements at west boundary for each element at east boundary
      int cnt=0;
      for (int jelem=0; jelem<nb_east_bdry_elements; ++jelem)
      {
        if (east_bdry_elements[jelem].f == f)
        {
          int east_elem=east_bdry_elements[jelem].e;
          elem_proc(new_elem) = elem_proc(east_elem);
          elem_glb_idx(new_elem) = ++elements_max_glb_idx;
          elem_master_glb_idx(new_elem) = elem_glb_idx(east_elem);
          for (int n1=0; n1<nb_nodes_per_elem[f]; ++n1)
          {
            int east_node = elem_nodes[f](east_elem,n1);

            // if east-node is on boundary element but not ON east boundary
            if ( std::abs( coords(east_node,XX)-max[XX] ) > tol )
            {
              // Create new ghost nodes at west side with east-node as master
              elem_nodes[f](new_elem,n1) = -glb_idx(east_node);
              new_nodes.insert(-glb_idx(east_node));
            }
            else // if east-node is ON the east boundary
            {
              // Find matching node on west bdry element, and make west the master of east
              for (int ielem=0; ielem<nb_west_bdry_elements; ++ielem)
              {
                int func_space_idx = west_bdry_elements[ielem].f;
                int west_elem      = west_bdry_elements[ielem].e;
                for (int n2=0; n2<nb_nodes_per_elem[func_space_idx]; ++n2)
                {
                  int west_node = elem_nodes[func_space_idx](west_elem,n2);
                  if ( std::abs(coords(west_node,YY)-coords(east_node,YY))<tol && std::abs(coords(west_node,XX)-min[XX])<tol )
                  {
                    // Found matching node --> this east node has to become a ghost node
                    transform_to_ghost.insert(east_node);
                    master_glb_idx(east_node) = glb_idx(west_node);
                    elem_nodes[f](new_elem,n1) = west_node;
                  }
                }
              }
            }
          }
          ++new_elem;
        }
      }
    }
  }

  //std::cout << "going to move ghost nodes to the back" << std::endl;
  // Move the transformed ghost nodes to the back

  if(0) {
    for (int f=0; f<mesh.nb_function_spaces(); ++f)
    {
      if (mesh.function_space(f).metadata<int>("type") == Entity::ELEMS)
      {
        FunctionSpace& elements = mesh.function_space(f);
        IndexView<int,2> elem_nodes( elements.field("nodes") );
        int nb_elems = elem_nodes.extents()[0];
        int nb_nodes_per_elem = elem_nodes.extents()[1];
        for (int elem=0; elem<nb_elems; ++elem)
        {
          for (int n=0; n<nb_nodes_per_elem; ++n)
          {
            int node = elem_nodes(elem,n);
            if (elem_nodes(elem,n) > 0 )
              elem_nodes(elem,n) = glb_idx(node);
          }
        }
      }
    }
  }


  std::vector< int > orig_glb_idx(glb_idx.size());
  std::copy( glb_idx.data(), glb_idx.data()+glb_idx.size(), orig_glb_idx.begin() );
  
  std::vector< int > orig_master_glb_idx( glb_idx.size() );
  std::copy( master_glb_idx.data(), master_glb_idx.data()+master_glb_idx.size(), orig_master_glb_idx.begin() );

  std::map<int,int> node_orig_glb_to_loc;
  for (int node=0; node<nb_nodes; ++node)
    node_orig_glb_to_loc[ orig_glb_idx[node] ] = node;

  int nodes_max_glb_idx = nodes.metadata<int>("max_glb_idx");

  int decrease=0;
  int ghost_glb_idx=nodes_max_glb_idx-transform_to_ghost.size();
  for (int node=0; node<nb_nodes; ++node)
  {
    glb_idx( node ) = orig_glb_idx[ node ] + decrease;
    if( transform_to_ghost.count(node) )
    {
      ++ghost_glb_idx;
      if( glb_idx(node) != ghost_glb_idx)
      {
        glb_idx( node ) = ghost_glb_idx;
        --decrease;
      }
    }
  }

  for (int node=0; node<nb_nodes; ++node)
  {
    master_glb_idx( node ) = glb_idx( node_orig_glb_to_loc[ orig_master_glb_idx[ node ] ] );
  }


  // Now add new nodes
  int nb_ghost_nodes = new_nodes.size();
  nodes.metadata().set("nb_ghost_nodes",nb_ghost_nodes);
  std::vector<int> nodes_extents = nodes.extents();
  nodes_extents[0] = nb_nodes+nb_ghost_nodes;
  nodes.resize( nodes_extents ); // WARNING! ArrayViews no longer valid
  coords         = ArrayView<double,2>( nodes.field("coordinates"   ) );
  glb_idx        = ArrayView<int,   1>( nodes.field("glb_idx"       ) );
  master_glb_idx = ArrayView<int,   1>( nodes.field("master_glb_idx") );
  proc           = ArrayView<int,   1>( nodes.field("proc"          ) );
  
  std::map<int,int> node_mapping;
  int new_node = nb_nodes;
  for( std::set<int>::iterator it = new_nodes.begin(); it!=new_nodes.end(); ++it)
  {
    int orig_gid = -(*it);
    int master_node = node_orig_glb_to_loc[ orig_gid ];
    //std::cout << "new node " << new_node << std::endl;
    //std::cout << "orig_gid = " << orig_gid << std::endl;
    //std::cout << "master_node = " << master_node << std::endl;
    coords(new_node,XX) = coords(master_node,XX) - dx;
    coords(new_node,YY) = coords(master_node,YY);
    glb_idx(new_node) = ++nodes_max_glb_idx;
    master_glb_idx(new_node) = glb_idx(master_node);
    proc(new_node) = proc(master_node);
    node_mapping[ -orig_gid ] = new_node;
    ++new_node;
  }
  nodes.metadata().set("max_glb_idx",nodes_max_glb_idx);

  // Now fix element connectivity
  for (int f=0; f<mesh.nb_function_spaces(); ++f)
  {
    if (mesh.function_space(f).metadata<int>("type") == Entity::ELEMS)
    {
      FunctionSpace& elements = mesh.function_space(f);
      for (int elem=0; elem<elements.extents()[0]; ++elem)
      {
        for (int n=0; n<nb_nodes_per_elem[f]; ++n)
        {
          int gid = elem_nodes[f](elem,n);
          if ( gid < 0 )
          {
            elem_nodes[f](elem,n) = node_mapping[ gid ];
          }
        }
      }
    }
  }

  for (int f=0; f<mesh.nb_function_spaces(); ++f)
  {
    if (mesh.function_space(f).metadata<int>("type") == Entity::ELEMS)
      mesh.function_space(f).metadata().set("max_glb_idx",elements_max_glb_idx);
    if (mesh.function_space(f).metadata<int>("type") == Entity::FACES)
      mesh.function_space(f).metadata().set("max_glb_idx",elements_max_glb_idx);
  }


#if 0
  FunctionSpace& edges   = mesh.function_space( "edges" );
  FieldT<int>& edge_to_elem   = edges.field<int>( "to_elem" );
  FieldT<int>& edge_nodes   = edges.field<int>( "nodes" );
  int nb_edges = edges.extents()[0];

  double min[2];
  double max[2];

  for (int node=0; node<nb_nodes; ++node)
  {
    min[XX] = std::min( min[XX], coords(XX,node) );
    min[YY] = std::min( min[YY], coords(YY,node) );
    max[XX] = std::max( max[XX], coords(XX,node) );
    max[YY] = std::max( max[YY], coords(YY,node) );
  }
  double dx = max[XX]-min[XX];

  std::vector<int> bdry_edges; bdry_edges.reserve(nb_nodes);

  // assemble boundary edges
  for (int edge=0; edge<nb_edges; ++edge)
  {
    if ( edge_to_elem(3,edge) < 0 )
    {
      bdry_edges.push_back(edge);
    }
  }
  int nb_bdry_edges = bdry_edges.size();

  std::vector<int> west_bdry_edges;  west_bdry_edges.reserve(nb_bdry_edges/2);
  std::vector<int> east_bdry_edges;  east_bdry_edges.reserve(nb_bdry_edges/2);
  std::vector<int> north_bdry_edges; north_bdry_edges.reserve(nb_bdry_edges/4);
  std::vector<int> south_bdry_edges; south_bdry_edges.reserve(nb_bdry_edges/4);

  std::cout << "dx = " << dx << std::endl;
  std::cout << "nb_bdry_edges = " << nb_bdry_edges << std::endl;

  double tol( 1.e-6 );

  for (int jedge=0; jedge<nb_bdry_edges; ++jedge)
  {
    int edge = bdry_edges[jedge];
    if      (   std::abs( coords(XX,edge_nodes(0,edge)) - min[XX] ) < tol
             && std::abs( coords(XX,edge_nodes(1,edge)) - min[XX] ) < tol )
      west_bdry_edges.push_back(edge);
    else if (   std::abs( coords(XX,edge_nodes(0,edge)) - max[XX] ) < tol
             && std::abs( coords(XX,edge_nodes(1,edge)) - max[XX] ) < tol )
      east_bdry_edges.push_back(edge);
    else if (   std::abs( coords(YY,edge_nodes(0,edge)) - max[YY] ) < tol
             && std::abs( coords(YY,edge_nodes(1,edge)) - max[YY] ) < tol )
      north_bdry_edges.push_back(edge);
    else if (   std::abs( coords(YY,edge_nodes(0,edge)) - min[YY] ) < tol
             && std::abs( coords(YY,edge_nodes(1,edge)) - min[YY] ) < tol )
      south_bdry_edges.push_back(edge);
    else
      throw std::runtime_error("Unidentified boundary edge");
  }


  std::vector<int> periodic_edges(nb_edges);
  for (int jedge=0; jedge<nb_edges; ++jedge)
    periodic_edges[jedge] = jedge;

  int nb_east_bdry_edges = east_bdry_edges.size();
  int nb_west_bdry_edges = west_bdry_edges.size();
  for (int jedge=0; jedge<nb_east_bdry_edges; ++jedge)
  {
    int east_edge = east_bdry_edges[jedge];
    for (int jnode=0; jnode<2; ++jnode)
    {
      int east_node = edge_nodes(jnode,east_edge);
      if( master_glb_idx(0,east_node) == glb_idx(0,east_node) )
      {
        double y0 = coords(YY,east_node);
        for (int iedge=0; iedge<nb_west_bdry_edges; ++iedge)
        {
          int west_edge = west_bdry_edges[iedge];
          for (int inode=0; inode<2; ++inode)
          {
            int west_node = edge_nodes(inode,west_edge);
            double y1 = coords(YY,west_node);
            if ( std::abs(y0-y1) < tol )
            {
              master_glb_idx(0,east_node) = glb_idx(0,west_node);
              proc(0,east_node) = proc(0,west_node);
              //periodic_edges[east_edge] = west_edge;
            }
          }
        }
      }
      if( master_glb_idx(0,east_node) == glb_idx(0,east_node) )
      {
        throw std::runtime_error("Could not find a periodic node! Should not happen");
      }
    }
  }





#endif

#if 0

  // Recreate all its edges at the west side
  // Create ghost-edges/nodes to make seamless connectivity

  std::vector<int> edge_bounds( edges.bounds() );
  edge_bounds[1] = nb_edges + nb_east_bdry_edges;
  edges.resize(edge_bounds);

  for (int jedge=0; jedge<nb_east_bdry_edges; ++jedge)
  {
    int east_edge = east_bdry_edges[jedge];
    int west_edge = periodic_edges[east_edge];


    // Find the element of this edge
    int func_space_idx = edge_to_elem(0,east_edge);
    int elem           = edge_to_elem(1,east_edge);

    FunctionSpace& func_space = mesh.function_space(func_space_idx);

    for (int jnode=0; jnode<2; ++jnode)
    {
      int east_node = edge_nodes(jnode,east_edge);
      std::cout << "moving edge for node " << east_node << std::endl;
      std::cout << "  connected to element " << func_space.name() << " ["<<elem<<"]"<<std::endl;

      FieldT<int>& elem_edges = func_space.field<int>("to_edge");
      for (int j=0; j<elem_edges.bounds()[0]; ++j)
      {
        int elem_edge = elem_edges(j,elem);
        if( elem_edge != east_edge )
        {
          if ( edge_nodes(0,elem_edge) == east_node )
          {
            edge_nodes(0,elem_edge) = master_glb_idx(0,east_node);
          }
          else if ( edge_nodes(1,elem_edge) == east_node )
          {
            edge_nodes(1,elem_edge) = master_glb_idx(1,east_node);
          }
        }
      }
    }
  }

  std::cout << "periodic_edges" << std::endl;
  for (int edge=0; edge<nb_edges; ++edge)
  {
    if (periodic_edges[edge] >= 0 )
    {
      std::cout << edge << " : " << periodic_edges[edge] << std::endl;
    }
  }
#endif
}


// ------------------------------------------------------------------
// C wrapper interfaces to C++ routines

void atlas__build_periodic_boundaries ( Mesh* mesh) {
  build_periodic_boundaries(*mesh);
}

// ------------------------------------------------------------------


} // namespace atlas

