#include "Gmsh.hpp"
#include "Mesh.hpp"
#include "FunctionSpace.hpp"
#include "Field.hpp"
#include <iostream>
#include <fstream>
#include <stdexcept>
namespace ecmwf {


void scan_function_space(
    FunctionSpace& func_space,
    std::vector< std::vector<int> >& node_to_face,
    std::vector<int>& tmp_edges,
    int& nb_faces,
    int& nb_inner_faces )
{
  FieldT<int>& nodes = func_space.field<int>( "nodes" );
  int nb_elems = func_space.bounds()[1];
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
        face_nodes[jnode] = nodes(face_node_numbering[f][jnode],e);

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
            ++nb_inner_faces;
            found_face = true;
            break;
          }
        }
      }

      if (found_face == false)
      {
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

void build_edges( Mesh& mesh )
{
  FunctionSpace& nodes_2d   = mesh.function_space( "nodes_2d" );
  FieldT<double>& coords    = nodes_2d.field<double>( "coords" );
  FieldT<int>& glb_idx      = nodes_2d.field<int>( "glb_idx" );
  int nb_nodes = nodes_2d.bounds()[1];

  int gid=1;
  FunctionSpace& quads       = mesh.function_space( "quads" );
  gid += quads.bounds()[1];
  FunctionSpace& triags      = mesh.function_space( "triags" );
  gid += triags.bounds()[1];

  std::vector< std::vector<int> > node_to_face(nb_nodes);
  std::vector< int > tmp_edges; tmp_edges.reserve(4*nb_nodes);
  int nb_faces = 0;
  int nb_inner_faces = 0;

  scan_function_space(quads, node_to_face,tmp_edges,nb_faces,nb_inner_faces);
  scan_function_space(triags,node_to_face,tmp_edges,nb_faces,nb_inner_faces);

  int nb_edges = nb_faces;
  std::vector<int> bounds(2);
  bounds[0] = Field::NB_VARS;
  bounds[1] = nb_edges;
  FunctionSpace& edges       = mesh.add_function_space( new FunctionSpace("edges", "LagrangeP1", bounds) );
  FieldT<int>& edge_nodes    = edges.create_field<int>( "nodes", 2 );
  FieldT<int>& edge_glb_idx  = edges.create_field<int>( "glb_idx", 1 );

  int cnt=0;
  for(int edge=0; edge<nb_edges; ++edge)
  {
    edge_glb_idx(0,edge) = gid++;
    edge_nodes(0,edge) = tmp_edges[cnt++];
    edge_nodes(1,edge) = tmp_edges[cnt++];
  }
}


// ------------------------------------------------------------------
// C wrapper interfaces to C++ routines

void ecmwf__build_edges ( Mesh* mesh) {
  build_edges(*mesh);
}

// ------------------------------------------------------------------


} // namespace ecmwf

