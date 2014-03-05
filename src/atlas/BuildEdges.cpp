#include <iostream>
#include <stdexcept>
#include <cmath>
#include <set>
#include "Mesh.hpp"
#include "FunctionSpace.hpp"
#include "Field.hpp"
#include "BuildEdges.hpp"
#include "Parameters.hpp"

namespace atlas {


void scan_function_space(
    FunctionSpace& func_space,
    std::vector< std::vector<int> >& node_to_face,
    std::vector<int>& tmp_edges,
    std::vector< std::vector<int> >& connectivity_edge_to_elem,
    int& nb_faces,
    int& nb_inner_faces )
{
  FieldT<int>& elem_nodes   = func_space.field<int>( "nodes" );
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
        face_nodes[jnode] = C_IDX( elem_nodes(face_node_numbering[f][jnode],e) );

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
            connectivity_edge_to_elem[face][2] = func_space.index();
            connectivity_edge_to_elem[face][3] = e;
            ++nb_inner_faces;
            found_face = true;
            break;
          }
        }
      }

      if (found_face == false)
      {
        connectivity_edge_to_elem.push_back( std::vector<int>(4,-1) );
        // if 3rd or 4th element stays negative, it is a bdry face
        connectivity_edge_to_elem[nb_faces][0] = func_space.index();
        connectivity_edge_to_elem[nb_faces][1] = e;
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

void build_element_to_edge_connectivity( Mesh& mesh, FieldT<int>& edge_to_elem )
{
  std::vector< FieldT<int>* > elem_to_edge( mesh.nb_function_spaces() );
  std::vector< std::vector<int> > edge_cnt( mesh.nb_function_spaces() );

  for( int func_space_idx=0; func_space_idx<mesh.nb_function_spaces(); ++func_space_idx)
  {
    FunctionSpace& func_space = mesh.function_space(func_space_idx);
    if( func_space.metadata<int>("type") == ELEMS )
    {
      int nb_edges_per_elem;
      if (func_space.name() == "quads") nb_edges_per_elem = 4;
      if (func_space.name() == "triags") nb_edges_per_elem = 3;
      elem_to_edge[func_space_idx] = &func_space.create_field<int>("to_edge",nb_edges_per_elem);
      edge_cnt[func_space_idx].resize( func_space.bounds()[1], 0);
    }
  }

  int nb_edges = edge_to_elem.bounds()[1];
  for( int edge=0; edge<nb_edges; ++edge)
  {
    for( int j=0; j<2; ++j)
    {
      int func_space_idx = C_IDX( edge_to_elem(0+2*j,edge) );
      int elem           = C_IDX( edge_to_elem(1+2*j,edge) );
      if ( elem >= 0 )
      {
        (*elem_to_edge[func_space_idx])(edge_cnt[func_space_idx][elem]++, elem) = F_IDX(edge);
      }
    }
  }
}


void build_pole_edges( Mesh& mesh, std::vector<int>& pole_edge_nodes, int& nb_pole_edges )
{
  FunctionSpace& nodes_2d   = mesh.function_space( "nodes_2d" );
  FieldT<double>& coords    = nodes_2d.field<double>( "coordinates" );
  FieldT<int>& glb_idx    = nodes_2d.field<int>( "glb_idx" );
  int nb_nodes = nodes_2d.bounds()[1];
  FunctionSpace& edges   = mesh.function_space( "edges" );
  FieldT<int>& edge_nodes    = edges.field<int>( "nodes" );
  FieldT<int>& edge_glb_idx    = edges.field<int>( "glb_idx" );
  int nb_edges = edges.bounds()[1];

  double ymin = nodes_2d.metadata<double>("ymin");
  double ymax = nodes_2d.metadata<double>("ymax");
  double xmin = nodes_2d.metadata<double>("xmin");
  double xmax = nodes_2d.metadata<double>("xmax");
  double dx   = xmax-xmin;
  double tol = 1e-6;
  std::vector<int> north_pole_edges;
  std::vector<int> south_pole_edges;

  std::cout << "ymax = " << ymax << std::endl;
  std::cout << "ymin = " << ymin << std::endl;
  std::set<int> north_nodes;
  std::set<int> south_nodes;

  for (int node=0; node<nb_nodes; ++node)
  {
    //std::cout << "node " << node << "   " << std::abs(coords(YY,node)-ymax) << std::endl;
    if ( std::abs(coords(YY,node)-ymax)<tol )
    {
      std::cout << "add north node " << node << "  coord " << coords(XX,node) << std::endl;
      north_nodes.insert(node);
    }
    else if ( std::abs(coords(YY,node)-ymin)<tol )
    {
      std::cout << "add south node " << node << "  coord " << coords(XX,node) << std::endl;
      south_nodes.insert(node);
    }
  }
  double pi = std::acos(-1.);
  std::cout << "\nnorth: " << north_nodes.size() << std::endl;
  for( std::set<int>::iterator it=north_nodes.begin(); it!=north_nodes.end(); ++it)
  {
    int node = *it;
    double x1 = coords(XX,node);
    double x2 = coords(XX,node) + 0.5*dx;
    double dist = 2.*pi;
    if( x1>=xmin-tol && x1<=(xmax-xmin)*0.5-tol )
    {
      int recip_node;
      for( std::set<int>::iterator itr=north_nodes.begin(); itr!=north_nodes.end(); ++itr)
      {
        int other_node = *itr;
        if( std::abs(coords(XX,other_node)-x2)<dist )
        {
          dist = std::abs(coords(XX,other_node)-x2);
          recip_node = other_node;
        }
      }
      pole_edge_nodes.push_back(node);
      pole_edge_nodes.push_back(recip_node);
      ++nb_pole_edges;
      std::cout << "connect " << glb_idx(node) << " and " << glb_idx(recip_node) << std::endl;
    }
  }

  std::cout << "\nsouth: " << south_nodes.size() << std::endl;
  for( std::set<int>::iterator it=south_nodes.begin(); it!=south_nodes.end(); ++it)
  {
    int node = *it;
    double x1 = coords(XX,node);
    double x2 = coords(XX,node) + 0.5*dx;
    double dist = 2.*pi;
    if( x1>=xmin-tol && x1<=(xmax-xmin)*0.5-tol )

    {
      int recip_node;
      for( std::set<int>::iterator itr=south_nodes.begin(); itr!=south_nodes.end(); ++itr)
      {
        int other_node = *itr;
        if( std::abs(coords(XX,other_node)-x2)<dist )
        {
          dist = std::abs(coords(XX,other_node)-x2);
          recip_node = other_node;
        }
      }
      if (coords(XX,node) >= xmin)
      {
        pole_edge_nodes.push_back(node);
        pole_edge_nodes.push_back(recip_node);
        ++nb_pole_edges;
        std::cout << "connect " << glb_idx(node) << " and " << glb_idx(recip_node) << std::endl;
      }
    }
  }
}


void build_edges( Mesh& mesh )
{
  FunctionSpace& nodes_2d   = mesh.function_space( "nodes_2d" );
  FieldT<double>& coords    = nodes_2d.field<double>( "coordinates" );
  FieldT<int>& glb_idx      = nodes_2d.field<int>( "glb_idx" );
  FieldT<int>& master_glb_idx      = nodes_2d.field<int>( "master_glb_idx" );
  int nb_nodes = nodes_2d.bounds()[1];

  FunctionSpace& quads       = mesh.function_space( "quads" );
  FunctionSpace& triags      = mesh.function_space( "triags" );

  std::vector< std::vector<int> > node_to_face(nb_nodes);
  std::vector< int > tmp_edges; tmp_edges.reserve(4*nb_nodes);
  std::vector< std::vector<int> > to_elem;
  to_elem.reserve(8*nb_nodes);
  int nb_faces = 0;
  int nb_inner_faces = 0;


  scan_function_space(quads, node_to_face,tmp_edges,to_elem,nb_faces,nb_inner_faces);
  scan_function_space(triags,node_to_face,tmp_edges,to_elem,nb_faces,nb_inner_faces);


  int nb_edges = nb_faces;
  std::vector<int> bounds(2);
  bounds[0] = Field::NB_VARS;
  bounds[1] = nb_edges;
  FunctionSpace& edges       = mesh.function_space("edges");
  edges.resize(bounds);
  FieldT<int>& edge_nodes    = edges.field<int>( "nodes" );
  FieldT<int>& edge_glb_idx  = edges.field<int>( "glb_idx" );
  FieldT<int>& edge_master_glb_idx  = edges.field<int>( "master_glb_idx" );
  FieldT<int>& edge_proc  = edges.field<int>( "proc" );
  FieldT<int>& edge_to_elem  = edges.create_field<int>( "to_elem", 4 );

  std::map<int,std::vector<int> > node_glb_to_edge;
  int gid=edges.metadata<int>("max_glb_idx");
  int cnt=0;
  for(int edge=0; edge<nb_edges; ++edge)
  {
    edge_glb_idx(0,edge) = ++gid;
    edge_master_glb_idx(0,edge) = gid;
    edge_proc(0,edge) = 0;
    edge_nodes(0,edge)   = F_IDX( tmp_edges[cnt++] );
    edge_nodes(1,edge)   = F_IDX( tmp_edges[cnt++] );
    edge_to_elem(0,edge) = F_IDX( to_elem[edge][0] );
    edge_to_elem(1,edge) = F_IDX( to_elem[edge][1] );
    edge_to_elem(2,edge) = F_IDX( to_elem[edge][2] );
    edge_to_elem(3,edge) = F_IDX( to_elem[edge][3] );
    node_glb_to_edge[ glb_idx( C_IDX( edge_nodes(0,edge) ) ) ].push_back(edge);
    node_glb_to_edge[ glb_idx( C_IDX( edge_nodes(1,edge) ) ) ].push_back(edge);
  }

  // Element to edge connectivity
  build_element_to_edge_connectivity(mesh,edge_to_elem);


  // Element glb to loc
  std::vector< std::map<int,int> > elem_glb_to_loc(mesh.nb_function_spaces());
  for (int f=0; f<mesh.nb_function_spaces(); ++f)
  {
    if (mesh.function_space(f).metadata<int>("type") == ELEMS)
    {
      FieldT<int>& elem_glb_idx = mesh.function_space(f).field<int>("glb_idx");
      int nb_elems = elem_glb_idx.bounds()[1];
      for (int elem=0; elem<nb_elems; ++elem)
      {
        elem_glb_to_loc[f][ elem_glb_idx(elem) ] = elem;
      }
    }
  }

  // Mark certain edges as ghost and fix node ordering
  for (int edge=0; edge<nb_edges; ++edge)
  {
    int node0 = C_IDX( edge_nodes(0,edge) );
    int node1 = C_IDX( edge_nodes(1,edge) );
    int node0_master_gid = master_glb_idx(node0);
    int node1_master_gid = master_glb_idx(node1);
    bool node0_ghost = ( glb_idx(node0) != master_glb_idx(node0) );
    bool node1_ghost = ( glb_idx(node1) != master_glb_idx(node1) );
    if ( node0_ghost && node1_ghost )
    {
      //std::cout << edge_glb_idx(edge) << " is ghost" << std::endl;
      for (int jedge=0; jedge<node_glb_to_edge[node0_master_gid].size(); ++jedge)
      {
        int master_edge = node_glb_to_edge[node0_master_gid][jedge];
        int master_node0_master_gid = master_glb_idx( C_IDX( edge_nodes(0,master_edge) ) );
        int master_node1_master_gid = master_glb_idx( C_IDX( edge_nodes(1,master_edge) ) );
        if (node0_master_gid == master_node1_master_gid && node1_master_gid == master_node0_master_gid)
        {
          // wrong order
          //std::cout << "fixed wrong order for edge " << edge_glb_idx(edge) << std::endl;
          edge_nodes(0,edge) = F_IDX(node1);
          edge_nodes(1,edge) = F_IDX(node0);
          edge_master_glb_idx(edge) = edge_glb_idx(master_edge);
          break;
        }
        else if (node0_master_gid == master_node0_master_gid && node1_master_gid == master_node1_master_gid)
        {
          edge_master_glb_idx(edge) = edge_glb_idx(master_edge);
          break;
        }
      }
      //std::cout << "  master = " << edge_master_glb_idx(edge) << std::endl;
    }
    else if ( node0_ghost || node1_ghost )
    {
      ElementRef left_element( C_IDX(edge_to_elem(0,edge)), C_IDX(edge_to_elem(1,edge)));
      FunctionSpace& left_func_space = mesh.function_space(left_element.f);
      FieldT<int>& left_glb_idx = left_func_space.field<int>("glb_idx");
      FieldT<int>& left_master_glb_idx = left_func_space.field<int>("master_glb_idx");
      FieldT<int>& left_nodes = left_func_space.field<int>("nodes");
      FieldT<int>& left_to_edges = left_func_space.field<int>("to_edge");
      int nb_edges_per_elem = left_to_edges.bounds()[0];
      int nb_nodes_per_elem = left_nodes.bounds()[0];
      if ( left_glb_idx(left_element.e) != left_master_glb_idx(left_element.e) )
      {
        // this is a ghost element, and its edges must be also.
        //std::cout << edge_glb_idx(edge) << " is ghost" << std::endl;
        int master_elem_gid = left_master_glb_idx(left_element.e);
        int master_elem = elem_glb_to_loc[left_element.f][master_elem_gid];
        for (int jedge=0; jedge<nb_edges_per_elem; ++jedge)
        {
          int master_edge = left_to_edges(jedge,master_elem);
          if (master_glb_idx(node0) == master_glb_idx( C_IDX(edge_nodes(0,master_edge) ) ) )
          {
            if (master_glb_idx(node1) == master_glb_idx( C_IDX(edge_nodes(1,master_edge) ) ) )
            {
              edge_master_glb_idx(edge) = edge_glb_idx(master_edge);
              break;
            }
          }
          else if (master_glb_idx(node0) == master_glb_idx( C_IDX( edge_nodes(1,master_edge) ) ) )
          {
            if (master_glb_idx(node1) == master_glb_idx( C_IDX( edge_nodes(0,master_edge) ) ) )
            {
              edge_master_glb_idx(edge) = edge_glb_idx(master_edge);
              break;
            }
          }
        }
        //stdcout << "  master = " << edge_master_glb_idx(edge) << std::endl;
      }
    }
  }

  int nb_pole_edges(0);
  std::vector<int> pole_edge_nodes;
  build_pole_edges( mesh, pole_edge_nodes, nb_pole_edges );
  bounds[1] += nb_pole_edges;
  edges.resize(bounds);
  cnt=0;
  for(int edge=nb_edges; edge<nb_edges+nb_pole_edges; ++edge)
  {
    edge_glb_idx(0,edge) = ++gid;
    edge_master_glb_idx(0,edge) = gid;
    edge_proc(0,edge) = 0;
    edge_nodes(0,edge)   = F_IDX(pole_edge_nodes[cnt++]);
    edge_nodes(1,edge)   = F_IDX(pole_edge_nodes[cnt++]);
    edge_to_elem(0,edge) = -1;
    edge_to_elem(1,edge) = -1;
    edge_to_elem(2,edge) = -1;
    edge_to_elem(3,edge) = -1;
    node_glb_to_edge[ glb_idx( C_IDX(edge_nodes(0,edge)) ) ].push_back(edge);
    node_glb_to_edge[ glb_idx( C_IDX(edge_nodes(1,edge)) ) ].push_back(edge);
  }


}


// ------------------------------------------------------------------
// C wrapper interfaces to C++ routines

void atlas__build_edges ( Mesh* mesh) {
  build_edges(*mesh);
}

// ------------------------------------------------------------------


} // namespace atlas

