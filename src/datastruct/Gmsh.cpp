#include "Gmsh.hpp"
#include "Mesh.hpp"
#include "FunctionSpace.hpp"
#include "Field.hpp"
#include <iostream>
#include <fstream>
#include <stdexcept>
namespace ecmwf {

Gmsh::~Gmsh()
{ 
}

Mesh& Gmsh::read(const std::string& file_path)
{
  Mesh* mesh = new Mesh();
  std::ifstream file;
  file.open( file_path.c_str() , std::ios::in );
  std::string line;
  for (int i=0; i<4; ++i)
    std::getline(file,line);

  // Create nodes
  int nb_nodes;
  file >> nb_nodes;
  std::vector<int> bounds(2);
  bounds[0] = Field::NB_VARS;
  bounds[1] = nb_nodes;
  FunctionSpace& nodes_2d    = mesh->add_function_space( new FunctionSpace( "nodes_2d", "Lagrange_P1", bounds ) );
  FieldT<double>& coords = nodes_2d.create_field<double>("coords",2);
  FieldT<int>& glb_idx = nodes_2d.create_field<int>("glb_idx",1);
  std::map<int,int> glb_to_loc;
  int g;
  double z;
  for (int n=0; n<nb_nodes; ++n)
  {
    file >> g >> coords(XX,n) >> coords(YY,n) >> z;
    glb_idx(0,n) = g;
    glb_to_loc[ g ] = n;
  }
  for (int i=0; i<3; ++i)
    std::getline(file,line);

  // Find out which element types are inside
  int nb_elements;
  file >> nb_elements;
  int position = file.tellg();
  std::vector<int> nb_etype(20,0);
  int etype;
  for (int e=0; e<nb_elements; ++e)
  {
    file >> g >> etype;
    ++nb_etype[etype];
    std::getline(file,line); // finish line
  }

  // Allocate data structures for quads, triags, edges

  int nb_quads = nb_etype[QUAD];
  bounds[1] = nb_quads;
  FunctionSpace& quads      = mesh->add_function_space( new FunctionSpace( "quads", "Lagrange_P1", bounds ) );
  FieldT<int>& quad_nodes   = quads.create_field<int>("nodes",4);
  FieldT<int>& quad_glb_idx = quads.create_field<int>("glb_idx",1);

  int nb_triags = nb_etype[TRIAG];
  bounds[1] = nb_triags;
  FunctionSpace& triags      = mesh->add_function_space( new FunctionSpace( "triags", "Lagrange_P1", bounds ) );
  FieldT<int>& triag_nodes   = triags.create_field<int>("nodes",3);
  FieldT<int>& triag_glb_idx = triags.create_field<int>("glb_idx",1);

  int nb_edges = nb_etype[LINE];
  bounds[1] = nb_edges;
  FunctionSpace& edges      = mesh->add_function_space( new FunctionSpace( "edges", "Lagrange_P1", bounds ) );
  FieldT<int>& edge_nodes   = edges.create_field<int>("nodes",2);
  FieldT<int>& edge_glb_idx = edges.create_field<int>("glb_idx",1);


  // Now read all elements
  file.seekg(position,std::ios::beg);
  int dummy, gn0, gn1, gn2, gn3;
  int quad=0, triag=0, edge=0;
  for (int e=0; e<nb_elements; ++e)
  {
    file >> g >> etype >> dummy >> dummy >> dummy;
    switch( etype )
    {
      case(QUAD):
        file >> gn0 >> gn1 >> gn2 >> gn3;
        quad_glb_idx(0,quad) = g;
        quad_nodes(0,quad) = glb_to_loc[gn0];
        quad_nodes(1,quad) = glb_to_loc[gn1];
        quad_nodes(2,quad) = glb_to_loc[gn2];
        quad_nodes(3,quad) = glb_to_loc[gn3];
        ++quad;
        break;
      case(TRIAG):
        file >> gn0 >> gn1 >> gn2;
        triag_glb_idx(0,triag) = g;
        triag_nodes(0,triag) = glb_to_loc[gn0];
        triag_nodes(1,triag) = glb_to_loc[gn1];
        triag_nodes(2,triag) = glb_to_loc[gn2];
        ++triag;
        break;
      case(LINE):
        file >> gn0 >> gn1;
        edge_glb_idx(0,edge) = g;
        edge_nodes(0,edge) = glb_to_loc[gn0];
        edge_nodes(1,edge) = glb_to_loc[gn1];
        ++edge;
        break;
      default:
        throw std::runtime_error("ERROR: element type not supported");
    }
  }
  file.close();
  return *mesh;
}

void Gmsh::write(Mesh& mesh, const std::string& file_path)
{
  FunctionSpace& nodes_2d   = mesh.function_space( "nodes_2d" );
  FieldT<double>& coords    = nodes_2d.field<double>( "coords" );
  FieldT<int>& glb_idx      = nodes_2d.field<int>( "glb_idx" );
  int nb_nodes = nodes_2d.bounds()[1];

  FunctionSpace& quads       = mesh.function_space( "quads" );
  FieldT<int>& quad_nodes    = quads.field<int>( "nodes" );
  FieldT<int>& quad_glb_idx  = quads.field<int>( "glb_idx" );
  int nb_quads = quads.bounds()[1];

  FunctionSpace& triags      = mesh.function_space( "triags" );
  FieldT<int>& triag_nodes   = triags.field<int>( "nodes" );
  FieldT<int>& triag_glb_idx = triags.field<int>( "glb_idx" );
  int nb_triags = triags.bounds()[1];

  FunctionSpace& edges       = mesh.function_space( "edges" );
  FieldT<int>& edge_nodes    = edges.field<int>( "nodes" );
  FieldT<int>& edge_glb_idx  = edges.field<int>( "glb_idx" );
  int nb_edges = edges.bounds()[1];


  std::ofstream file;
  file.open( file_path.c_str(), std::ios::out );

  file << "$MeshFormat\n";
  file << "2.2 0 8\n";
  file << "$EndMeshFormat\n";
  file << "$Nodes\n";
  file << nb_nodes << "\n";
  for( int n=0; n<nb_nodes; ++n)
    file << glb_idx(0,n) << " " << coords(XX,n) << " " << coords(YY,n) << " " << 0. << "\n";
  file << "$EndNodes\n";
  file << "$Elements\n";
  file << nb_quads+nb_triags+nb_edges << "\n";
  for( int e=0; e<nb_quads; ++e)
  {
    file << quad_glb_idx(0,e) << " 3 2 1 1";
    for( int n=0; n<4; ++n )
      file << " " << glb_idx(0,quad_nodes(n,e));
    file << "\n";
  }
  for( int e=0; e<nb_triags; ++e)
  {
    file << triag_glb_idx(0,e) << " 2 2 1 1";
    for( int n=0; n<3; ++n )
      file << " " << glb_idx(0,triag_nodes(n,e));
    file << "\n";
  }
  for( int e=0; e<nb_edges; ++e)
  {
    file << edge_glb_idx(0,e) << " 1 2 2 1";
    for( int n=0; n<2; ++n )
      file << " " << glb_idx(0,edge_nodes(n,e));
    file << "\n";
  }
  file << "$EndElements\n";
  file << std::flush;
}

// ------------------------------------------------------------------
// C wrapper interfaces to C++ routines

Gmsh* ecmwf__Gmsh__new () {
  return new Gmsh();
}

void ecmwf__Gmsh__delete (Gmsh* This) {
  delete This;
}

Mesh* ecmwf__Gmsh__read (Gmsh* This, char* file_path) {
  return &This->read( std::string(file_path) );
}

void ecmwf__Gmsh__write (Gmsh* This, Mesh* mesh, char* file_path) {
  This->write( *mesh, std::string(file_path) );
}

// ------------------------------------------------------------------


} // namespace ecmwf

