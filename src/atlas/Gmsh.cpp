#include <iostream>
#include <fstream>
#include <stdexcept>
#include <cmath>
#include "Gmsh.hpp"
#include "Mesh.hpp"
#include "FunctionSpace.hpp"
#include "Field.hpp"
#include "Parameters.hpp"
namespace atlas {

double pi = std::acos(-1.);
double scaling = pi/180.;

Gmsh::~Gmsh()
{ 
}

Mesh& Gmsh::read(const std::string& file_path)
{
  Mesh* mesh = new Mesh();
  std::ifstream file;
  file.open( file_path.c_str() , std::ios::in );
  if( !file.is_open() )
    throw std::runtime_error("Could not open file "+file_path);

  std::string line;
  for (int i=0; i<4; ++i)
    std::getline(file,line);

  // Create nodes
  int nb_nodes;
  file >> nb_nodes;
  std::cout << "nb_nodes = " << nb_nodes << std::endl;
  std::vector<int> bounds(2);
  bounds[0] = Field::NB_VARS;
  bounds[1] = nb_nodes;
  FunctionSpace& nodes_2d    = mesh->add_function_space( new FunctionSpace( "nodes_2d", "Lagrange_P0", bounds ) );
  nodes_2d.metadata().set("type",static_cast<int>(NODES));
  FieldT<double>& coords = nodes_2d.create_field<double>("coordinates",2);
  FieldT<int>& glb_idx = nodes_2d.create_field<int>("glb_idx",1);
  FieldT<int>& master_glb_idx = nodes_2d.create_field<int>("master_glb_idx",1);
  FieldT<int>& proc = nodes_2d.create_field<int>("proc",1);

  std::map<int,int> glb_to_loc;
  int g;
  double x,y,z;
  int max_glb_idx=0;
  for (int n=0; n<nb_nodes; ++n)
  {
    file >> g >> x >> y >> z;
    glb_idx(n) = g;
    coords(XX,n) = x*scaling;
    coords(YY,n) = y*scaling;
    glb_to_loc[ g ] = n;
    master_glb_idx(n) = g;
    proc(n) = 0;
    max_glb_idx = std::max(max_glb_idx, g);
  }
  for (int i=0; i<3; ++i)
    std::getline(file,line);
  nodes_2d.metadata().set("nb_owned",nb_nodes);
  nodes_2d.metadata().set("max_glb_idx",max_glb_idx);

  // Find out which element types are inside
  int nb_elements;
  file >> nb_elements;
  int position = file.tellg();
  std::vector<int> nb_etype(20,0);
  int elements_max_glb_idx(0);
  int etype;
  for (int e=0; e<nb_elements; ++e)
  {
    file >> g >> etype; std::getline(file,line); // finish line
    ++nb_etype[etype];
    elements_max_glb_idx = std::max(elements_max_glb_idx,g);
  }

  // Allocate data structures for quads, triags, edges

  int nb_quads = nb_etype[QUAD];
  bounds[1] = nb_quads;
  FunctionSpace& quads      = mesh->add_function_space( new FunctionSpace( "quads", "Lagrange_P1", bounds ) );
  quads.metadata().set("type",static_cast<int>(ELEMS));
  FieldT<int>& quad_nodes   = quads.create_field<int>("nodes",4);
  FieldT<int>& quad_glb_idx = quads.create_field<int>("glb_idx",1);
  FieldT<int>& quad_master_glb_idx = quads.create_field<int>("master_glb_idx",1);
  FieldT<int>& quad_proc = quads.create_field<int>("proc",1);

  int nb_triags = nb_etype[TRIAG];
  bounds[1] = nb_triags;
  FunctionSpace& triags      = mesh->add_function_space( new FunctionSpace( "triags", "Lagrange_P1", bounds ) );
  triags.metadata().set("type",static_cast<int>(ELEMS));
  FieldT<int>& triag_nodes   = triags.create_field<int>("nodes",3);
  FieldT<int>& triag_glb_idx = triags.create_field<int>("glb_idx",1);
  FieldT<int>& triag_master_glb_idx = triags.create_field<int>("master_glb_idx",1);
  FieldT<int>& triag_proc = triags.create_field<int>("proc",1);

  int nb_edges = nb_etype[LINE];
  nb_edges = 0;
  bounds[1] = nb_edges;
  FunctionSpace& edges      = mesh->add_function_space( new FunctionSpace( "edges", "Lagrange_P1", bounds ) );
  edges.metadata().set("type",static_cast<int>(FACES));
  FieldT<int>& edge_nodes   = edges.create_field<int>("nodes",2);
  FieldT<int>& edge_glb_idx = edges.create_field<int>("glb_idx",1);
  FieldT<int>& edge_master_glb_idx = edges.create_field<int>("master_glb_idx",1);
  FieldT<int>& edge_proc = edges.create_field<int>("proc",1);



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
        quad_glb_idx(quad) = g;
        quad_master_glb_idx(quad) = g;
        quad_proc(quad) = 0;
        quad_nodes(0,quad) = F_IDX( glb_to_loc[gn0] );
        quad_nodes(1,quad) = F_IDX( glb_to_loc[gn1] );
        quad_nodes(2,quad) = F_IDX( glb_to_loc[gn2] );
        quad_nodes(3,quad) = F_IDX( glb_to_loc[gn3] );
        ++quad;
        break;
      case(TRIAG):
        file >> gn0 >> gn1 >> gn2;
        triag_glb_idx(triag) = g;
        triag_master_glb_idx(triag) = g;
        triag_proc(triag) = 0;
        triag_nodes(0,triag) = F_IDX( glb_to_loc[gn0] );
        triag_nodes(1,triag) = F_IDX( glb_to_loc[gn1] );
        triag_nodes(2,triag) = F_IDX( glb_to_loc[gn2] );
        ++triag;
        break;
      case(LINE):
        file >> gn0 >> gn1;
        //edge_glb_idx(edge) = g;
        //edge_master_glb_idx(edge) = g;
        //edge_nodes(0,edge) = F_IDX( glb_to_loc[gn0] );
        //edge_nodes(1,edge) = F_IDX( glb_to_loc[gn1] );
        //edge_proc(proc) = 0;
        //++edge;
        break;
      case(POINT):
        file >> gn0;
        break;
      default:
        std::cout << "etype " << etype << std::endl;
        throw std::runtime_error("ERROR: element type not supported");
    }
  }
  quads.metadata().set("nb_owned",nb_etype[QUAD]);
  triags.metadata().set("nb_owned",nb_etype[TRIAG]);
  edges.metadata().set("nb_owned",nb_etype[LINE]);
  quads.metadata().set("max_glb_idx",elements_max_glb_idx);
  triags.metadata().set("max_glb_idx",elements_max_glb_idx);
  edges.metadata().set("max_glb_idx",elements_max_glb_idx);
  file.close();
  return *mesh;
}

void Gmsh::write(Mesh& mesh, const std::string& file_path)
{
  bool spherical=true;
  bool include_ghost_elements = false;

  FunctionSpace& nodes_2d   = mesh.function_space( "nodes_2d" );
  FieldT<double>& coords    = nodes_2d.field<double>( "coordinates" );
  FieldT<int>& glb_idx      = nodes_2d.field<int>( "glb_idx" );
  int nb_nodes = nodes_2d.metadata<int>("nb_owned");
  nb_nodes = nodes_2d.bounds()[1];

  FunctionSpace& quads       = mesh.function_space( "quads" );
  FieldT<int>& quad_nodes    = quads.field<int>( "nodes" );
  FieldT<int>& quad_glb_idx  = quads.field<int>( "glb_idx" );
  int nb_quads = quads.metadata<int>("nb_owned");

  FunctionSpace& triags      = mesh.function_space( "triags" );
  FieldT<int>& triag_nodes   = triags.field<int>( "nodes" );
  FieldT<int>& triag_glb_idx = triags.field<int>( "glb_idx" );
  int nb_triags = triags.metadata<int>("nb_owned");

  FunctionSpace& edges       = mesh.function_space( "edges" );
  FieldT<int>& edge_nodes    = edges.field<int>( "nodes" );
  FieldT<int>& edge_glb_idx  = edges.field<int>( "glb_idx" );
  int nb_edges = edges.metadata<int>("nb_owned");
  nb_edges = edges.bounds()[1];

  if( include_ghost_elements == true )
  {
    nb_quads = quads.bounds()[1];
    nb_triags = triags.bounds()[1];
  }

  std::string ext = "";
  if( spherical ) ext = ".sphere";
  std::cout << "writing file " << file_path+ext << std::endl;
  std::ofstream file;
  file.open( (file_path+ext).c_str(), std::ios::out );

  file << "$MeshFormat\n";
  file << "2.2 0 8\n";
  file << "$EndMeshFormat\n";
  file << "$Nodes\n";
  file << nb_nodes << "\n";
  for( int n=0; n<nb_nodes; ++n)
  {
    double r     = 1.;
    double lon   = coords(XX,n);
    double lat   = coords(YY,n);

    if(spherical)
    {
      double x = r*std::cos(lat)*std::cos(lon);
      double y = r*std::cos(lat)*std::sin(lon);
      double z = r*std::sin(lat);
      file << glb_idx(n) << " " << x << " " << y << " " << z << "\n";
    }
    else
    {
      file << glb_idx(n) << " " << lon/scaling << " " << lat/scaling << " " << 0. << "\n";
    }
  }
  file << "$EndNodes\n";
  file << "$Elements\n";
  file << nb_quads+nb_triags+nb_edges << "\n";
  for( int e=0; e<nb_quads; ++e)
  {
    file << quad_glb_idx(e) << " 3 2 1 1";
    for( int n=0; n<4; ++n )
      file << " " << glb_idx( C_IDX( quad_nodes(n,e) ) );
    file << "\n";
  }
  for( int e=0; e<nb_triags; ++e)
  {
    file << triag_glb_idx(e) << " 2 2 1 1";
    for( int n=0; n<3; ++n )
      file << " " << glb_idx( C_IDX( triag_nodes(n,e) ) );
    file << "\n";
  }
  for( int e=0; e<nb_edges; ++e)
  {
   file << edge_glb_idx(e) << " 1 2 2 1";
    for( int n=0; n<2; ++n )
      file << " " << glb_idx( C_IDX( edge_nodes(n,e) ) );
    file << "\n";
  }
  file << "$EndElements\n";
  file << std::flush;
  file.close();


  if (nodes_2d.has_field("dual_volumes"))
  {
    FieldT<double>& field = nodes_2d.field<double>("dual_volumes");
    file.open( (field.name()+".msh").c_str() , std::ios::out );
    file << "$MeshFormat\n";
    file << "2.2 0 8\n";
    file << "$EndMeshFormat\n";
    file << "$NodeData\n";
    file << "1\n";
    file << "\""+field.name()+"\"\n";
    file << "1\n";
    file << "0.\n";
    file << "3\n";
    file << "0\n";
    file << "1\n";
    file << nb_nodes << "\n";
    for (int n=0; n<nb_nodes; ++n)
      file << glb_idx(n) << " " << field(n)<<"\n";
    file << "$EndNodeData\n";
    file << std::flush;
    file.close();
  }

  if (edges.has_field("dual_normals"))
  {
    FieldT<double>& field = edges.field<double>("dual_normals");
    file.open( (field.name()+".msh").c_str() , std::ios::out );
    file << "$MeshFormat\n";
    file << "2.2 0 8\n";
    file << "$EndMeshFormat\n";
    file << "$ElementNodeData\n";
    file << "1\n";
    file << "\""+field.name()+"\"\n";
    file << "1\n";
    file << "0.\n";
    file << "3\n";
    file << "0\n";
    file << "3\n";
    file << nb_edges << "\n";
    for (int edge=0; edge<nb_edges; ++edge)
    {
      file << edge_glb_idx(edge) << " 2";
      for (int n=0; n<2; ++n)
        file << " " << field(XX,edge) << " " << field(YY,edge) << " 0";
      file <<"\n";
    }
    file << "$EndElementNodeData\n";
    file << std::flush;
    file.close();
  }
}

// ------------------------------------------------------------------
// C wrapper interfaces to C++ routines

Gmsh* atlas__Gmsh__new () {
  return new Gmsh();
}

void atlas__Gmsh__delete (Gmsh* This) {
  delete This;
}

Mesh* atlas__Gmsh__read (Gmsh* This, char* file_path) {
  return &This->read( std::string(file_path) );
}

void atlas__Gmsh__write (Gmsh* This, Mesh* mesh, char* file_path) {
  This->write( *mesh, std::string(file_path) );
}

Mesh* atlas__read_gmsh (char* file_path)
{
  return &Gmsh::read(std::string(file_path));
}

void atlas__write_gmsh (Mesh* mesh, char* file_path) {
  Gmsh::write( *mesh, std::string(file_path) );
}


// ------------------------------------------------------------------


} // namespace atlas

