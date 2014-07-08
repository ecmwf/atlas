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
#include <fstream>
#include <stdexcept>
#include <cmath>
#include <limits>

#include "atlas/io/Gmsh.hpp"
#include "atlas/mesh/Mesh.hpp"
#include "atlas/mesh/FunctionSpace.hpp"
#include "atlas/mesh/Field.hpp"
#include "atlas/util/ArrayView.hpp"
#include "atlas/util/IndexView.hpp"
#include "atlas/mesh/Parameters.hpp"

namespace atlas {

 double deg_to_rad = M_PI / 180.;
 double rad_to_deg = 180. * M_1_PI;

 Gmsh::~Gmsh()
{ 
}

Mesh* Gmsh::read(const std::string& file_path)
{
    Mesh* mesh = new Mesh();
    Gmsh::read(file_path,*mesh);
    return mesh;
}

void Gmsh::read(const std::string& file_path, Mesh& mesh )
{
  std::ifstream file;
  file.open( file_path.c_str() , std::ios::in );
  if( !file.is_open() )
    throw std::runtime_error("Could not open file "+file_path);

  std::string line;
  while(line != "$Nodes")
    std::getline(file,line);

  // Create nodes
  int nb_nodes;
  file >> nb_nodes;

  std::vector<int> extents(2);

  extents[0] = nb_nodes;
  extents[1] = Field::UNDEF_VARS;
  
  if( mesh.has_function_space("nodes") )
  {
      if( mesh.function_space("nodes").extents()[0] != nb_nodes )
          throw std::runtime_error("existing nodes function space has incompatible number of nodes");
  }
  else
  {
    mesh.add_function_space( new FunctionSpace( "nodes", "Lagrange_P0", extents ) )
            .metadata().set("type",static_cast<int>(Entity::NODES));
  }

  FunctionSpace& nodes = mesh.function_space("nodes");

  if( ! nodes.has_field("coordinates") )
      nodes.create_field<double>("coordinates",3);
  if( ! nodes.has_field("glb_idx") )
      nodes.create_field<int>("glb_idx",1);
  if( ! nodes.has_field("partition") )
      nodes.create_field<int>("partition",1);

  ArrayView<double,2> coords         ( nodes.field("coordinates")    );
  ArrayView<int,   1> glb_idx        ( nodes.field("glb_idx")        );
  ArrayView<int,   1> part           ( nodes.field("partition")      );

  std::map<int,int> glb_to_loc;
  int g;
  double x,y,z;
  double xmax = -std::numeric_limits<double>::max();
  double zmax = -std::numeric_limits<double>::max();
  int max_glb_idx=0;
  for( size_t n = 0; n < nb_nodes; ++n )
  {
    file >> g >> x >> y >> z;
    glb_idx(n) = g;
    coords(n,XX) = x;
    coords(n,YY) = y;
    coords(n,ZZ) = z;
    glb_to_loc[ g ] = n;
    part(n) = 0;
    max_glb_idx = std::max(max_glb_idx, g);
    xmax = std::max(x,xmax);
    zmax = std::max(z,zmax);
  }
  if( xmax > 4*M_PI && zmax == 0. )
  {
    for( size_t n = 0; n < nb_nodes; ++n )
    {
      coords(n,XX) *= deg_to_rad;
      coords(n,YY) *= deg_to_rad;
    }
  }
  for (int i=0; i<3; ++i)
    std::getline(file,line);
  nodes.metadata().set("nb_owned",nb_nodes);
  nodes.metadata().set("max_glb_idx",max_glb_idx);

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
  extents[0] = nb_quads;
  FunctionSpace& quads      = mesh.add_function_space( new FunctionSpace( "quads", "Lagrange_P1", extents ) );
  quads.metadata().set("type",static_cast<int>(Entity::ELEMS));
  IndexView<int,2> quad_nodes          ( quads.create_field<int>("nodes",         4) );
  ArrayView<int,1> quad_glb_idx        ( quads.create_field<int>("glb_idx",       1) );
  ArrayView<int,1> quad_part           ( quads.create_field<int>("partition",     1) );

  int nb_triags = nb_etype[TRIAG];
  extents[0] = nb_triags;
  FunctionSpace& triags      = mesh.add_function_space( new FunctionSpace( "triags", "Lagrange_P1", extents ) );
  triags.metadata().set("type",static_cast<int>(Entity::ELEMS));
  IndexView<int,2> triag_nodes          ( triags.create_field<int>("nodes",         3) );
  ArrayView<int,1> triag_glb_idx        ( triags.create_field<int>("glb_idx",       1) );
  ArrayView<int,1> triag_part           ( triags.create_field<int>("partition",     1) );

  int nb_edges = nb_etype[LINE];
  IndexView<int,2> edge_nodes;
  ArrayView<int,1> edge_glb_idx;
  ArrayView<int,1> edge_part;

  if( nb_edges > 0 )
  {
    extents[0] = nb_edges;
    FunctionSpace& edges      = mesh.add_function_space( new FunctionSpace( "edges", "Lagrange_P1", extents ) );
    edges.metadata().set("type",static_cast<int>(Entity::FACES));
    edge_nodes   = IndexView<int,2> ( edges.create_field<int>("nodes",         2) );
    edge_glb_idx = ArrayView<int,1> ( edges.create_field<int>("glb_idx",       1) );
    edge_part    = ArrayView<int,1> ( edges.create_field<int>("partition",     1) );
  }

  // Now read all elements
  file.seekg(position,std::ios::beg);
  int dummy, gn0, gn1, gn2, gn3;
  int quad=0, triag=0, edge=0;
  int ntags, tags[100];
  for (int e=0; e<nb_elements; ++e)
  {
    file >> g >> etype >> ntags;
    for( int t=0; t<ntags; ++t ) file >> tags[t];
    int part=0;
    if( ntags > 3 ) part = std::max( part, *std::max_element(tags+3,tags+ntags-1) ); // one positive, others negative
    switch( etype )
    {
      case(QUAD):
        file >> gn0 >> gn1 >> gn2 >> gn3;
        quad_glb_idx(quad) = g;
        quad_part(quad) = part;
        quad_nodes(quad,0) = glb_to_loc[gn0];
        quad_nodes(quad,1) = glb_to_loc[gn1];
        quad_nodes(quad,2) = glb_to_loc[gn2];
        quad_nodes(quad,3) = glb_to_loc[gn3];
        ++quad;
        break;
      case(TRIAG):
        file >> gn0 >> gn1 >> gn2;
        triag_glb_idx(triag) = g;
        triag_part(triag) = part;
        triag_nodes(triag,0) = glb_to_loc[gn0];
        triag_nodes(triag,1) = glb_to_loc[gn1];
        triag_nodes(triag,2) = glb_to_loc[gn2];
        ++triag;
        break;
      case(LINE):
        file >> gn0 >> gn1;
        edge_glb_idx(edge) = g;
        edge_part(edge) = part;
        edge_nodes(edge,0) = glb_to_loc[gn0];
        edge_nodes(edge,1) = glb_to_loc[gn1];
        ++edge;
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
  quads.metadata().set("max_glb_idx",elements_max_glb_idx);
  triags.metadata().set("max_glb_idx",elements_max_glb_idx);

  if( nb_edges > 0 )
  {
    mesh.function_space("edges").metadata().set("nb_owned",nb_etype[LINE]);
    mesh.function_space("edges").metadata().set("max_glb_idx",elements_max_glb_idx);
  }

  file.close();
}

void Gmsh::write(Mesh& mesh, const std::string& file_path)
{

  bool include_ghost_elements = true;

  FunctionSpace& nodes    = mesh.function_space( "nodes" );
  ArrayView<double,2> coords  ( nodes.field( "coordinates" ) );
  ArrayView<int,   1> glb_idx ( nodes.field( "glb_idx" ) );
  int nb_nodes = nodes.metadata<int>("nb_owned");
  nb_nodes = nodes.extents()[0];

  FunctionSpace& quads       = mesh.function_space( "quads" );
  IndexView<int,2> quad_nodes   ( quads.field( "nodes" ) );
  ArrayView<int,1> quad_glb_idx ( quads.field( "glb_idx" ) );
  ArrayView<int,1> quad_part    ( quads.field( "partition" ) );
  int nb_quads = quads.metadata<int>("nb_owned");

  FunctionSpace& triags      = mesh.function_space( "triags" );
  IndexView<int,2> triag_nodes   ( triags.field( "nodes" ) );
  ArrayView<int,1> triag_glb_idx ( triags.field( "glb_idx" ) );
  ArrayView<int,1> triag_part    ( triags.field( "partition" ) );
  int nb_triags = triags.metadata<int>("nb_owned");

  int nb_edges(0);
  if( mesh.has_function_space("edges") )
  {
    FunctionSpace& edges       = mesh.function_space( "edges" );
    nb_edges = edges.extents()[0];
  }

  if( include_ghost_elements == true )
  {
    nb_quads = quads.extents()[0];
    nb_triags = triags.extents()[0];
  }

  for( int spherical=0; spherical<2; ++spherical )
  {
//  bool spherical=true;

  std::string ext = "";
  if( spherical ) ext = ".sphere";
  if( MPL::rank() == 0 )
    std::cout << "writing file " << file_path+ext << std::endl;
  std::ofstream file;
  file.open( (file_path+ext).c_str(), std::ios::out );

  file << "$MeshFormat\n";
  file << "2.2 0 8\n";
  file << "$EndMeshFormat\n";
  file << "$Nodes\n";
  file << nb_nodes << "\n";
  for( size_t n = 0; n < nb_nodes; ++n )
  {
    double r     = 1.;
    double lon   = coords(n,XX);
    double lat   = coords(n,YY);

    if(spherical)
    {
      double x = r*std::cos(lat)*std::cos(lon);
      double y = r*std::cos(lat)*std::sin(lon);
      double z = r*std::sin(lat);
      file << glb_idx(n) << " " << x << " " << y << " " << z << "\n";
    }
    else
    {
      file << glb_idx(n) << " " << lon << " " << lat << " " << 0. << "\n";
    }
  }
  file << "$EndNodes\n";
  file << "$Elements\n";
  file << nb_quads+nb_triags+nb_edges << "\n";
  for( int e=0; e<nb_quads; ++e)
  {
    file << quad_glb_idx(e) << " 3 4 1 1 1 " << quad_part(e);
    for( int n=0; n<4; ++n )
      file << " " << glb_idx( quad_nodes(e,n) );
    file << "\n";
  }
  for( int e=0; e<nb_triags; ++e)
  {
    file << triag_glb_idx(e) << " 2 4 1 1 1 " << triag_part(e);
    for( int n=0; n<3; ++n )
      file << " " << glb_idx( triag_nodes(e,n) );
    file << "\n";
  }

  if( mesh.has_function_space("edges") )
  {
    FunctionSpace& edges       = mesh.function_space( "edges" );
    IndexView<int,2> edge_nodes   ( edges.field( "nodes" ) );
    ArrayView<int,1> edge_glb_idx ( edges.field( "glb_idx" ) );
    if( edges.has_field("partition") )
    {
      ArrayView<int,1> edge_part ( edges.field( "partition" ) );
      for( int e=0; e<nb_edges; ++e)
      {
        file << edge_glb_idx(e) << " 1 4 1 1 1 " << edge_part(e);
        for( int n=0; n<2; ++n )
          file << " " << glb_idx( edge_nodes(e,n) );
        file << "\n";
      }
    }
    else
    {
      for( int e=0; e<nb_edges; ++e)
      {
        file << edge_glb_idx(e) << " 1 2 1 1";
        for( int n=0; n<2; ++n )
          file << " " << glb_idx( edge_nodes(e,n) );
        file << "\n";
      }

    }
  }
  file << "$EndElements\n";
  file << std::flush;
  file.close();

  }
  if (nodes.has_field("dual_volumes"))
  {
    std::ofstream file;
    ArrayView<double,1> field( nodes.field("dual_volumes") );
    std::stringstream filename; filename << "dual_volumes_p"<<MPL::rank()<<".msh";
    file.open( filename.str().c_str() , std::ios::out );
    file << "$MeshFormat\n";
    file << "2.2 0 8\n";
    file << "$EndMeshFormat\n";
    file << "$NodeData\n";
    file << "1\n";
    file << "\"dual_volumes\"\n";
    file << "1\n";
    file << "0.\n";
    file << "3\n";
    file << "0\n";
    file << "1\n";
    file << nb_nodes << "\n";
    for( size_t n = 0; n < nb_nodes; ++n )
      file << glb_idx(n) << " " << field(n)<<"\n";
    file << "$EndNodeData\n";
    file << std::flush;
    file.close();
  }

  if( mesh.has_function_space("edges") )
  {
    FunctionSpace& edges       = mesh.function_space( "edges" );
    ArrayView<int,1> edge_glb_idx ( edges.field( "glb_idx" ) );
    nb_edges = edges.extents()[0];

    if (edges.has_field("dual_normals"))
    {
      std::ofstream file;
      ArrayView<double,2> field ( edges.field("dual_normals") );
      std::stringstream filename; filename << "dual_normals_p"<<MPL::rank()<<".msh";
      file.open( filename.str().c_str() , std::ios::out );
      file << "$MeshFormat\n";
      file << "2.2 0 8\n";
      file << "$EndMeshFormat\n";
      file << "$ElementNodeData\n";
      file << "1\n";
      file << "\"dual_normals\"\n";
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
          file << " " << field(edge,XX) << " " << field(edge,YY) << " 0";
        file <<"\n";
      }
      file << "$EndElementNodeData\n";
      file << std::flush;
      file.close();
    }

    if (edges.has_field("skewness"))
    {
      std::ofstream file;
      ArrayView<double,1> field( edges.field("skewness") );
      file.open( "skewness.msh" , std::ios::out );
      file << "$MeshFormat\n";
      file << "2.2 0 8\n";
      file << "$EndMeshFormat\n";
      file << "$ElementNodeData\n";
      file << "1\n";
      file << "\"skewness\"\n";
      file << "1\n";
      file << "0.\n";
      file << "3\n";
      file << "0\n";
      file << "1\n";
      file << nb_edges << "\n";
      for (int edge=0; edge<nb_edges; ++edge)
      {
        file << edge_glb_idx(edge) << " 2";
        for (int n=0; n<2; ++n)
          file << " " << field(edge);
        file <<"\n";
      }
      file << "$EndElementNodeData\n";
      file << std::flush;
      file.close();
    }
  }

}

void Gmsh::write3dsurf(Mesh &mesh, const std::string& file_path)
{
    size_t nb_nodes  = 0;
    size_t nb_triags = 0;
    size_t nb_quads  = 0;
    size_t nb_edges  = 0;

    std::ofstream file;
    file.open( (file_path).c_str(), std::ios::out );

    // header

    file << "$MeshFormat\n";
    file << "2.2 0 8\n";
    file << "$EndMeshFormat\n";

    // nodes

    FunctionSpace& nodes   = mesh.function_space( "nodes" );
    ArrayView<double,2> coords  ( nodes.field( "coordinates" ) );
    ArrayView<int,   1> glb_idx ( nodes.field( "glb_idx" ) );

    nb_nodes = nodes.extents()[0];

    file << "$Nodes\n";
    file << nb_nodes << "\n";

    for( size_t n = 0; n < nb_nodes; ++n )
    {
        const double x = coords(n,XX);
        const double y = coords(n,YY);
        const double z = coords(n,ZZ);

        file << glb_idx(n) << " " << x << " " << y << " " << z << "\n";

    }

    file << "$EndNodes\n";
    file << std::flush;

    file << "$Elements\n";

    if( mesh.has_function_space("triags") ) nb_triags = mesh.function_space( "triags" ).extents()[0];
    if( mesh.has_function_space("quads") )  nb_quads = mesh.function_space( "quads" ).extents()[0];
    if( mesh.has_function_space("edges") )  nb_edges = mesh.function_space( "edges" ).extents()[0];

    file << nb_triags + nb_quads + nb_edges << "\n";

    // triags

    if( mesh.has_function_space("triags") )
    {
        FunctionSpace& triags      = mesh.function_space( "triags" );
        IndexView<int,2> triag_nodes ( triags.field( "nodes" ) );

        for( size_t e = 0; e < nb_triags; ++e )
        {
            file << e << " 2 2 1 1";
            for( int n=0; n<3; ++n )
                file << " " << triag_nodes(e,n);
            file << "\n";
        }
    }

    if( mesh.has_function_space("quads") )
    {
        FunctionSpace& quads      = mesh.function_space( "quads" );
        IndexView<int,2> quad_nodes ( quads.field( "nodes" ) );

        for( int e=0; e<nb_quads; ++e)
        {
          file << e << " 3 2 1 1";
          for( int n=0; n<4; ++n )
            file << " " << quad_nodes(e,n);
          file << "\n";
        }
    }

    if( mesh.has_function_space("edges") )
    {
        FunctionSpace& edges      = mesh.function_space( "edges" );
        IndexView<int,2> edge_nodes ( edges.field( "nodes" ) );

        for( int e=0; e<nb_edges; ++e)
        {
            file << e << " 1 2 2 1";
            for( int n=0; n<2; ++n )
                file << " " << edge_nodes(e,n);
            file << "\n";
        }
    }

    file << "$EndElements\n";
    file << std::flush;

    // nodal data
    for( size_t fidx = 0; fidx < nodes.nb_fields(); ++fidx )
    {
        Field& field = nodes.field(fidx);

        if( field.data_type() == "real64" )
        {
            ArrayView<double,2> f ( field );

            for( size_t idx = 0; idx < field.extents()[1]; ++idx )
            {
                file << "$NodeData\n";
                file << "1\n";
                file << "\""+field.name() << "[" << idx << "]" << "\"\n";
                file << "1\n";
                file << "0.\n";
                file << "3\n";
                file << "0\n";
                file << "1\n";
                file << nb_nodes << "\n";
                for( size_t n = 0; n < nb_nodes; ++n )
                  file << glb_idx(n) << " " << f(n,idx) <<"\n";
                file << "$EndNodeData\n";
                file << std::flush;
            }
        }
    }
    file.close();
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
  return This->read( std::string(file_path) );
}

void atlas__Gmsh__write (Gmsh* This, Mesh* mesh, char* file_path) {
  This->write( *mesh, std::string(file_path) );
}

Mesh* atlas__read_gmsh (char* file_path)
{
  return Gmsh::read(std::string(file_path));
}

void atlas__write_gmsh (Mesh* mesh, char* file_path) {
  Gmsh::write( *mesh, std::string(file_path) );
}


// ------------------------------------------------------------------


} // namespace atlas

