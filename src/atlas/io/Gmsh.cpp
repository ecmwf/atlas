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

#include <eckit/filesystem/LocalPathName.h>

#include "atlas/io/Gmsh.hpp"
#include "atlas/mpl/GatherScatter.hpp"
#include "atlas/util/Array.hpp"
#include "atlas/util/ArrayView.hpp"
#include "atlas/util/IndexView.hpp"
#include "atlas/mesh/Mesh.hpp"
#include "atlas/mesh/FunctionSpace.hpp"
#include "atlas/mesh/Field.hpp"
#include "atlas/mesh/FieldSet.hpp"
#include "atlas/mesh/Parameters.hpp"

namespace atlas {

namespace {

double deg_to_rad = M_PI / 180.;
double rad_to_deg = 180. * M_1_PI;

class GmshFile : public std::ofstream
{
public:
  GmshFile(const std::string& file_path, std::ios_base::openmode mode, int part=MPL::rank())
  {
    eckit::LocalPathName par_path(file_path);
    bool is_new_file = (mode != std::ios_base::app || !par_path.exists() );
    if( MPL::size() == 1)
    {
      std::ofstream::open(par_path.c_str(), mode );
    }
    else
    {
      eckit::Translator<int,std::string> to_str;
      if( MPL::rank() == 0 )
      {
        eckit::LocalPathName par_path(file_path);
        std::ofstream par_file( par_path.c_str(), std::ios_base::out );
        for (int p=0; p<MPL::size(); ++p)
        {
          eckit::LocalPathName loc_path(file_path);
          loc_path = loc_path.dirName()+"/"+loc_path.baseName(false)+"_p"+to_str(p)+".msh";
          par_file << "Merge \"" << loc_path << "\";" << std::endl;
        }
        par_file.close();
      }
      eckit::LocalPathName path(file_path);
      path = path.dirName()+"/"+path.baseName(false)+"_p"+to_str(part)+".msh";
      std::ofstream::open(path.c_str(), mode );
    }
  }
};

enum GmshElementTypes { LINE=1, TRIAG=2, QUAD=3, POINT=15 };

void write_header(std::ostream& out)
{
  out << "$MeshFormat\n";
  out << "2.2 0 8\n";
  out << "$EndMeshFormat\n";
}

template< typename DATA_TYPE >
void write_field_nodes(const Gmsh& gmsh, Field& field, std::ostream& out)
{
  FunctionSpace& function_space = field.function_space();
  bool gather( gmsh.options.get<bool>("gather") );

  int ndata = field.extents()[0];
  int nvars = field.extents()[1];
  ArrayView<int,1    > gidx ( function_space.field( "glb_idx" ) );
  ArrayView<DATA_TYPE> data ( field );
  if( gather )
  {
    mpl::GatherScatter::Ptr gather_scatter = function_space.gather_scatter();
    function_space.parallelise();
    ndata = gather_scatter->glb_dof();
    DEBUG_VAR(gather_scatter->glb_dof());
    DEBUG_VAR(gather_scatter->loc_dof());
    Array<DATA_TYPE> field_glb_arr(Extents(ndata,nvars));
    Array<int      > gidx_glb_arr (Extents(ndata));
    ArrayView<DATA_TYPE> data_glb( field_glb_arr );
    ArrayView<int,1> gidx_glb( gidx_glb_arr );
    gather_scatter->gather( gidx, gidx_glb );
    gather_scatter->gather( data, data_glb );
    gidx = gidx_glb;
    data = data_glb;
  }

  if( gather && MPL::rank() == 0 || !gather )
  {
    double time   = field.metadata().has<double>("time") ? field.metadata().get<double>("time") : 0. ;
    int time_step = field.metadata().has<int>("time_step") ? field.metadata().get<int>("time_step") : 0. ;
    out << "$NodeData\n";
    out << "1\n";
    out << "\"" << field.name() << "\"\n";
    out << "1\n";
    out << time << "\n";
    out << "4\n";
    out << time_step << "\n";
    if     ( nvars == 1 ) out << nvars << "\n";
    else if( nvars <= 3 ) out << 3     << "\n";
    out << ndata << "\n";
    out << MPL::rank() << "\n";

    if( nvars == 1)
    {
      for( int n = 0; n < ndata; ++n )
      {
        out << gidx(n) << " " << data(n,0) << "\n";
      }
    }
    else if( nvars <= 3 )
    {
      std::vector<DATA_TYPE> data_vec(3,0.);
      for( size_t n = 0; n < ndata; ++n )
      {
        out << gidx(n);
        for( int v=0; v<nvars; ++v)
          data_vec[v] = data(n,v);
        for( int v=0; v<3; ++v)
          out << " " << data_vec[v];
        out << "\n";
      }
    }
    out << "$EndNodeData\n";
  }
}

template< typename DATA_TYPE >
void write_field_elems(const Gmsh& gmsh, Field& field, std::ostream& out)
{
  FunctionSpace& function_space = field.function_space();
  int ndata = field.extents()[0];
  int nvars = field.extents()[1];
  double time   = field.metadata().has<double>("time") ? field.metadata().get<double>("time") : 0. ;
  int time_step = field.metadata().has<int>("time_step") ? field.metadata().get<int>("time_step") : 0. ;
  int nnodes = IndexView<int,2>( function_space.field("nodes") ).extents()[1];
  out << "$ElementNodeData\n";
  out << "1\n";
  out << "\"" << field.name() << "\"\n";
  out << "1\n";
  out << time << "\n";
  out << "4\n";
  out << time_step << "\n";
  if     ( nvars == 1 ) out << nvars << "\n";
  else if( nvars <= 3 ) out << 3     << "\n";
  out << ndata << "\n";
  out << MPL::rank() << "\n";

  ArrayView<int,1> gidx ( function_space.field( "glb_idx" ) );
  if( nvars == 1)
  {
    ArrayView<DATA_TYPE,1> data( field );
    for (size_t jelem=0; jelem<ndata; ++jelem)
    {
      out << gidx(jelem) << " " << nnodes;
      for (size_t n=0; n<nnodes; ++n)
        out << " " << data(jelem);
      out <<"\n";
    }
  }
  else if( nvars <= 3 )
  {
    std::vector<DATA_TYPE> data_vec(3,0.);
    ArrayView<DATA_TYPE,2> data( field );
    for (size_t jelem=0; jelem<ndata; ++jelem)
    {
      out << gidx(jelem) << " " << nnodes;
      for (size_t n=0; n<nnodes; ++n)
      {
        for( int v=0; v<nvars; ++v)
          data_vec[v] = data(jelem,v);
        for( int v=0; v<3; ++v)
          out << " " << data_vec[v];
      }
      out <<"\n";
    }
  }
  out << "$EndElementNodeData\n";
}

} // end anonymous namespace

Gmsh::Gmsh()
{
  options.set<bool>("gather",true);
}

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


void Gmsh::write(Mesh& mesh, const std::string& file_path) const
{
  int part = MPL::rank();
  if( mesh.metadata().has<int>("part") )
    part = mesh.metadata().get<int>("part");

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
    eckit::LocalPathName path(file_path);
    if( spherical )
      path = path.dirName()+"/"+path.baseName(false)+"_sphere.msh";

    if( MPL::rank() == 0 ) std::cout << "writing mesh to gmsh file " << path << std::endl;
    GmshFile file(path,std::ios_base::out,part);

    // Header
    write_header(file);

    // Nodes
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

    // Elements
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

  eckit::LocalPathName mesh_info(file_path);
  mesh_info = mesh_info.dirName()+"/"+mesh_info.baseName(false)+"_info.msh";

  if (nodes.has_field("partition"))
  {
    write(nodes.field("partition"),mesh_info,std::ios_base::out);
  }

  if (nodes.has_field("dual_volumes"))
  {
    write(nodes.field("dual_volumes"),mesh_info,std::ios_base::app);
  }

  if( mesh.has_function_space("edges") )
  {
    FunctionSpace& edges       = mesh.function_space( "edges" );
    if (edges.has_field("dual_normals"))
    {
      write(edges.field("dual_normals"),mesh_info,std::ios_base::app);
    }

    if (edges.has_field("skewness"))
    {
      write(edges.field("skewness"),mesh_info,std::ios_base::app);
    }
  }
}

void Gmsh::write(FieldSet& fieldset, const std::string& file_path, openmode mode) const
{
  eckit::LocalPathName path(file_path);
  bool is_new_file = (mode != std::ios_base::app || !path.exists() );
  GmshFile file(path,mode);

  if( MPL::rank() == 0) std::cout << "writing fieldset " << fieldset.name() << " to gmsh file " << path << std::endl;

  // Header
  if( is_new_file )
    write_header(file);

  // Fields
  for( int field_idx=0; field_idx<fieldset.size(); ++field_idx )
  {
    Field& field = fieldset.field(field_idx);
    FunctionSpace& function_space = field.function_space();
    if( function_space.metadata().get<int>("type") == Entity::NODES )
    {
      if     ( field.data_type() == "int32"  ) {  write_field_nodes<int   >(*this,field,file); }
      else if( field.data_type() == "real32" ) {  write_field_nodes<float >(*this,field,file); }
      else if( field.data_type() == "real64" ) {  write_field_nodes<double>(*this,field,file); }
    }
    else if( function_space.metadata().get<int>("type") == Entity::ELEMS )
    {
      if     ( field.data_type() == "int32"  ) {  write_field_elems<int   >(*this,field,file); }
      else if( field.data_type() == "real32" ) {  write_field_elems<float >(*this,field,file); }
      else if( field.data_type() == "real64" ) {  write_field_elems<double>(*this,field,file); }
    }
    file << std::flush;
  }

  file.close();
}

void Gmsh::write(Field& field, const std::string& file_path, openmode mode) const
{
  eckit::LocalPathName path(file_path);
  bool is_new_file = (mode != std::ios_base::app || !path.exists() );

  GmshFile file(path,mode);

  if( MPL::rank() == 0) std::cout << "writing field " << field.name() << " to gmsh file " << path << std::endl;

  // Header
  if( is_new_file )
    write_header(file);

  // Field
  FunctionSpace& function_space = field.function_space();

  if( function_space.metadata().get<int>("type") == Entity::NODES )
  {
    if     ( field.data_type() == "int32"  ) {  write_field_nodes<int   >(*this,field,file); }
    else if( field.data_type() == "real32" ) {  write_field_nodes<float >(*this,field,file); }
    else if( field.data_type() == "real64" ) {  write_field_nodes<double>(*this,field,file); }
  }
  else if( function_space.metadata().get<int>("type") == Entity::ELEMS ||
           function_space.metadata().get<int>("type") == Entity::FACES )
  {
    if     ( field.data_type() == "int32"  ) {  write_field_elems<int   >(*this,field,file); }
    else if( field.data_type() == "real32" ) {  write_field_elems<float >(*this,field,file); }
    else if( field.data_type() == "real64" ) {  write_field_elems<double>(*this,field,file); }
  }
  file << std::flush;
  file.close();
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
    write_header(file);

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

void atlas__write_gmsh_mesh (Mesh* mesh, char* file_path) {
  Gmsh writer;
  writer.write( *mesh, std::string(file_path) );
}

void atlas__write_gmsh_fieldset (FieldSet* fieldset, char* file_path, int mode) {
  Gmsh writer;
  writer.write( *fieldset, std::string(file_path) );
}

void atlas__write_gmsh_field (Field* field, char* file_path, int mode) {
  Gmsh writer;
  writer.write( *field, std::string(file_path) );
}

// ------------------------------------------------------------------


} // namespace atlas

