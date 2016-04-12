/*
 * (C) Copyright 1996-2016 ECMWF.
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
#include "eckit/filesystem/PathName.h"
#include "eckit/exception/Exceptions.h"
#include "eckit/config/Resource.h"
#include "eckit/runtime/Context.h"
#include "atlas/mesh/Nodes.h"
#include "atlas/mesh/HybridElements.h"
#include "atlas/mesh/ElementType.h"
#include "atlas/mesh/Elements.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/field/Field.h"
#include "atlas/field/FieldSet.h"
#include "atlas/functionspace/FunctionSpace.h"
#include "atlas/util/Constants.h"
#include "atlas/internals/Parameters.h"
#include "atlas/util/io/Gmsh.h"
#include "atlas/parallel/GatherScatter.h"
#include "atlas/array/Array.h"
#include "atlas/array/ArrayView.h"
#include "atlas/array/IndexView.h"
#include "atlas/runtime/Log.h"

using namespace eckit;
using atlas::functionspace::NodeColumns;

namespace atlas {
namespace util {
namespace io {

namespace {

static double deg = util::Constants::radiansToDegrees();
static double rad = util::Constants::degreesToRadians();

class GmshFile : public std::ofstream {
public:
  GmshFile(const PathName& file_path, std::ios_base::openmode mode, int part = eckit::mpi::rank())
  {
    PathName par_path(file_path);
    if (eckit::mpi::size() == 1 || part == -1) {
      std::ofstream::open(par_path.localPath(), mode);
    } else {
      Translator<int, std::string> to_str;
      if (eckit::mpi::rank() == 0) {
        PathName par_path(file_path);
        std::ofstream par_file(par_path.localPath(), std::ios_base::out);
        for(size_t p = 0; p < eckit::mpi::size(); ++p) {
          PathName loc_path(file_path);
          loc_path = loc_path.baseName(false) + "_p" + to_str(p) + ".msh";
          par_file << "Merge \"" << loc_path << "\";" << std::endl;
        }
        par_file.close();
      }
      PathName path(file_path);
      path = path.dirName() + "/" + path.baseName(false) + "_p" + to_str(part) + ".msh";
      std::ofstream::open(path.localPath(), mode);
    }
  }
};

enum GmshElementTypes { LINE=1, TRIAG=2, QUAD=3, POINT=15 };



// ----------------------------------------------------------------------------
void write_header_ascii(std::ostream& out)
{
  out << "$MeshFormat\n";
  out << "2.2 0 "<<sizeof(double)<<"\n";
  out << "$EndMeshFormat\n";
}
// ----------------------------------------------------------------------------



// ----------------------------------------------------------------------------
void write_header_binary(std::ostream& out)
{
  out << "$MeshFormat\n";
  out << "2.2 1 "<<sizeof(double)<<"\n";
  int one = 1;
  out.write(reinterpret_cast<const char*>(&one),sizeof(int));
  out << "\n$EndMeshFormat\n";
}
// ----------------------------------------------------------------------------



// ----------------------------------------------------------------------------
template< typename DATATYPE >
void write_field_nodes(const Gmsh& gmsh, const functionspace::NodeColumns& function_space, const field::Field& field, std::ostream& out)
{
  Log::info() << "writing field " << field.name() << "..." << std::endl;

  bool gather( gmsh.options.get<bool>("gather") );
  bool binary( !gmsh.options.get<bool>("ascii") );
  int nlev  = field.levels();
  int ndata = field.shape(0);
  int nvars = field.stride(0)/nlev;
  array::ArrayView<gidx_t,1> gidx ( function_space.nodes().global_index() );
  array::ArrayView<DATATYPE,2> data ( field.data<DATATYPE>(), array::make_shape(field.shape(0),field.stride(0)) );
  field::Field::Ptr gidx_glb;
  field::Field::Ptr data_glb;
  if( gather )
  {
    gidx_glb.reset( function_space.createGlobalField( "gidx_glb", function_space.nodes().global_index() ) );
    function_space.gather(function_space.nodes().global_index(),*gidx_glb);
    gidx = array::ArrayView<gidx_t,1>( *gidx_glb );

    data_glb.reset( function_space.createGlobalField( "glb_field",field ) );
    function_space.gather(field,*data_glb);
    data = array::ArrayView<DATATYPE,2>( data_glb->data<DATATYPE>(), array::make_shape(data_glb->shape(0),data_glb->stride(0)) );
  }

  std::vector<long> lev;
  std::vector<long> gmsh_levels;
  gmsh.options.get("levels",gmsh_levels);
  if( gmsh_levels.empty() || nlev == 1 )
  {
    lev.resize(nlev);
    for (int ilev=0; ilev<nlev; ++ilev)
      lev[ilev] = ilev;
  }
  else
  {
    lev = gmsh_levels;
  }
  for (int ilev=0; ilev<lev.size(); ++ilev)
  {
    int jlev = lev[ilev];
    if( ( gather && eckit::mpi::rank() == 0 ) || !gather )
    {
      char field_lev[6] = {0, 0, 0, 0, 0, 0};
      if( field.has_levels() )
        std::sprintf(field_lev, "[%03d]",jlev);
      double time = field.metadata().has("time") ? field.metadata().get<double>("time") : 0.;
      int step = field.metadata().has("step") ? field.metadata().get<size_t>("step") : 0 ;
      out << "$NodeData\n";
      out << "1\n";
      out << "\"" << field.name() << field_lev << "\"\n";
      out << "1\n";
      out << time << "\n";
      out << "4\n";
      out << step << "\n";
      if     ( nvars == 1 ) out << nvars << "\n";
      else if( nvars <= 3 ) out << 3     << "\n";
      out << ndata << "\n";
      out << eckit::mpi::rank() << "\n";

      if( binary )
      {
        if( nvars == 1)
        {
          double value;
          for( int n = 0; n < ndata; ++n )
          {
            out.write(reinterpret_cast<const char*>(&gidx(n)),sizeof(int));
            value = data(n,jlev*nvars+0);
            out.write(reinterpret_cast<const char*>(&value),sizeof(double));
          }
        }
        else if( nvars <= 3 )
        {
          double value[3] = {0,0,0};
          for( size_t n = 0; n < ndata; ++n )
          {
            out.write(reinterpret_cast<const char*>(&gidx(n)),sizeof(int));
            for( int v=0; v<nvars; ++v)
              value[v] = data(n,jlev*nvars+v);
            out.write(reinterpret_cast<const char*>(&value),sizeof(double)*3);
          }
        }
        out << "\n";
      }
      else
      {
        if( nvars == 1)
        {
          for( int n = 0; n < ndata; ++n )
          {
            ASSERT( jlev*nvars < data.shape(1) );
            ASSERT( n < gidx.shape(0) );
            out << gidx(n) << " " << data(n,jlev*nvars+0) << "\n";
          }
        }
        else if( nvars <= 3 )
        {
          std::vector<DATATYPE> data_vec(3,0.);
          for( size_t n = 0; n < ndata; ++n )
          {
            out << gidx(n);
            for( int v=0; v<nvars; ++v)
              data_vec[v] = data(n,jlev*nvars+v);
            for( int v=0; v<3; ++v)
              out << " " << data_vec[v];
            out << "\n";
          }
        }

      }
      out << "$EndNodeData\n";
    }
  }
}
// ----------------------------------------------------------------------------



// ----------------------------------------------------------------------------
template< typename DATATYPE >
void write_field_nodes(
    const Gmsh&                            gmsh,
    const functionspace::StructuredColumns& function_space,
    const field::Field&                           field,
    std::ostream&                          out)
{
    Log::info() << "writing field " << field.name() << "..." << std::endl;

    //bool gather(gmsh.options.get<bool>("gather"));
    bool binary(!gmsh.options.get<bool>("ascii"));

    int nlev  = field.levels();
    int nvars = field.stride(0) / nlev;

    //field::Field::Ptr gidxField(function_space.createField<gidx_t>("gidx"));
    //array::ArrayView<gidx_t,1>   gidx(gidxField);
    array::ArrayView<DATATYPE,2> data(field.data<DATATYPE>(),
                               array::make_shape(field.shape(0),
                                          field.stride(0)));

    field::Field::Ptr gidx_glb;
    field::Field::Ptr data_glb;

    //gidx_glb.reset(function_space.createGlobalField(
    //    "gidx_glb", function_space.nodes().global_index()));

    //function_space.gather(
    //    field.global_index(), *gidx_glb);
    //array::ArrayView<gidx_t,1> gidx = array::ArrayView<gidx_t,1>(*gidx_glb);

    data_glb = function_space.createGlobalField<double>("glb_field");
    function_space.gather(field, *data_glb);
    data = array::ArrayView<DATATYPE,2>(data_glb->data<DATATYPE>(),
                                 array::make_shape(data_glb->shape(0),
                                            data_glb->stride(0)));

    int ndata = data_glb->shape(0);

    std::vector<long> lev;
    std::vector<long> gmsh_levels;
    gmsh.options.get("levels", gmsh_levels);

    if (gmsh_levels.empty() || nlev == 1)
    {
        lev.resize(nlev);
        for (int ilev=0; ilev<nlev; ++ilev)
        {
            lev[ilev] = ilev;
        }
    }
    else
    {
        lev = gmsh_levels;
    }

    if (eckit::mpi::rank() == 0)
    {
        for (int ilev = 0; ilev < lev.size(); ++ilev)
        {
            int jlev = lev[ilev];
            char field_lev[6] = {0, 0, 0, 0, 0, 0};

            if (field.has_levels())
            {
                std::sprintf(field_lev, "[%03d]", jlev);
            }

            double time = field.metadata().has("time") ?
                field.metadata().get<double>("time") : 0.;

            int step = field.metadata().has("step") ?
                field.metadata().get<size_t>("step") : 0 ;

            out << "$NodeData\n";
            out << "1\n";
            out << "\"" << field.name() << field_lev << "\"\n";
            out << "1\n";
            out << time << "\n";
            out << "4\n";
            out << step << "\n";
            if     ( nvars == 1 ) out << nvars << "\n";
            else if( nvars <= 3 ) out << 3     << "\n";
            out << ndata << "\n";
            out << eckit::mpi::rank() << "\n";

            if (binary)
            {
                if (nvars == 1)
                {
                    double value;
                    for (int n = 0; n < ndata; ++n)
                    {
                        out.write(reinterpret_cast<const char*>(n+1),
                                  sizeof(int));

                        value = data(n,jlev*nvars+0);

                        out.write(reinterpret_cast<const char*>(&value),
                                  sizeof(double));
                    }
                }
                else if (nvars <= 3)
                {
                    double value[3] = {0,0,0};
                    for (size_t n = 0; n < ndata; ++n)
                    {
                        out.write(reinterpret_cast<const char*>(n+1),
                                  sizeof(int));
                        for (int v = 0; v < nvars; ++v)
                        {
                            value[v] = data(n,jlev*nvars+v);
                        }
                        out.write(reinterpret_cast<const char*>(&value),
                                  sizeof(double)*3);
                    }
                }
                out << "\n";
            }
            else
            {
                if (nvars == 1)
                {
                    for (int n = 0; n < ndata; ++n)
                    {
                        ASSERT(jlev*nvars < data_glb->shape(1));
                        ASSERT(n < data_glb->shape(0));
                        out << n+1 << " "
                            << data(n, jlev*nvars+0) << "\n";
                    }
                }
                else if (nvars <= 3)
                {
                    std::vector<DATATYPE> data_vec(3,0.);
                    for (size_t n = 0; n < ndata; ++n)
                    {
                        out << n+1;
                        for (int v = 0; v < nvars; ++v)
                        {
                            data_vec[v] = data(n, jlev*nvars+v);
                        }
                        for (int v = 0; v < 3; ++v)
                        {
                            out << " " << data_vec[v];
                        }
                        out << "\n";
                    }
                }
            }
            out << "$EndNodeData\n";
        }
    }
}
// ----------------------------------------------------------------------------



// ----------------------------------------------------------------------------
#if 0
template< typename DATA_TYPE >
void write_field_elems(const Gmsh& gmsh, const FunctionSpace& function_space, const field::Field& field, std::ostream& out)
{
  Log::info() << "writing field " << field.name() << "..." << std::endl;
  bool gather( gmsh.options.get<bool>("gather") );
  bool binary( !gmsh.options.get<bool>("ascii") );
  int nlev = field.metadata().has("nb_levels") ? field.metadata().get<size_t>("nb_levels") : 1;
  int ndata = field.shape(0);
  int nvars = field.shape(1)/nlev;
  array::ArrayView<gidx_t,1    > gidx ( function_space.field( "glb_idx" ) );
  array::ArrayView<DATA_TYPE> data ( field );
  array::ArrayT<DATA_TYPE> field_glb_arr;
  array::ArrayT<gidx_t   > gidx_glb_arr;
  if( gather )
  {
    mpl::GatherScatter& fullgather = function_space.fullgather();
    ndata = fullgather.glb_dof();
    field_glb_arr.resize(ndata,field.shape(1));
    gidx_glb_arr.resize(ndata);
    array::ArrayView<DATA_TYPE> data_glb( field_glb_arr );
    array::ArrayView<gidx_t,1> gidx_glb( gidx_glb_arr );
    fullgather.gather( gidx, gidx_glb );
    fullgather.gather( data, data_glb );
    gidx = array::ArrayView<gidx_t,1>( gidx_glb_arr );
    data = data_glb;
  }

  double time = field.metadata().has("time") ? field.metadata().get<double>("time") : 0.;
  size_t step = field.metadata().has("step") ? field.metadata().get<size_t>("step") : 0 ;

  int nnodes = IndexView<int,2>( function_space.field("nodes") ).shape(1);

  for (int jlev=0; jlev<nlev; ++jlev)
  {
    char field_lev[6] = {0, 0, 0, 0, 0, 0};
    if( field.metadata().has("nb_levels") )
      std::sprintf(field_lev, "[%03d]",jlev);

    out << "$ElementNodeData\n";
    out << "1\n";
    out << "\"" << field.name() << field_lev << "\"\n";
    out << "1\n";
    out << time << "\n";
    out << "4\n";
    out << step << "\n";
    if     ( nvars == 1 ) out << nvars << "\n";
    else if( nvars <= 3 ) out << 3     << "\n";
    out << ndata << "\n";
    out << eckit::mpi::rank() << "\n";

    if( binary )
    {
      if( nvars == 1)
      {
        double value;
        for (size_t jelem=0; jelem<ndata; ++jelem)
        {
          out.write(reinterpret_cast<const char*>(&gidx(jelem)),sizeof(int));
          out.write(reinterpret_cast<const char*>(&nnodes),sizeof(int));
          for (size_t n=0; n<nnodes; ++n)
          {
            value = data(jelem,jlev);
            out.write(reinterpret_cast<const char*>(&value),sizeof(double));
          }
        }
      }
      else if( nvars <= 3 )
      {
        double value[3] = {0,0,0};
        for (size_t jelem=0; jelem<ndata; ++jelem)
        {
          out << gidx(jelem) << " " << nnodes;
          for (size_t n=0; n<nnodes; ++n)
          {
            for( int v=0; v<nvars; ++v)
              value[v] = data(jelem,jlev*nvars+v);
            out.write(reinterpret_cast<const char*>(&value),sizeof(double)*3);
          }
        }
      }
      out <<"\n";
    }
    else
    {
      if( nvars == 1)
      {
        for (size_t jelem=0; jelem<ndata; ++jelem)
        {
          out << gidx(jelem) << " " << nnodes;
          for (size_t n=0; n<nnodes; ++n)
            out << " " << data(jelem,jlev);
          out <<"\n";
        }
      }
      else if( nvars <= 3 )
      {
        std::vector<DATA_TYPE> data_vec(3,0.);
        for (size_t jelem=0; jelem<ndata; ++jelem)
        {
          out << gidx(jelem) << " " << nnodes;
          for (size_t n=0; n<nnodes; ++n)
          {
            for( int v=0; v<nvars; ++v)
              data_vec[v] = data(jelem,jlev*nvars+v);
            for( int v=0; v<3; ++v)
              out << " " << data_vec[v];
          }
          out <<"\n";
        }
      }
    }
    out << "$EndElementNodeData\n";
  }
}
#endif
// ----------------------------------------------------------------------------



// ----------------------------------------------------------------------------
// Unused private function, in case for big-endian
// ----------------------------------------------------------------------------
#if 0
void swap_bytes(char *array, int size, int n)
{
  char *x = new char[size];
  for(int i = 0; i < n; i++) {
    char *a = &array[i * size];
    memcpy(x, a, size);
    for(int c = 0; c < size; c++)
      a[size - 1 - c] = x[c];
  }
  delete [] x;
}
#endif
// ----------------------------------------------------------------------------
} // end anonymous namespace



// ----------------------------------------------------------------------------
Gmsh::Gmsh()
{
  // which field holds the Nodes
  options.set<std::string>("nodes", Resource<std::string>("atlas.gmsh.nodes", "lonlat"));

  // Gather fields to one proc before writing
  options.set<bool>("gather",  Resource<bool>("atlas.gmsh.gather",  false));

  // Output of ghost nodes / elements
  options.set<bool>("ghost",   Resource<bool>("atlas.gmsh.ghost",   false));

  // ASCII format (true) or binary (false)
  options.set<bool>("ascii",   Resource<bool>("atlas.gmsh.ascii",   true ));

  // Output of elements
  options.set<bool>("elements",Resource<bool>("atlas.gmsh.elements",true ));

  // Output of edges
  options.set<bool>("edges",   Resource<bool>("atlas.gmsh.edges",   true ));

  // Radius of the planet
  options.set<double>("radius",   Resource<bool>("atlas.gmsh.radius", 1.0 ));

  // Levels of fields to use
  options.set< std::vector<long> >("levels", Resource< std::vector<long> >("atlas.gmsh.levels", std::vector<long>() ) );
}

Gmsh::~Gmsh()
{
}

mesh::Mesh* Gmsh::read(const PathName& file_path) const
{
  mesh::Mesh* mesh = new mesh::Mesh();
  Gmsh::read(file_path,*mesh);
  return mesh;
}

namespace {
mesh::ElementType* make_element_type(int type)
{
  if( type == QUAD  ) return new mesh::temporary::Quadrilateral();
  if( type == TRIAG ) return new mesh::temporary::Triangle();
  if( type == LINE  ) return new mesh::temporary::Line();
  throw eckit::SeriousBug("Element type not supported", Here());
  return 0;
}
}

void Gmsh::read(const PathName& file_path, mesh::Mesh& mesh ) const
{
  std::ifstream file;
  file.open( file_path.localPath() , std::ios::in | std::ios::binary );
  if( !file.is_open() )
    throw CantOpenFile(file_path);

  std::string line;

  while(line != "$MeshFormat")
    std::getline(file,line);
  double version;
  int binary;
  int size_of_real;
  file >> version >> binary >> size_of_real;

  while(line != "$Nodes")
    std::getline(file,line);

  // Create nodes
    size_t nb_nodes;
  file >> nb_nodes;

  mesh.nodes().resize(nb_nodes);

  mesh::Nodes& nodes = mesh.nodes();

  nodes.add( field::Field::create<double>("xyz",array::make_shape(nb_nodes,3) ) );

  array::ArrayView<double,2> coords         ( nodes.field("xyz")    );
  array::ArrayView<gidx_t,1> glb_idx        ( nodes.global_index()  );
  array::ArrayView<int,   1> part           ( nodes.partition()     );

  std::map<int,int> glb_to_loc;
  int g;
  double x,y,z;
  double xyz[3];
  double xmax = -std::numeric_limits<double>::max();
  double zmax = -std::numeric_limits<double>::max();
  gidx_t max_glb_idx=0;
  while(binary && file.peek()=='\n') file.get();
  for( size_t n = 0; n < nb_nodes; ++n )
  {
    if( binary )
    {
      file.read(reinterpret_cast<char*>(&g), sizeof(int));
      file.read(reinterpret_cast<char*>(&xyz), sizeof(double)*3);
      x = xyz[internals::XX];
      y = xyz[internals::YY];
      z = xyz[internals::ZZ];
    }
    else
    {
      file >> g >> x >> y >> z;
    }
    glb_idx(n) = g;
    coords(n,internals::XX) = x;
    coords(n,internals::YY) = y;
    coords(n,internals::ZZ) = z;
    glb_to_loc[ g ] = n;
    part(n) = 0;
    max_glb_idx = std::max(max_glb_idx, static_cast<gidx_t>(g));
    xmax = std::max(x,xmax);
    zmax = std::max(z,zmax);
  }
  if( xmax < 4*M_PI && zmax == 0. )
  {
    for( size_t n = 0; n < nb_nodes; ++n )
    {
      coords(n,internals::XX) *= deg;
      coords(n,internals::YY) *= deg;
    }
  }
  for (int i=0; i<3; ++i)
    std::getline(file,line);

  int nb_elements=0;

  while(line != "$Elements")
    std::getline(file,line);

  file >> nb_elements;

  if( binary )
  {
    while(file.peek()=='\n') file.get();
    int accounted_elems = 0;
    while( accounted_elems < nb_elements )
    {
      int header[3];
      int data[100];
      file.read(reinterpret_cast<char*>(&header),sizeof(int)*3);

      int etype = header[0];
      size_t netype = header[1];
      size_t ntags = header[2];
      accounted_elems += netype;
      mesh::Elements* elements;
      if( etype == LINE ) {
        size_t jtype = mesh.edges().add(make_element_type(etype),netype);
        elements = &mesh.edges().elements(jtype);
      } else {
        size_t jtype = mesh.cells().add(make_element_type(etype),netype);
        elements = &mesh.edges().elements(jtype);
      }

      size_t nnodes_per_elem = elements->element_type().nb_nodes();
      mesh::Elements::Connectivity& conn = elements->node_connectivity();
      array::ArrayView<gidx_t,1> egidx ( elements->global_index() );
      array::ArrayView<int,1> epart    ( elements->partition() );

      size_t dsize = 1+ntags+nnodes_per_elem;
      int part;
      for (size_t e=0; e<netype; ++e)
      {
        file.read(reinterpret_cast<char*>(&data),sizeof(int)*dsize);
        part = 0;
        egidx(e) = data[0];
        epart(e) = part;
        for(size_t n=0; n<nnodes_per_elem; ++n)
          conn.set(e,n,glb_to_loc[data[1+ntags+n]]);
      }
    }
  }
  else
  {
    // Find out which element types are inside
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

    int nb_quads  = nb_etype[QUAD];
    int nb_triags = nb_etype[TRIAG];
    int nb_edges = nb_etype[LINE];


    mesh::Elements& quads  = mesh.cells().elements( mesh.cells().add( make_element_type(QUAD),  nb_quads  ) );
    mesh::Elements& triags = mesh.cells().elements( mesh.cells().add( make_element_type(TRIAG), nb_triags ) );
    mesh::Elements& edges  = mesh.edges().elements( mesh.edges().add( make_element_type(LINE),  nb_edges ) );

    mesh::Elements::Connectivity& quad_nodes  = quads.node_connectivity();
    mesh::Elements::Connectivity& triag_nodes = triags.node_connectivity();
    mesh::Elements::Connectivity& edge_nodes  = edges.node_connectivity();

    array::ArrayView<gidx_t,1> quad_glb_idx   ( quads.global_index() );
    array::ArrayView<int,1> quad_part         ( quads.partition()    );

    array::ArrayView<gidx_t,1> triag_glb_idx  ( triags.global_index() );
    array::ArrayView<int,1> triag_part        ( triags.partition()    );

    array::ArrayView<gidx_t,1> edge_glb_idx   ( edges.global_index() );
    array::ArrayView<int,1> edge_part         ( edges.partition()    );

    // Now read all elements
    file.seekg(position,std::ios::beg);
    int gn0, gn1, gn2, gn3;
    int quad=0, triag=0, edge=0;
    int ntags, tags[100];
    for (int e=0; e<nb_elements; ++e)
    {
      file >> g >> etype >> ntags;
      for( int t=0; t<ntags; ++t ) file >> tags[t];
      int part=0;
      if( ntags > 3 ) part = std::max( part, *std::max_element(tags+3,tags+ntags-1) ); // one positive, others negative

      int enodes[4];

      switch( etype )
      {
        case(QUAD):
          file >> gn0 >> gn1 >> gn2 >> gn3;
          quad_glb_idx(quad) = g;
          quad_part(quad) = part;
          enodes[0] = glb_to_loc[gn0];
          enodes[1] = glb_to_loc[gn1];
          enodes[2] = glb_to_loc[gn2];
          enodes[3] = glb_to_loc[gn3];
          quad_nodes.set(quad,enodes);
          ++quad;
          break;
        case(TRIAG):
          file >> gn0 >> gn1 >> gn2;
          triag_glb_idx(triag) = g;
          triag_part(triag) = part;
          enodes[0] = glb_to_loc[gn0];
          enodes[1] = glb_to_loc[gn1];
          enodes[2] = glb_to_loc[gn2];
          triag_nodes.set(triag,enodes);
          ++triag;
          break;
        case(LINE):
          file >> gn0 >> gn1;
          edge_glb_idx(edge) = g;
          edge_part(edge) = part;
          enodes[0] = glb_to_loc[gn0];
          enodes[1] = glb_to_loc[gn1];
          edge_nodes.set(edge,enodes);
          ++edge;
          break;
        case(POINT):
          file >> gn0;
          break;
        default:
          std::cout << "etype " << etype << std::endl;
          throw Exception("ERROR: element type not supported",Here());
      }
    }
  }

  file.close();
}


void Gmsh::write(const mesh::Mesh& mesh, const PathName& file_path) const
{
  int part = mesh.metadata().has("part") ? mesh.metadata().get<size_t>("part") : eckit::mpi::rank();
  bool include_ghost = options.get<bool>("ghost") && options.get<bool>("elements");

  std::string nodes_field = options.get<std::string>("nodes");

  const mesh::Nodes& nodes    = mesh.nodes();
  array::ArrayView<double,2> coords  ( nodes.field( nodes_field ) );
  array::ArrayView<gidx_t,   1> glb_idx ( nodes.global_index() );

  const size_t surfdim = coords.shape(1); // nb of variables in coords

  ASSERT(surfdim == 2 || surfdim == 3);

  Log::info() << "writing mesh to gmsh file " << file_path << std::endl;

  bool binary = !options.get<bool>("ascii");

  openmode mode = std::ios::out;
  if( binary )
    mode = std::ios::out | std::ios::binary;
  GmshFile file(file_path,mode,part);

  // Header
  if( binary )
    write_header_binary(file);
  else
    write_header_ascii(file);

  // Nodes
  const size_t nb_nodes = nodes.size();
  file << "$Nodes\n";
  file << nb_nodes << "\n";
  double xyz[3] = {0.,0.,0.};
  for( size_t n = 0; n < nb_nodes; ++n )
  {
    int g = glb_idx(n);

    for(size_t d = 0; d < surfdim; ++d)
        xyz[d] = coords(n,d);

    if( binary )
    {
        file.write(reinterpret_cast<const char*>(&g), sizeof(int));
        file.write(reinterpret_cast<const char*>(&xyz), sizeof(double)*3 );
    }
    else
    {
        file << g << " " << xyz[internals::XX] << " " << xyz[internals::YY] << " " << xyz[internals::ZZ] << "\n";
    }
  }
  if( binary ) file << "\n";
  file << "$EndNodes\n";

  // Elements
  file << "$Elements\n";
  {
    std::vector<const mesh::HybridElements*> grouped_elements;
    if( options.get<bool>("elements") )
      grouped_elements.push_back(&mesh.cells());
    if( options.get<bool>("edges") )
      grouped_elements.push_back(&mesh.edges());

    size_t nb_elements(0);
    for( size_t jgroup=0; jgroup<grouped_elements.size(); ++jgroup )
    {
      const mesh::HybridElements& hybrid = *grouped_elements[jgroup];
      nb_elements += hybrid.size();
      if( !include_ghost )
      {
        const array::ArrayView<int,1> hybrid_halo( hybrid.halo() );
        for( size_t e=0; e<hybrid.size(); ++e )
        {
          if( hybrid_halo(e) ) --nb_elements;
        }
      }
    }

    file << nb_elements << "\n";

    for( size_t jgroup=0; jgroup<grouped_elements.size(); ++jgroup )
    {
      const mesh::HybridElements& hybrid = *grouped_elements[jgroup];
      for( size_t etype=0; etype<hybrid.nb_types(); ++etype )
      {
        const mesh::Elements& elements = hybrid.elements(etype);
        const mesh::ElementType& element_type = elements.element_type();
        int gmsh_elem_type;
        if( element_type.name() == "Line" )
          gmsh_elem_type = 1;
        else if( element_type.name() == "Triangle" )
          gmsh_elem_type = 2;
        else if( element_type.name() == "Quadrilateral" )
          gmsh_elem_type = 3;
        else
          NOTIMP;

        const mesh::Elements::Connectivity& node_connectivity = elements.node_connectivity();
        const array::ArrayView<gidx_t,1> elems_glb_idx = elements.view<gidx_t,1>( elements.global_index() );
        const array::ArrayView<int,1> elems_partition = elements.view<int,1>( elements.partition() );
        const array::ArrayView<int,1> elems_halo = elements.view<int,1>( elements.halo() );

        if( binary )
        {
          size_t nb_elems=elements.size();
          if( !include_ghost )
          {
            for(size_t elem=0; elem<elements.size(); ++elem )
            {
              if( elems_halo(elem) ) --nb_elems;
            }
          }

          int header[3];
          int data[9];
          header[0] = gmsh_elem_type;
          header[1] = nb_elems;
          header[2] = 4; // nb_tags
          file.write(reinterpret_cast<const char*>(&header), sizeof(int)*3 );
          data[1]=1;
          data[2]=1;
          data[3]=1;
          size_t datasize = sizeof(int)*(5+node_connectivity.cols());
          for(size_t elem=0; elem<elements.size(); ++elem )
          {
            if( include_ghost || !elems_halo(elem) )
            {
              data[0] = elems_glb_idx(elem);
              data[4] = elems_partition(elem);
              for( int n=0; n<node_connectivity.cols(); ++n )
                data[5+n] = glb_idx(node_connectivity(elem,n));
              file.write(reinterpret_cast<const char*>(&data), datasize );
            }
          }
        }
        else
        {
          std::stringstream ss_elem_info; ss_elem_info << " " << gmsh_elem_type << " 4 1 1 1 ";
          std::string elem_info = ss_elem_info.str();
          for(size_t elem=0; elem<elements.size(); ++elem )
          {
            if( include_ghost || !elems_halo(elem) )
            {
              file << elems_glb_idx(elem) << elem_info << elems_partition(+elem);
              for( int n=0; n<node_connectivity.cols(); ++n ) {
                file << " " << glb_idx(node_connectivity(elem,n));
              }
              file << "\n";
            }
          }
        }
      }
    }
  }
  if( binary ) file << "\n";
  file << "$EndElements\n";
  file << std::flush;
  file.close();

  // Optional mesh information file
  if( options.has("info") && options.get<bool>("info") )
  {
    PathName mesh_info(file_path);
    mesh_info = mesh_info.dirName()+"/"+mesh_info.baseName(false)+"_info.msh";

    //[next]  make NodesFunctionSpace accept const mesh
    eckit::SharedPtr<functionspace::NodeColumns> function_space( new functionspace::NodeColumns(const_cast<mesh::Mesh&>(mesh)) );

    write(nodes.partition(),*function_space,mesh_info,std::ios_base::out);

    if (nodes.has_field("dual_volumes"))
    {
      write(nodes.field("dual_volumes"),*function_space,mesh_info,std::ios_base::app);
    }

    if (nodes.has_field("dual_delta_sph"))
    {
      write(nodes.field("dual_delta_sph"),*function_space,mesh_info,std::ios_base::app);
    }

    //[next] if( mesh.has_function_space("edges") )
    //[next] {
    //[next]   FunctionSpace& edges = mesh.function_space( "edges" );

    //[next]   if (edges.has_field("dual_normals"))
    //[next]   {
    //[next]     write(edges.field("dual_normals"),mesh_info,std::ios_base::app);
    //[next]   }

    //[next]   if (edges.has_field("skewness"))
    //[next]   {
    //[next]     write(edges.field("skewness"),mesh_info,std::ios_base::app);
    //[next]   }

    //[next]   if (edges.has_field("arc_length"))
    //[next]   {
    //[next]     write(edges.field("arc_length"),mesh_info,std::ios_base::app);
    //[next]   }
    //[next] }
  }

}

// ----------------------------------------------------------------------------
void Gmsh::write(
    const field::Field&    field,
    const PathName& file_path,
    openmode        mode) const
{
    if (!field.functionspace())
    {
        std::stringstream msg;
        msg << "Field ["<<field.name()<<"] has no functionspace";

        throw eckit::AssertionFailed(msg.str(), Here());
    }

    if ( field.functionspace().cast<functionspace::NodeColumns>() )
    {
        const functionspace::NodeColumns* functionspace =
            field.functionspace().cast<functionspace::NodeColumns>();

        field::FieldSet fieldset;
        fieldset.add(field);
        write(fieldset, *functionspace, file_path, mode);
    }
    else if ( field.functionspace().cast<functionspace::StructuredColumns>() )
    {
        const functionspace::StructuredColumns* functionspace =
            field.functionspace().cast<functionspace::StructuredColumns>();

        field::FieldSet fieldset;
        fieldset.add(field);
        write(fieldset, *functionspace, file_path, mode);
    }
    else
    {
        std::stringstream msg;
        msg << "Field ["<<field.name()
            <<"] has functionspace ["
            << field.functionspace().name()
            << "] but requires a [functionspace::NodeColumns "
            << "or functionspace::StructuredColumns]";

        throw eckit::AssertionFailed(msg.str(), Here());
    }
}
// ----------------------------------------------------------------------------



// ----------------------------------------------------------------------------
void Gmsh::write(
    const field::Field&                field,
    const functionspace::NodeColumns& functionspace,
    const PathName&             file_path,
    openmode                    mode) const
{
  field::FieldSet fieldset;
  fieldset.add(field);
  write(fieldset, functionspace, file_path, mode);
}
// ----------------------------------------------------------------------------



// ----------------------------------------------------------------------------
void Gmsh::write(
    const field::Field&                           field,
    const functionspace::StructuredColumns& functionspace,
    const PathName&                        file_path,
    openmode                               mode) const
{
    field::FieldSet fieldset;
    fieldset.add(field);
    write(fieldset, functionspace, file_path, mode);
}
// ----------------------------------------------------------------------------



// ----------------------------------------------------------------------------
void Gmsh::write(
    const field::FieldSet&             fieldset,
    const functionspace::NodeColumns& functionspace,
    const PathName&             file_path,
    openmode                    mode) const
{
    bool is_new_file = (mode != std::ios_base::app || !file_path.exists() );
    bool binary( !options.get<bool>("ascii") );
    if ( binary ) mode |= std::ios_base::binary;
    bool gather = options.has("gather") ? options.get<bool>("gather") : false;
    GmshFile file(file_path,mode,gather?-1:eckit::mpi::rank());

    // Header
    if (is_new_file)
    {
        write_header_ascii(file);
    }

    // field::Fields
    for(size_t field_idx = 0; field_idx < fieldset.size(); ++field_idx)
    {
        const field::Field& field = fieldset[field_idx];
        Log::info() << "writing field " << field.name()
                    << " to gmsh file " << file_path << std::endl;

        if (field.datatype() == array::DataType::int32())
        {
            write_field_nodes<int   >(*this,functionspace,field,file);
        }
        else if (field.datatype() == array::DataType::int64())
        {
            write_field_nodes<long  >(*this,functionspace,field,file);
        }
        else if (field.datatype() == array::DataType::real32())
        {
            write_field_nodes<float >(*this,functionspace,field,file);
        }
        else if (field.datatype() == array::DataType::real64())
        {
            write_field_nodes<double>(*this,functionspace,field,file);
        }

        file << std::flush;
    }
    file.close();
}
// ----------------------------------------------------------------------------



// ----------------------------------------------------------------------------
void Gmsh::write(
    const field::FieldSet& fieldset,
    const functionspace::StructuredColumns& functionspace,
    const PathName& file_path, openmode mode) const
{
    bool is_new_file = (mode != std::ios_base::app || !file_path.exists());
    bool binary(!options.get<bool>("ascii"));

    if (binary) mode |= std::ios_base::binary;

    bool gather = options.has("gather") ? options.get<bool>("gather") : false;

    GmshFile file(file_path,mode,gather?-1:eckit::mpi::rank());

    // Header
    if (is_new_file)
        write_header_ascii(file);

    // field::Fields
    for (size_t field_idx = 0; field_idx < fieldset.size(); ++field_idx)
    {
        const field::Field& field = fieldset[field_idx];
        Log::info() << "writing field " << field.name()
                    << " to gmsh file " << file_path << std::endl;

        if (field.datatype() == array::DataType::int32())
        {
            write_field_nodes<int   >(*this, functionspace, field, file);
        }
        else if (field.datatype() == array::DataType::int64())
        {
            write_field_nodes<long  >(*this, functionspace, field, file);
        }
        else if (field.datatype() == array::DataType::real32())
        {
            write_field_nodes<float >(*this, functionspace, field, file);
        }
        else if (field.datatype() == array::DataType::real64())
        {
            write_field_nodes<double>(*this, functionspace, field, file);
        }

        file << std::flush;
    }

    file.close();
}
// ----------------------------------------------------------------------------




// ----------------------------------------------------------------------------
// C wrapper interfaces to C++ routines
// ----------------------------------------------------------------------------
Gmsh* atlas__Gmsh__new () {
  return new Gmsh();
}

void atlas__Gmsh__delete (Gmsh* This) {
  delete This;
}

mesh::Mesh* atlas__Gmsh__read (Gmsh* This, char* file_path) {
  return This->read( PathName(file_path) );
}

void atlas__Gmsh__write (Gmsh* This, mesh::Mesh* mesh, char* file_path) {
  This->write( *mesh, PathName(file_path) );
}

mesh::Mesh* atlas__read_gmsh (char* file_path)
{
  return Gmsh().read(PathName(file_path));
}

void atlas__write_gmsh_mesh (mesh::Mesh* mesh, char* file_path) {
  Gmsh writer;
  writer.write( *mesh, PathName(file_path) );
}

void atlas__write_gmsh_fieldset (field::FieldSet* fieldset, functionspace::NodeColumns* functionspace, char* file_path, int mode) {
  Gmsh writer;
  writer.write( *fieldset, *functionspace, PathName(file_path) );
}

void atlas__write_gmsh_field (field::Field* field, functionspace::NodeColumns* functionspace, char* file_path, int mode) {
  Gmsh writer;
  writer.write( *field, *functionspace, PathName(file_path) );
}
// ----------------------------------------------------------------------------

} // namespace io
} // namespace util
} // namespace atlas

