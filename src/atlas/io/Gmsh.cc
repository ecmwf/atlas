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

#include "atlas/io/Gmsh.h"
#include "atlas/mpl/GatherScatter.h"
#include "atlas/util/Array.h"
#include "atlas/util/ArrayView.h"
#include "atlas/util/IndexView.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/mesh/FunctionSpace.h"
#include "atlas/mesh/Field.h"
#include "atlas/mesh/FieldSet.h"
#include "atlas/mesh/Parameters.h"

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

void write_header_ascii(std::ostream& out)
{
  out << "$MeshFormat\n";
  out << "2.2 0 "<<sizeof(double)<<"\n";
  out << "$EndMeshFormat\n";
}
void write_header_binary(std::ostream& out)
{
  out << "$MeshFormat\n";
  out << "2.2 1 "<<sizeof(double)<<"\n";
	int one = 1;
	out.write(reinterpret_cast<const char*>(&one),sizeof(int));
  out << "\n$EndMeshFormat\n";
}

template< typename DATA_TYPE >
void write_field_nodes(const Gmsh& gmsh, Field& field, std::ostream& out)
{
  FunctionSpace& function_space = field.function_space();
  bool gather( gmsh.options.get<bool>("gather") );
	bool binary( !gmsh.options.get<bool>("ascii") );

  int ndata = field.shape(0);
  int nvars = field.shape(1);
  ArrayView<int,1    > gidx ( function_space.field( "glb_idx" ) );
  ArrayView<DATA_TYPE> data ( field );
  if( gather )
  {
    mpl::GatherScatter& gather_scatter = function_space.gather_scatter();
    ndata = gather_scatter.glb_dof();
    Array<DATA_TYPE> field_glb_arr(ndata,nvars);
    Array<int      > gidx_glb_arr (ndata);
    ArrayView<DATA_TYPE> data_glb( field_glb_arr );
    ArrayView<int,1> gidx_glb( gidx_glb_arr );
    gather_scatter.gather( gidx, gidx_glb );
    gather_scatter.gather( data, data_glb );
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

		if( binary )
		{
			if( nvars == 1)
	    {
				double value;
	      for( int n = 0; n < ndata; ++n )
	      {
					out.write(reinterpret_cast<const char*>(&gidx(n)),sizeof(int));
					value = data(n,0);
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
	          value[v] = data(n,v);
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

		}
    out << "$EndNodeData\n";
  }
}

template< typename DATA_TYPE >
void write_field_elems(const Gmsh& gmsh, Field& field, std::ostream& out)
{
	bool binary( !gmsh.options.get<bool>("ascii") );

	FunctionSpace& function_space = field.function_space();
	int ndata = field.shape(0);
	int nvars = field.shape(1);
	double time   = field.metadata().has<double>("time") ? field.metadata().get<double>("time") : 0. ;
	int time_step = field.metadata().has<int>("time_step") ? field.metadata().get<int>("time_step") : 0. ;
	int nnodes = IndexView<int,2>( function_space.field("nodes") ).shape(1);
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

	if( binary )
	{
		if( nvars == 1)
		{
			ArrayView<DATA_TYPE,1> data( field );
			double value;
			for (size_t jelem=0; jelem<ndata; ++jelem)
			{
				out.write(reinterpret_cast<const char*>(&gidx(jelem)),sizeof(int));
				out.write(reinterpret_cast<const char*>(&nnodes),sizeof(int));
				for (size_t n=0; n<nnodes; ++n)
				{
					value = data(jelem);
					out.write(reinterpret_cast<const char*>(&value),sizeof(double));
				}
			}
		}
		else if( nvars <= 3 )
		{
			double value[3] = {0,0,0};
			ArrayView<DATA_TYPE,2> data( field );
			for (size_t jelem=0; jelem<ndata; ++jelem)
			{
				out << gidx(jelem) << " " << nnodes;
				for (size_t n=0; n<nnodes; ++n)
				{
					for( int v=0; v<nvars; ++v)
						value[v] = data(jelem,v);
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
	}
	out << "$EndElementNodeData\n";
}

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


} // end anonymous namespace

Gmsh::Gmsh()
{
  options.set<bool>("gather",false);
  options.set<int>("surfdim",2);    // lonlat
  options.set<bool>("ghost",true);
	options.set<bool>("ascii",true);
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
	file.open( file_path.c_str() , std::ios::in | std::ios::binary );
	if( !file.is_open() )
		throw std::runtime_error("Could not open file "+file_path);

	std::string line;

	while(line != "$MeshFormat")
		std::getline(file,line);
	double version;
	int binary;
	int size_of_real;
	file >> version >> binary >> size_of_real;

//	if( binary )
//	{
//		throw eckit::NotImplemented("Reading of binary format not implemented",Here());
//	}

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
			if( mesh.function_space("nodes").shape(0)!= nb_nodes )
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
	double xyz[3];
	double xmax = -std::numeric_limits<double>::max();
	double zmax = -std::numeric_limits<double>::max();
	int max_glb_idx=0;
	bool swap = true;
	while(binary && file.peek()=='\n') file.get();
	for( size_t n = 0; n < nb_nodes; ++n )
	{
		if( binary )
		{
			file.read(reinterpret_cast<char*>(&g), sizeof(int));
			file.read(reinterpret_cast<char*>(&xyz), sizeof(double)*3);
			x = xyz[XX];
			y = xyz[YY];
			z = xyz[ZZ];
		}
		else
		{
			file >> g >> x >> y >> z;
		}
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

	int nb_elements=0;
	int nb_quads=0;
	int nb_triags=0;
	int nb_edges=0;

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
			int netype = header[1];
			int ntags = header[2];
			accounted_elems += netype;
			std::string name;
			int nnodes_per_elem;
			switch( etype ) {
			case(QUAD):
				nnodes_per_elem = 4;
				name = "quads";
				nb_quads = netype;
				break;
			case(TRIAG):
				nnodes_per_elem = 3;
				name = "triags";
				nb_triags=netype;
				break;
			case(LINE):
				nnodes_per_elem = 3;
				name = "edges";
				nb_edges = netype;
				break;
			default:
				std::cout << "etype " << etype << std::endl;
				throw std::runtime_error("ERROR: element type not supported");
			}
			extents[0] = netype;
			FunctionSpace& funcspace = mesh.add_function_space( new FunctionSpace( name, "Lagrange_P1", extents ) );
			funcspace.metadata().set("type",static_cast<int>(Entity::ELEMS));
			IndexView<int,2> conn  ( funcspace.create_field<int>("nodes",         4) );
			ArrayView<int,1> egidx ( funcspace.create_field<int>("glb_idx",       1) );
			ArrayView<int,1> epart ( funcspace.create_field<int>("partition",     1) );
			int dsize = 1+ntags+nnodes_per_elem;
			int part;
			for (int e=0; e<netype; ++e)
			{
				file.read(reinterpret_cast<char*>(&data),sizeof(int)*dsize);
				part = 0;
				egidx(e) = data[0];
				epart(e) = part;
				for(int n=0; n<nnodes_per_elem; ++n)
					conn(e,n) = glb_to_loc[data[1+ntags+n]];
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

	}

	file.close();
}


void Gmsh::write(Mesh& mesh, const std::string& file_path) const
{
	int part = MPL::rank();
	if( mesh.metadata().has<int>("part") )
		part = mesh.metadata().get<int>("part");

	bool include_ghost_elements = options.get<bool>("ghost");
	int surfdim = options.get<int>("surfdim");

	FunctionSpace& nodes    = mesh.function_space( "nodes" );
	ArrayView<double,2> coords  ( nodes.field( "coordinates" ) );
	ArrayView<int,   1> glb_idx ( nodes.field( "glb_idx" ) );
	int nb_nodes;
	if( nodes.metadata().has<int>("nb_owned") )
		nb_nodes = nodes.metadata().get<int>("nb_owned");
	nb_nodes = nodes.shape(0);

	FunctionSpace& quads       = mesh.function_space( "quads" );
	IndexView<int,2> quad_nodes   ( quads.field( "nodes" ) );
	ArrayView<int,1> quad_glb_idx ( quads.field( "glb_idx" ) );
	ArrayView<int,1> quad_part    ( quads.field( "partition" ) );
	int nb_quads = quads.shape(0);
	if( quads.metadata().has<int>("nb_owned") )
		nb_quads = quads.metadata().get<int>("nb_owned");

	FunctionSpace& triags      = mesh.function_space( "triags" );
	IndexView<int,2> triag_nodes   ( triags.field( "nodes" ) );
	ArrayView<int,1> triag_glb_idx ( triags.field( "glb_idx" ) );
	ArrayView<int,1> triag_part    ( triags.field( "partition" ) );
	int nb_triags = quads.shape(0);
	if( triags.metadata().has<int>("nb_owned") )
		nb_triags = triags.metadata().get<int>("nb_owned");

	int nb_edges(0);
	if( mesh.has_function_space("edges") )
	{
		FunctionSpace& edges       = mesh.function_space( "edges" );
		nb_edges = edges.shape(0);
	}

	if( include_ghost_elements == true )
	{
		nb_quads = quads.shape(0);
		nb_triags = triags.shape(0);
	}

	eckit::LocalPathName path(file_path);

	if( MPL::rank() == 0 ) std::cout << "writing mesh to gmsh file " << path << std::endl;

	bool binary = !options.get<bool>("ascii");

	openmode mode = std::ios::out;
	if( binary )
		mode = std::ios::out | std::ios::binary;
	GmshFile file(path,mode,part);

	// Header
	if( binary )
		write_header_binary(file);
	else
		write_header_ascii(file);

	// Nodes
	file << "$Nodes\n";
	file << nb_nodes << "\n";
	double xyz[3] = {0.,0.,0.};
	double r      = 1.;
	for( size_t n = 0; n < nb_nodes; ++n )
	{
		double lon = coords(n,XX);
		double lat = coords(n,YY);

		if(surfdim == 2)
		{
			xyz[XX] = lon;
			xyz[YY] = lat;
		}
		if(surfdim == 3)
		{
			xyz[XX] = r*std::cos(lat)*std::cos(lon);
			xyz[YY] = r*std::cos(lat)*std::sin(lon);
			xyz[ZZ] = r*std::sin(lat);
		}
		if( binary )
		{
			file.write(reinterpret_cast<const char*>(&glb_idx(n)), sizeof(int));
			file.write(reinterpret_cast<const char*>(&xyz), sizeof(double)*3 );
		}
		else
		{
			file << glb_idx(n) << " " << xyz[XX] << " " << xyz[YY] << " " << xyz[ZZ] << "\n";
		}
	}
	if( binary ) file << "\n";
	file << "$EndNodes\n";

	// Elements
	file << "$Elements\n";

	if( binary)
	{
		file << nb_quads+nb_triags+nb_edges << "\n";
		int header[3];
		int data[9];
		header[0] = 3;         // elm_type = QUAD
		header[1] = nb_quads;  // nb_elems
		header[2] = 4;         // nb_tags
		file.write(reinterpret_cast<const char*>(&header), sizeof(int)*3 );
		data[1]=1;
		data[2]=1;
		data[3]=1;
		for( int e=0; e<nb_quads; ++e)
		{
			data[0] = quad_glb_idx(e);
			data[4] = quad_part(e);
			for( int n=0; n<4; ++n )
				data[5+n] = glb_idx( quad_nodes(e,n) );
			file.write(reinterpret_cast<const char*>(&data), sizeof(int)*9 );
		}
		header[0] = 2;         // elm_type = TRIAG
		header[1] = nb_triags; // nb_elems
		header[2] = 4;         // nb_tags
		file.write(reinterpret_cast<const char*>(&header), sizeof(int)*3 );
		data[1]=1;
		data[2]=1;
		data[3]=1;
		for( int e=0; e<nb_triags; ++e)
		{
			data[0] = triag_glb_idx(e);
			data[4] = triag_part(e);
			for( int n=0; n<3; ++n )
				data[5+n] = glb_idx( triag_nodes(e,n) );
			file.write(reinterpret_cast<const char*>(&data), sizeof(int)*8 );
		}
		if( mesh.has_function_space("edges") )
		{
			FunctionSpace& edges       = mesh.function_space( "edges" );
			IndexView<int,2> edge_nodes   ( edges.field( "nodes" ) );
			ArrayView<int,1> edge_glb_idx ( edges.field( "glb_idx" ) );
			header[0] = 1;         // elm_type = LINE
			header[1] = nb_edges;  // nb_elems
			if( edges.has_field("partition") )
			{
				header[2] = 4;         // nb_tags
				data[1]=1;
				data[2]=1;
				data[3]=1;
				ArrayView<int,1> edge_part ( edges.field( "partition" ) );
				for( int e=0; e<nb_edges; ++e)
				{
					data[0] = edge_glb_idx(e);
					data[4] = edge_part(e);
					for( int n=0; n<2; ++n )
						data[5+n] = glb_idx(edge_nodes(e,n) );
					file.write(reinterpret_cast<const char*>(&data), sizeof(int)*7 );
				}
			}
			else
			{
				header[2] = 2;         // nb_tags
				data[1]=1;
				data[2]=1;
				for( int e=0; e<nb_edges; ++e)
				{
					data[0] = edge_glb_idx(e);
					file << edge_glb_idx(e) << " 1 2 1 1";
					for( int n=0; n<2; ++n )
						data[3+n] = glb_idx(edge_nodes(e,n) );
					file.write(reinterpret_cast<const char*>(&data), sizeof(int)*5 );
				}
			}
		}

		file << "\n";
	}
	else
	{
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
	}
	file << "$EndElements\n";
	file << std::flush;
	file.close();

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
		FunctionSpace& edges = mesh.function_space( "edges" );

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
		write_header_ascii(file);

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
	bool binary( !options.get<bool>("ascii") );
	if ( binary ) mode |= std::ios_base::binary;
	GmshFile file(path,mode);

	if( MPL::rank() == 0) std::cout << "writing field " << field.name() << " to gmsh file " << path << std::endl;

	// Header
	if( is_new_file )
	{
		if( binary )
			write_header_binary(file);
		else
			write_header_ascii(file);
	}

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
	if( MPL::rank() == 0) std::cout << "done writing field " << field.name() << " to gmsh file " << path << std::endl;
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
    write_header_ascii(file);

    // nodes

    FunctionSpace& nodes   = mesh.function_space( "nodes" );
    ArrayView<double,2> coords  ( nodes.field( "coordinates" ) );
    ArrayView<int,   1> glb_idx ( nodes.field( "glb_idx" ) );

    nb_nodes = nodes.shape(0);

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

    if( mesh.has_function_space("triags") ) nb_triags = mesh.function_space( "triags" ).shape(0);
    if( mesh.has_function_space("quads") )  nb_quads = mesh.function_space( "quads" ).shape(0);
    if( mesh.has_function_space("edges") )  nb_edges = mesh.function_space( "edges" ).shape(0);

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

            for( size_t idx = 0; idx < field.shape(1); ++idx )
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

