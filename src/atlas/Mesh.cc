/*
 * (C) Copyright 1996-2015 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "atlas/Mesh.h"

#include <sstream>
#include <stdexcept>

#include "eckit/exception/Exceptions.h"
#include "atlas/runtime/Log.h"

#include "atlas/runtime/ErrorHandling.h"
#include "atlas/FunctionSpace.h"
#include "atlas/Field.h"
#include "atlas/Grid.h"
#include "atlas/Parameters.h"
#include "atlas/mesh/Nodes.h"
#include "atlas/mesh/HybridElements.h"

namespace atlas {

//----------------------------------------------------------------------------------------------------------------------

Mesh* Mesh::create( const eckit::Parametrisation& params )
{
  return new Mesh(params);
}

Mesh* Mesh::create( const Grid& grid, const eckit::Parametrisation& params )
{
  return new Mesh(grid,params);
}

Mesh::Mesh( const eckit::Parametrisation& ):
  grid_(NULL), dimensionality_(2)
{
  createElements();
}

Mesh::Mesh(const Grid& grid, const eckit::Parametrisation& ) :
  grid_(&grid), dimensionality_(2)
{
  createNodes(grid);
  createElements();
}

Mesh::~Mesh()
{
}

bool Mesh::has_function_space(const std::string& name) const
{
  ASSERT( name != "nodes" );
  return function_spaces_.has(name);
}

FunctionSpace& Mesh::create_function_space(const std::string& name, const std::string& shape_func, const std::vector<size_t>& shape)
{
  ASSERT( name != "nodes" );
	if( has_function_space(name) )
	{
		throw eckit::Exception( "Functionspace '" + name + "' already exists", Here() );
	}

	FunctionSpace::Ptr fs (new FunctionSpace(name,shape_func,shape,*this) );

  function_spaces_.insert(name,fs);
  function_spaces_.sort();

  fs->set_index( function_spaces_.size() - 1 ); ///< @todo revisit this once we can remove functionspaces

	return *fs;
}

FunctionSpace& Mesh::function_space(const std::string& name) const
{
  ASSERT( name != "nodes" );
	if( ! has_function_space(name) )
	{
		std::stringstream msg;
		msg << "Could not find FunctionSpace '" << name << "' in mesh";
		throw eckit::OutOfRange(msg.str(),Here());
	}
	return *( function_spaces_[ name ] );
}

FunctionSpace& Mesh::function_space( size_t idx) const
{
	if( idx >= function_spaces_.size() )
		throw eckit::OutOfRange(idx,function_spaces_.size(),Here());
	return *function_spaces_[ idx ];
}


mesh::Nodes& Mesh::createNodes(const Grid& g)
{
  set_grid(g);
  size_t nb_nodes = g.npts();
  mesh::Nodes& nodes = createNodes(nb_nodes);

  g.fillLonLat(nodes.lonlat().data<double>(), nb_nodes*2);
  return nodes;
}

void Mesh::prettyPrint(std::ostream& os) const
{
    os << "Mesh:" << std::endl;
    for( size_t i = 0; i < nb_function_spaces(); ++i )
    {
        os << function_space(i) << std::endl;
    }
}
void Mesh::print(std::ostream& os) const
{
    os << "Mesh[";
    for( size_t i = 0; i < nb_function_spaces(); ++i )
    {
        os << function_space(i);
    }
    os << "]";
}

void Mesh::createElements()
{
  cells_.reset( new mesh::HybridElements() );
  facets_.reset( new mesh::HybridElements() );
  ridges_.reset( new mesh::HybridElements() );
  peaks_.reset( new mesh::HybridElements() );
  if( dimensionality_ == 2 )
    edges_ = facets_;
  else if( dimensionality_ == 3)
    edges_ = ridges_;
  else
    throw eckit::Exception("Invalid Mesh dimensionality",Here());
}

mesh::Nodes& Mesh::createNodes( size_t size )
{
  if( nodes_ )
  {
    Log::error() << "ERROR: Re-creating nodes.\n"
                        << "This error can be ignored in MIR/prodgen\n"
                        << "until MIR stops using deprecated atlas::Grid::mesh() function."
                        << std::endl;
    //throw eckit::SeriousBug("Cannot create nodes in mesh as they already exist.", Here());
  }
  nodes_.reset( new mesh::Nodes(size) );
  return *nodes_;
}


size_t Mesh::nb_function_spaces() const
{
  return function_spaces_.size();
}


//----------------------------------------------------------------------------------------------------------------------

// C wrapper interfaces to C++ routines

Mesh* atlas__Mesh__new () {
	return new Mesh();
}

void atlas__Mesh__delete (Mesh* This) {
	delete This;
}

void atlas__Mesh__create_function_space(Mesh* This, char* name,char* shape_func,int shape[], int shape_size, int fortran_ordering)
{
  ATLAS_ERROR_HANDLING(

    ASSERT( This != NULL );
    ASSERT( shape_size >= 0 );

    std::vector<size_t> vshape(shape_size);

    if( fortran_ordering )
    {
      size_t r=vshape.size()-1;
      for(size_t n=0; n<vshape.size(); ++n)
      {
        vshape[n]=shape[r];
        --r;
      };
    }
    else
    {
      for(size_t n=0; n<vshape.size(); ++n)
          vshape[n]=shape[n];
    }
    if( std::string(name) == "nodes" )
      This->createNodes(vshape[0]);
    else
      This->create_function_space(std::string(name), std::string(shape_func),vshape);
  );
}

FunctionSpace* atlas__Mesh__function_space (Mesh* This, char* name) {
  ATLAS_ERROR_HANDLING(
    ASSERT( This != NULL );
    if( std::string(name) == "nodes" )
      throw eckit::BadParameter("nodes is no longer a *old* FunctionSpace",Here());
    else
      return &This->function_space( std::string(name) );
  );
  return NULL;
}

mesh::Nodes* atlas__Mesh__create_nodes (Mesh* This, int nb_nodes)
{
  ATLAS_ERROR_HANDLING(
    ASSERT( This );
    return &This->createNodes(nb_nodes);
  );
  return NULL;
}

mesh::Nodes* atlas__Mesh__nodes (Mesh* This) {
  ATLAS_ERROR_HANDLING(
    ASSERT( This != NULL );
    return &This->nodes();
  );
  return NULL;
}

//----------------------------------------------------------------------------------------------------------------------

} // namespace atlas




#include "atlas/Mesh.h"
#include "atlas/FunctionSpace.h"
#include "atlas/Field.h"
#include "atlas/util/IndexView.h"
#include "atlas/util/Debug.h"
#include "atlas/Parameters.h"
#include "atlas/mesh/ElementType.h"
#include "atlas/mesh/HybridElements.h"
#include "atlas/mesh/Elements.h"
#include "atlas/Metadata.h"

namespace atlas {
namespace mesh {
namespace temporary {

void Convert::convertMesh( Mesh& mesh )
{
  mesh.createElements();
  createElements(mesh,Entity::ELEMS,mesh.cells());
  createElements(mesh,Entity::FACES,mesh.facets());
}
  
void Convert::createElements( const Mesh& mesh, Entity::Type type, HybridElements& elements )
{
  HybridElements* hybrid_elements = &elements;
  
  std::set<std::string> fields;
  for(size_t func_space_idx = 0; func_space_idx < mesh.nb_function_spaces(); ++func_space_idx)
  {
    FunctionSpace& functionspace = mesh.function_space(func_space_idx);
    if( functionspace.metadata().get<long>("type") == type )
    {
      size_t nb_elems = functionspace.field("nodes").shape(0);
      size_t nb_nodes_per_elem = functionspace.field("nodes").shape(1);
      ElementType* element_type(0);
      if( nb_nodes_per_elem==2 ) element_type = new Line();
      if( nb_nodes_per_elem==3 ) element_type = new Triangle();
      if( nb_nodes_per_elem==4 ) element_type = new Quadrilateral();
      if(!element_type) throw eckit::SeriousBug("element_type not supported",Here());
    
      bool fortran_array = true;
      size_t t = hybrid_elements->add(element_type,nb_elems,functionspace.field("nodes").data<idx_t>(),fortran_array);
      Elements& elements = hybrid_elements->elements(t);
      
      if( functionspace.has_field("to_edge") )
      {
        size_t nb_edges_per_elem = functionspace.field("to_edge").shape(1);
        hybrid_elements->edge_connectivity().add(nb_elems,nb_edges_per_elem,functionspace.field("to_edge").data<idx_t>(),fortran_array);
        ASSERT( nb_elems == functionspace.field("to_edge").shape(0) );
      }
      
      if( functionspace.has_field("to_elem") )
      {
        size_t nb_connected_cells = 2;
        std::vector<idx_t> to_elem(nb_elems*2);
        const IndexView<int,3> to_elem_field ( functionspace.field( "to_elem" ).data<int>(), make_shape(nb_elems,2,2) );
        for( size_t e=0; e<nb_elems; ++e )
        {
          for( size_t j=0; j<2; ++j )
          {
            int fs = to_elem_field(e,j,0);
            if( fs == -1 )
              to_elem[e*2+j] = hybrid_elements->cell_connectivity().missing_value();
            else if( fs == 0 )
              to_elem[e*2+j] = to_elem_field(e,j,1) + 0;
            else if( fs == 1 )
              to_elem[e*2+j] = to_elem_field(e,j,1) + mesh.function_space("quads").shape(0);
          }
        }
        hybrid_elements->cell_connectivity().add(nb_elems,nb_connected_cells,to_elem.data(),false);
        ASSERT( nb_elems == functionspace.field("to_elem").shape(0) );
      }
      
      std::map<std::string,std::string> rename;
      std::set<std::string> ignore_connectivities;
      rename["glb_idx"]    = hybrid_elements->global_index().name();
      rename["remote_idx"] = hybrid_elements->remote_index().name();
      rename["partition"]  = hybrid_elements->partition().name();
      ignore_connectivities.insert("nodes");
      ignore_connectivities.insert("to_edge");
      ignore_connectivities.insert("to_elem");
            
      for( size_t f=0; f<functionspace.nb_fields(); ++f )
      {
        const Field& field = functionspace.field(f);
        std::string fname = field.name();
        if( ignore_connectivities.find(fname) == ignore_connectivities.end() )
        {
          if( rename.find(fname) != rename.end() )
            fname = rename[fname];
          fields.insert(fname);
          
          eckit::Log::debug(2) << fname << "<"<< field.datatype().str() << ">" << std::endl;
          
          if( !hybrid_elements->has_field(fname) )
          {
            hybrid_elements->add( Field::create(fname, Array::create(field.array())) );
          }
          else
          {
            if( field.rank() == 1 )
            {
              if( field.datatype().kind() == DataType::KIND_REAL64 ) {
                ArrayView<double,1> old_data ( field );
                ArrayView<double,1> new_data ( hybrid_elements->field(fname) );
                for( size_t e=0; e<nb_elems; ++e ) new_data(e + elements.begin()) = old_data(e);
              }
              if( field.datatype().kind() == DataType::KIND_REAL32 ) {
                ArrayView<float,1> old_data ( field );
                ArrayView<float,1> new_data ( hybrid_elements->field(fname) );
                for( size_t e=0; e<nb_elems; ++e ) new_data(e + elements.begin()) = old_data(e);
              }
              if( field.datatype().kind() == DataType::KIND_INT64 ) {
                ArrayView<long,1> old_data ( field );
                ArrayView<long,1> new_data ( hybrid_elements->field(fname) );
                for( size_t e=0; e<nb_elems; ++e ) new_data(e + elements.begin()) = old_data(e);
              }
              if( field.datatype().kind() == DataType::KIND_INT32 ) {
                ArrayView<int,1> old_data ( field );
                ArrayView<int,1> new_data ( hybrid_elements->field(fname) );
                for( size_t e=0; e<nb_elems; ++e ) new_data(e + elements.begin()) = old_data(e);
              }
            }
            if( field.rank() == 2 )
            {
              if( field.datatype().kind() == DataType::KIND_REAL64 ) {
                ArrayView<double,2> old_data ( field );
                ArrayView<double,2> new_data ( hybrid_elements->field(fname) );
                for( size_t e=0; e<nb_elems; ++e ) {
                  for( size_t v=0; v<field.shape(1); ++v ) {
                    new_data(e + elements.begin(),v) = old_data(e,v);
                  }
                }
              }
              if( field.datatype().kind() == DataType::KIND_REAL32 ) {
                ArrayView<float,2> old_data ( field );
                ArrayView<float,2> new_data ( hybrid_elements->field(fname) );
                for( size_t e=0; e<nb_elems; ++e ) {
                  for( size_t v=0; v<field.shape(1); ++v ) {
                    new_data(e + elements.begin(),v) = old_data(e,v);
                  }
                }
              }
              if( field.datatype().kind() == DataType::KIND_INT64 ) {
                ArrayView<long,2> old_data ( field );
                ArrayView<long,2> new_data ( hybrid_elements->field(fname) );
                for( size_t e=0; e<nb_elems; ++e ) {
                  for( size_t v=0; v<field.shape(1); ++v ) {
                    new_data(e + elements.begin(),v) = old_data(e,v);
                  }
                }
              }
              if( field.datatype().kind() == DataType::KIND_INT32 ) {
                ArrayView<int,2> old_data ( field );
                ArrayView<int,2> new_data ( hybrid_elements->field(fname) );
                for( size_t e=0; e<nb_elems; ++e ) {
                  for( size_t v=0; v<field.shape(1); ++v ) {
                    new_data(e + elements.begin(),v) = old_data(e,v);
                  }
                }
              }
            }
          }
        }
      
      }
      
      if( functionspace.metadata().has("nb_owned") )
      {
        size_t nb_owned = functionspace.metadata().get<size_t>("nb_owned");
        ArrayView<int,1> halo ( hybrid_elements->halo() );
        for( size_t e=0; e<nb_owned; ++e )
          halo(elements.begin()+e) = 0;
        for( size_t e=nb_owned; e<functionspace.shape(0); ++e )
          halo(elements.begin()+e) = 1;
      }
      else
      {
        ArrayView<int,1> halo ( hybrid_elements->halo() );
        for( size_t e=0; e<functionspace.shape(0); ++e )
          halo(elements.begin()+e) = 0;
      }
    }
  }
}


}
}
}
