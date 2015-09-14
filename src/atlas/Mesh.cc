/*
 * (C) Copyright 1996-2014 ECMWF.
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
#include "eckit/log/Log.h"

#include "atlas/runtime/ErrorHandling.h"
#include "atlas/FunctionSpace.h"
#include "atlas/Field.h"
#include "atlas/Grid.h"
#include "atlas/Parameters.h"
#include "atlas/Nodes.h"

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
  grid_(NULL)
{
}

Mesh::Mesh(const Grid& grid, const eckit::Parametrisation& ) :
  grid_(&grid)
{
  createNodes(grid);
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


Nodes& Mesh::createNodes(const Grid& g)
{
  set_grid(g);
  size_t nb_nodes = g.npts();
  Nodes& nodes = createNodes(nb_nodes);

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

Nodes& Mesh::createNodes( size_t size )
{
  if( nodes_ )
  {
    eckit::Log::error() << "ERROR: Re-creating nodes.\n"
                        << "This error can be ignored in MIR/prodgen\n"
                        << "until MIR stops using deprecated atlas::Grid::mesh() function."
                        << std::endl;
    //throw eckit::SeriousBug("Cannot create nodes in mesh as they already exist.", Here());
  }
  nodes_.reset( new Nodes(size) );
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

Nodes* atlas__Mesh__create_nodes (Mesh* This, int nb_nodes)
{
  ATLAS_ERROR_HANDLING(
    ASSERT( This );
    return &This->createNodes(nb_nodes);
  );
  return NULL;
}

Nodes* atlas__Mesh__nodes (Mesh* This) {
  ATLAS_ERROR_HANDLING(
    ASSERT( This != NULL );
    return &This->nodes();
  );
  return NULL;
}

//----------------------------------------------------------------------------------------------------------------------

} // namespace atlas

