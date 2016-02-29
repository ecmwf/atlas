/*
 * (C) Copyright 1996-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include <sstream>
#include <stdexcept>
#include "eckit/exception/Exceptions.h"
#include "atlas/grid/Grid.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/mesh/Nodes.h"
#include "atlas/mesh/HybridElements.h"
#include "atlas/mesh/Elements.h"
#include "atlas/mesh/ElementType.h"
#include "atlas/field/Field.h"
#include "atlas/functionspace/FunctionSpace.h"
#include "atlas/internals/Parameters.h"
#include "atlas/util/runtime/Log.h"
#include "atlas/util/runtime/ErrorHandling.h"

namespace atlas {
namespace mesh {

//----------------------------------------------------------------------------------------------------------------------

Mesh* Mesh::create( const eckit::Parametrisation& params )
{
  return new Mesh(params);
}

Mesh* Mesh::create( const grid::Grid& grid, const eckit::Parametrisation& params )
{
  return new Mesh(grid,params);
}

Mesh::Mesh( const eckit::Parametrisation& ):
  grid_(NULL), dimensionality_(2)
{
  nodes_.reset( new mesh::Nodes() );
  createElements();
}

Mesh::Mesh(const grid::Grid& grid, const eckit::Parametrisation& ) :
  grid_(&grid), dimensionality_(2)
{
  nodes_.reset( new mesh::Nodes() );
  createNodes(grid);
  createElements();
}

Mesh::~Mesh()
{
}

mesh::Nodes& Mesh::createNodes(const grid::Grid& g)
{
  set_grid(g);
  size_t nb_nodes = g.npts();
  nodes().resize(nb_nodes);
  g.fillLonLat(nodes().lonlat().data<double>(), nb_nodes*2);
  return nodes();
}

void Mesh::prettyPrint(std::ostream& os) const
{
#if !DEPRECATE_OLD_FUNCTIONSPACE
    os << "Mesh:" << std::endl;
    for( size_t i = 0; i < nb_function_spaces(); ++i )
    {
        os << function_space(i) << std::endl;
    }
#endif
}
void Mesh::print(std::ostream& os) const
{
#if !DEPRECATE_OLD_FUNCTIONSPACE
    os << "Mesh[";
    for( size_t i = 0; i < nb_function_spaces(); ++i )
    {
        os << function_space(i);
    }
    os << "]";
#endif
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

  ASSERT( edges_.owners() == 2 );

#if !DEPRECATE_OLD_FUNCTIONSPACE
  // transitional
  edges_->mesh_ = this;
  edges_->type_ = Entity::FACES;
  cells_->mesh_ = this;
  cells_->type_ = Entity::ELEMS;
#endif
}

//----------------------------------------------------------------------------------------------------------------------

#if ! DEPRECATE_OLD_FUNCTIONSPACE

bool deprecated::FunctionSpaceContainer::has_function_space(const std::string& name) const
{
  return index_.count(name);
}

size_t deprecated::FunctionSpaceContainer::nb_function_spaces() const
{
  return function_spaces_.size();
}

deprecated::FunctionSpace& deprecated::FunctionSpaceContainer::create_function_space(const std::string& name, const std::string& shape_func, const std::vector<size_t>& shape)
{
  ASSERT( name != "nodes" );

  if( has_function_space(name) )
	{
    throw eckit::Exception( "Functionspace '" + name + "' already exists", Here() );
	}

  deprecated::FunctionSpace::Ptr fs (new deprecated::FunctionSpace(name,shape_func,shape, *(Mesh*)(this) ) );

  index_[name] = function_spaces_.size();
  function_spaces_.push_back( fs );

  fs->set_index( function_spaces_.size() - 1 );

	return *fs;
}

deprecated::FunctionSpace& deprecated::FunctionSpaceContainer::function_space(const std::string& name) const
{
  ASSERT( name != "nodes" );
	if( ! has_function_space(name) )
	{
		std::stringstream msg;
		msg << "Could not find FunctionSpace '" << name << "' in mesh";
		throw eckit::OutOfRange(msg.str(),Here());
	}
  return *function_spaces_[ index_.at(name) ];
}

deprecated::FunctionSpace& deprecated::FunctionSpaceContainer::function_space( size_t idx) const
{
	if( idx >= function_spaces_.size() )
		throw eckit::OutOfRange(idx,function_spaces_.size(),Here());
	return *function_spaces_[ idx ];
}
#endif

#if ! DEPRECATE_OLD_FUNCTIONSPACE
void Mesh::convert_new_to_old()
{
  std::vector<std::string> functionspace_names;
  std::vector<mesh::Elements*> elements_vec;

  elements_vec.push_back( &cells().elements(0) );
  elements_vec.push_back( &cells().elements(1) );
  functionspace_names.push_back("quads");
  functionspace_names.push_back("triags");

  if( edges().size() ) {
    elements_vec.push_back( &edges().elements(0) );
    functionspace_names.push_back("edges");
  }


  for( size_t jtype=0; jtype<elements_vec.size(); ++jtype )
  {
    const mesh::Elements &elements = *elements_vec[jtype];
    const std::string& name = functionspace_names[jtype];
    size_t nb_elems = elements.size();
    std::vector<size_t> shape = array::make_shape(nb_elems,deprecated::FunctionSpace::UNDEF_VARS);

    if( ! has_function_space(name) )
      create_function_space("name","bla",shape);
    else
      function_space(name).resize(shape);
    deprecated::FunctionSpace& fs = function_space(name);


    for( size_t f=0; f<elements.nb_fields(); ++f )
    {
      const Field& field = elements.field(f);
      std::string fname = field.name();

      Log::debug(2) << fname << "<"<< field.datatype().str() << ">" << std::endl;

      if( !fs.has_field(fname) )
      {
        fs.add( Field::create(fname, Array::create(field.array())) );
      }
      else
      {
        if( field.rank() == 1 )
        {
          if( field.datatype().kind() == DataType::KIND_REAL64 ) {
            array::ArrayView<double,1> data = elements.view<double,1>(field) ;
            array::ArrayView<double,1> old_data ( fs.field(fname) );
            for( size_t e=0; e<nb_elems; ++e ) old_data(e) = data(e);
          }
          if( field.datatype().kind() == DataType::KIND_REAL32 ) {
            array::ArrayView<float,1> data = elements.view<float,1>(field);
            array::ArrayView<float,1> old_data ( fs.field(fname) );
            for( size_t e=0; e<nb_elems; ++e ) old_data(e) = data(e);
          }
          if( field.datatype().kind() == DataType::KIND_INT64 ) {
            array::ArrayView<long,1> data = elements.view<long,1>(field);
            array::ArrayView<long,1> old_data ( fs.field(fname) );
            for( size_t e=0; e<nb_elems; ++e ) old_data(e) = data(e);
          }
          if( field.datatype().kind() == DataType::KIND_INT32 ) {
            array::ArrayView<int,1> data = elements.view<int,1>(field);
            array::ArrayView<int,1> old_data ( fs.field(fname) );
            for( size_t e=0; e<nb_elems; ++e ) old_data(e) = data(e);
          }
        }
        if( field.rank() == 2 )
        {
         if( field.datatype().kind() == DataType::KIND_REAL64 ) {
            array::ArrayView<double,2> data = elements.view<double,2>(field);
            array::ArrayView<double,2> old_data ( fs.field(fname) );
            for( size_t e=0; e<nb_elems; ++e ) {
              for( size_t v=0; v<field.shape(1); ++v ) {
                old_data(e,v) = data(e,v);
              }
            }
          }
          if( field.datatype().kind() == DataType::KIND_REAL32 ) {
            array::ArrayView<float,2> data = elements.view<float,2>(field);
            array::ArrayView<float,2> old_data ( fs.field(fname) );
            for( size_t e=0; e<nb_elems; ++e ) {
              for( size_t v=0; v<field.shape(1); ++v ) {
                old_data(e,v) = data(e,v);
              }
            }
          }
          if( field.datatype().kind() == DataType::KIND_INT64 ) {
            array::ArrayView<long,2> data = elements.view<long,2>(field);
            array::ArrayView<long,2> old_data ( fs.field(fname) );
            for( size_t e=0; e<nb_elems; ++e ) {
              for( size_t v=0; v<field.shape(1); ++v ) {
                old_data(e,v) = data(e,v);
              }
            }
          }
          if( field.datatype().kind() == DataType::KIND_INT32 ) {
            array::ArrayView<int,2> data = elements.view<int,2>(field);
            array::ArrayView<int,2> old_data ( fs.field(fname) );
            for( size_t e=0; e<nb_elems; ++e ) {
              for( size_t v=0; v<field.shape(1); ++v ) {
                old_data(e,v) = data(e,v);
              }
            }
          }
        }
      }

      if( elements.node_connectivity().rows() )
      {
        const mesh::Elements::Connectivity &node_connectivity = elements.node_connectivity();
        if( !fs.has_field("nodes") )
        {
          fs.create_field<int>("nodes",node_connectivity.cols());
        }
        array::IndexView<int,2> old_nodes ( fs.field("nodes") );
        for( size_t jelem=0; jelem<node_connectivity.rows(); ++jelem )
        for( size_t jnode=0; jnode<node_connectivity.cols(); ++jnode )
          old_nodes(jelem,jnode) = node_connectivity(jelem,jnode);
      }
      if( elements.edge_connectivity().rows() )
      {
        const mesh::Elements::Connectivity &edge_connectivity = elements.edge_connectivity();
        if( !fs.has_field("to_edge") )
        {
          fs.create_field<int>("to_edge",edge_connectivity.cols());
        }
        array::IndexView<int,2> old_edges ( fs.field("to_edge") );
        for( size_t jelem=0; jelem<edge_connectivity.rows(); ++jelem )
        for( size_t jedge=0; jedge<edge_connectivity.cols(); ++jedge )
          old_edges(jelem,jedge) = edge_connectivity(jelem,jedge);
      }
      if( elements.cell_connectivity().rows() )
      {
        const mesh::Elements::Connectivity &cell_connectivity = elements.cell_connectivity();
        if( !fs.has_field("to_elem") )
        {
          fs.create_field<int>("to_elem",cell_connectivity.cols());
        }
        array::IndexView<int,2> old_cells ( fs.field("to_elem") );
        for( size_t jelem=0; jelem<cell_connectivity.rows(); ++jelem )
        for( size_t jcell=0; jcell<cell_connectivity.cols(); ++jcell )
          old_cells(jelem,jcell) = cell_connectivity(jelem,jcell);
      }
    }
  }

  if( nodes().edge_connectivity().rows() )
  {
    const mesh::Nodes::Connectivity& node_edge_connectivity = nodes().edge_connectivity();

    if( ! nodes().has_field("to_edge_size") )
      nodes().add( Field::create<int>( "to_edge_size", array::make_shape(nodes().size(),1) ) );
    array::ArrayView<int,1> to_edge_size ( nodes().field("to_edge_size") );

    // Get max_edge_cnt
    int max_edge_cnt(0);
    for( size_t jnode=0; jnode<nodes().size(); ++jnode )
    {
      to_edge_size(jnode) = node_edge_connectivity.cols(jnode);
      max_edge_cnt = std::max(max_edge_cnt,to_edge_size(jnode));
    }

    ECKIT_MPI_CHECK_RESULT( MPI_Allreduce( MPI_IN_PLACE, &max_edge_cnt, 1, MPI_INT, MPI_MAX, eckit::mpi::comm() ) );

    if( ! nodes().has_field("to_edge") )
      nodes().add( Field::create<int>("to_edge",array::make_shape(nodes().size(),max_edge_cnt)));
    array::IndexView<int,2> node_to_edge ( nodes().field("to_edge") );

    for( size_t jnode=0; jnode<nodes().size(); ++jnode )
    {
      for( size_t jedge=0; jedge<node_edge_connectivity.cols(jnode); ++jedge )
      {
        node_to_edge(jnode,jedge) = node_edge_connectivity(jnode,jedge);
      }
    }
  }


}
#endif

//----------------------------------------------------------------------------------------------------------------------

// C wrapper interfaces to C++ routines

Mesh* atlas__Mesh__new () {
	return new Mesh();
}

void atlas__Mesh__delete (Mesh* This) {
	delete This;
}

#if ! DEPRECATE_OLD_FUNCTIONSPACE
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
      This->nodes().resize(vshape[0]);
    else
      This->create_function_space(std::string(name), std::string(shape_func),vshape);
  );
}

deprecated::FunctionSpace* atlas__Mesh__function_space (Mesh* This, char* name) {
  ATLAS_ERROR_HANDLING(
    ASSERT( This != NULL );
    if( std::string(name) == "nodes" )
      throw eckit::BadParameter("nodes is no longer a *old* FunctionSpace",Here());
    else
      return &This->function_space( std::string(name) );
  );
  return NULL;
}
#endif

mesh::Nodes* atlas__Mesh__create_nodes (Mesh* This, int nb_nodes)
{
  ATLAS_ERROR_HANDLING(
    ASSERT( This );
    This->nodes().resize(nb_nodes);
    return &This->nodes();
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

mesh::Edges* atlas__Mesh__edges (Mesh* This) {
  ATLAS_ERROR_HANDLING(
    ASSERT( This != NULL );
    return &This->edges();
  );
  return NULL;
}

mesh::Cells* atlas__Mesh__cells (Mesh* This) {
  ATLAS_ERROR_HANDLING(
    ASSERT( This != NULL );
    return &This->cells();
  );
  return NULL;
}


//----------------------------------------------------------------------------------------------------------------------

} // namespace mesh
} // namespace atlas
