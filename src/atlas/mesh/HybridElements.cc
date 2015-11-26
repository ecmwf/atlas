/*
 * (C) Copyright 1996-2015 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "atlas/runtime/ErrorHandling.h"
#include "atlas/mesh/ElementType.h"
#include "atlas/mesh/HybridElements.h"
#include "atlas/mesh/Elements.h"
#include "atlas/Field.h"
#include "atlas/atlas_config.h"
#include "atlas/atlas_defines.h"

#ifdef ATLAS_HAVE_FORTRAN
#define FORTRAN_BASE 1
#define TO_FORTRAN +1
#else
#define FORTRAN_BASE 0
#endif

namespace atlas {
namespace mesh {

//------------------------------------------------------------------------------------------------------

HybridElements::HybridElements() :
  size_(0),
  elements_size_(),
  elements_begin_(1,0ul),
  type_idx_()
{
  global_index_ = &add( Field::create<gidx_t>("glb_idx",   make_shape(size(),1)) );
  remote_index_ = &add( Field::create<int   >("remote_idx",make_shape(size(),1)) );
  partition_    = &add( Field::create<int   >("partition", make_shape(size(),1)) );
  ghost_        = &add( Field::create<int   >("ghost",     make_shape(size(),1)) );
  
  node_connectivity_ = &add( "node", new Connectivity() );  
  edge_connectivity_ = &add( "edge", new Connectivity() );  
  cell_connectivity_ = &add( "cell", new Connectivity() );  
}

HybridElements::~HybridElements()
{
}

Field& HybridElements::add( Field* field )
{
  ASSERT( field != NULL );
  ASSERT( ! field->name().empty() );

  if( has_field(field->name()) ) {
    std::stringstream msg;
    msg << "Trying to add field '"<<field->name()<<"' to Nodes, but Nodes already has a field with this name.";
    throw eckit::Exception(msg.str(),Here());
  }
  fields_[field->name()] = eckit::SharedPtr<Field>(field);
  return *field;
}

void HybridElements::resize( size_t size )
{
  size_ = size;
  for( FieldMap::iterator it = fields_.begin(); it != fields_.end(); ++it )
  {
    Field& field = *it->second;
    ArrayShape shape = field.shape();
    shape[0] = size_;
    field.resize(shape);
  }
}

void HybridElements::remove_field(const std::string& name)
{
  if( ! has_field(name) )
  {
    std::stringstream msg;
    msg << "Trying to remove field `"<<name<<"' in Nodes, but no field with this name is present in Nodes.";
    throw eckit::Exception(msg.str(),Here());
  }
  fields_.erase(name);
}


const Field& HybridElements::field(const std::string& name) const
{
  if( ! has_field(name) )
  {
    std::stringstream msg;
    msg << "Trying to access field `"<<name<<"' in Nodes, but no field with this name is present in Nodes.";
    throw eckit::Exception(msg.str(),Here());
  }
  return *fields_.find(name)->second;
}

Field& HybridElements::field(const std::string& name)
{
  return const_cast<Field&>(static_cast<const HybridElements*>(this)->field(name));
}

const Field& HybridElements::field(size_t idx) const
{
  ASSERT(idx < nb_fields());
  size_t c(0);
  for( FieldMap::const_iterator it = fields_.begin(); it != fields_.end(); ++it )
  {
    if( idx == c )
    {
      const Field& field = *it->second;
      return field;
    }
    c++;
  }
  eckit::SeriousBug("Should not be here!",Here());
  static Field* ret;
  return *ret;
}

Field& HybridElements::field(size_t idx)
{
  return const_cast<Field&>(static_cast<const HybridElements*>(this)->field(idx));
}


HybridElements::Connectivity& HybridElements::add( const std::string& name, Connectivity *connectivity )
{
  connectivities_[name] = connectivity;
  return *connectivity;
}


size_t HybridElements::add( const ElementType* element_type, size_t nb_elements, const std::vector<idx_t> &connectivity )
{
  return add(element_type,nb_elements,connectivity.data());
}

size_t HybridElements::add( const ElementType* element_type, size_t nb_elements, const idx_t connectivity[] )
{
  return add(element_type,nb_elements,connectivity,false);
}

size_t HybridElements::add( const ElementType* element_type, size_t nb_elements, const idx_t connectivity[], bool fortran_array )
{
  eckit::SharedPtr<const ElementType> etype ( element_type );

  size_t old_size=size();
  size_t new_size = old_size+nb_elements;

  size_t nb_nodes = etype->nb_nodes();
  size_t nb_edges = etype->nb_edges();

  type_idx_.resize(new_size);

  for( size_t e=old_size; e<new_size; ++e ) {
    type_idx_[e] = element_types_.size();
  }

  elements_begin_.push_back(new_size);
  elements_size_.push_back(nb_elements);

  element_types_.push_back( etype );
  elements_.resize(element_types_.size());
  for( size_t t=0; t<nb_types(); ++t ) {
    elements_[t].reset( new Elements(*this,t) );
  }
  
  node_connectivity_->add(nb_elements,nb_nodes,connectivity);
  resize( new_size );
  return element_types_.size()-1;
}

size_t HybridElements::add( const Elements& elems )
{
  bool fortran_array = true;
  return add( &elems.element_type(), elems.size(), elems.node_connectivity().data(), fortran_array );
}


const std::string& HybridElements::name(size_t elem_idx) const { return element_types_[type_idx_[elem_idx]]->name(); }

//-----------------------------------------------------------------------------

extern "C" {

}

//-----------------------------------------------------------------------------

}  // namespace mesh
}  // namespace atlas

#undef FORTRAN_BASE


#include "atlas/Mesh.h"
#include "atlas/FunctionSpace.h"
#include "atlas/util/IndexView.h"
#include "atlas/util/Debug.h"
#include "atlas/Parameters.h"

namespace atlas {
namespace mesh {
namespace temporary {

HybridElements* Convert::createCells( const Mesh& mesh )
{
  return createElements(mesh,Entity::ELEMS);
}
HybridElements* Convert::createFaces( const Mesh& mesh )
{
  return createElements(mesh,Entity::FACES);
}
HybridElements* Convert::createElements( const Mesh& mesh, Entity::Type type )
{
  HybridElements* hybrid_elements = new HybridElements();
  
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
          
          eckit::Log::info() << fname << "<"<< field.datatype().str() << ">" << std::endl;
          
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
    }
  }
  return hybrid_elements;
}


}
}
}



