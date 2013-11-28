#include "Mesh.hpp"
#include "FunctionSpace.hpp"
namespace ecmwf {

Mesh::~Mesh() 
{ 
  index_.clear();
  for( size_t f=0; f<function_spaces_.size(); ++f )
    if( function_spaces_[f] ) delete(function_spaces_[f]);
  function_spaces_.clear();
}

void Mesh::add_function_space( FunctionSpace* function_space )
{
  index_[function_space->name()] = function_spaces_.size();
  function_spaces_.push_back( function_space );
}

FunctionSpace& Mesh::function_space(const std::string& name)
{
  return *function_spaces_[ index_.at(name) ];
}

// ------------------------------------------------------------------
// C wrapper interfaces to C++ routines

Mesh* ecmwf__Mesh__new () { 
  return new Mesh(); 
}

void ecmwf__Mesh__delete (Mesh* This) {
  delete This;
}

void ecmwf__Mesh__add_function_space (Mesh* This, FunctionSpace* function_space) {
  This->add_function_space(function_space);
}

FunctionSpace* ecmwf__Mesh__function_space (Mesh* This, char* name) {
  return &This->function_space( std::string(name) );
}
// ------------------------------------------------------------------


} // namespace ecmwf

