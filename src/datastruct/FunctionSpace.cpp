#include "FunctionSpace.hpp"
#include <iostream>
#include "Field.hpp"

namespace ecmwf {

FunctionSpace::FunctionSpace(const std::string& name, const std::string& shape_func, size_t nb_nodes) :
  name_(name), bounds_(1,nb_nodes)
{ 
  //std::cout << "C++ : FunctionSpace Constructor" << std::endl;
}

FunctionSpace::FunctionSpace(const std::string& name, const std::string& shape_func, const std::vector<size_t>& bounds) :
  name_(name), bounds_(bounds)
{ 
  //std::cout << "C++ : FunctionSpace Constructor" << std::endl;
}

FunctionSpace::~FunctionSpace() 
{ 
  std::cout << "C++ : FunctionSpace Destructor ("<<name_<<")" << std::endl;
  index_.clear();
  for( size_t f=0; f<fields_.size(); ++f )
    if( fields_[f] ) delete(fields_[f]);
  fields_.clear();
}

void FunctionSpace::create_field(const std::string& name, size_t nb_cols)
{
  //std::cout << "C++ : Create field " << name << " with size " << size*nb_nodes_ << std::endl;
  index_[name] = fields_.size();
  fields_.push_back( new Field(name,nb_cols,*this) );
  fields_.back()->allocate();
}

void FunctionSpace::remove_field(const std::string& name)
{
  //std::cout << "C++ : Create field " << name << " with size " << size*nb_nodes_ << std::endl;
  delete( fields_[ index_.at(name) ] );
  fields_[ index_.at(name) ] = 0;
  index_.erase(name);
}

Field& FunctionSpace::field(const std::string& name)
{
  //std::cout << "C++ : Access field " << name << std::endl;
  return *fields_[ index_.at(name) ];
}

void FunctionSpace::print_field(const std::string& name)
{
  Field& field_obj = field( name );
  //std::cout << "C++ : " << name << " : ";
  for( size_t i=0; i<field_obj.size(); ++i)
    std::cout << field_obj.data()[i] << "  ";
  std::cout << std::endl;
}

PrismaticFunctionSpace::PrismaticFunctionSpace(const std::string& name, const std::string& shape_func, const std::vector<size_t>& bounds) :
  FunctionSpace(name, shape_func, bounds),
  nb_levels_(bounds[0]),
  nb_nodes_(bounds[1])
{
}

PrismaticFunctionSpace::~PrismaticFunctionSpace()
{ 
}

void PrismaticFunctionSpace::create_field(const std::string& name, size_t nb_cols)
{
  index_[name] = fields_.size();
  fields_.push_back( new PrismaticField(name,nb_cols,*this) );
  fields_.back()->allocate();
  std::cout << "C++ : Created prismatic field " << name << " with size " << fields_.back()->size() << std::endl;
}

PrismaticField& PrismaticFunctionSpace::field(const std::string& name)
{
  //std::cout << "C++ : Access field " << name << std::endl;
  return *dynamic_cast<PrismaticField*>(fields_[ index_.at(name) ]);
}
// ------------------------------------------------------------------
// C wrapper interfaces to C++ routines

FunctionSpace* ecmwf__FunctionSpace__new (char* name, char* shape_func, int nb_nodes) { 
  return new FunctionSpace( std::string(name), std::string(shape_func), nb_nodes ); 
}

PrismaticFunctionSpace* ecmwf__PrismaticFunctionSpace__new (char* name, char* shape_func, int nb_levels, int nb_nodes) { 
  std::cout << "creating new prismatic function space" << std::endl;
  std::vector<size_t> bounds(2);
  bounds[0] = nb_levels;
  bounds[1] = nb_nodes;
  return new PrismaticFunctionSpace( std::string(name), std::string(shape_func), bounds ); 
}

void ecmwf__FunctionSpace__create_field (FunctionSpace* This, char* name, int size) {
  This->create_field( std::string(name), size );
}

void ecmwf__FunctionSpace__remove_field (FunctionSpace* This, char* name ) {
  This->remove_field( std::string(name) );
}

const char* ecmwf__FunctionSpace__name (FunctionSpace* This) {
  return This->name().c_str();
}

Field* ecmwf__FunctionSpace__field (FunctionSpace* This, char* name) {
  return &This->field( std::string(name) );
} 

void ecmwf__FunctionSpace__field_data (FunctionSpace* This, char* name, double* &field_data, int &field_size) { 
  Field& field = This->field( std::string(name) );
  field_data = &field.data()[0];
  field_size = field.size();
}

void ecmwf__FunctionSpace__print_field (FunctionSpace* This, char* name) {
  This->print_field( std::string(name) );
}

void ecmwf__FunctionSpace__delete (FunctionSpace* This)  {
  delete This;
}
// ------------------------------------------------------------------

} // namespace ecmwf

