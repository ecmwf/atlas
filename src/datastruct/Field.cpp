#include "Field.hpp"
#include "FunctionSpace.hpp"
#include <iostream>

namespace ecmwf {

Field::Field(const std::string& name, size_t nb_cols, FunctionSpace& function_space) :
  name_(name), nb_cols_(nb_cols), data_(0), function_space_(function_space)
{
}

Field::~Field()
{
  std::cout << "C++ : Field Destructor (" << name_ << ")" << std::endl;

  name_ = "";
  data_.clear();
}

void Field::allocate()
{
  bounds_.clear();
  bounds_.reserve(function_space_.bounds().size()+1);
  size_t size(1);
  for( int i=0; i<function_space_.bounds().size(); ++i ) 
  {
    const size_t& bound = function_space_.bounds()[i];
    size *= bound;
    bounds_.push_back(bound);
  }
  size *= nb_cols_;
  bounds_.push_back(nb_cols_);
  data_.resize(size);
}

PrismaticField::PrismaticField(const std::string& name, size_t nb_cols, FunctionSpace& function_space)
  : Field(name, nb_cols, function_space)
{
  PrismaticFunctionSpace& fs = this->function_space();
  nb_levels_ = fs.nb_levels();
  nb_nodes_  = fs.nb_nodes();
}

PrismaticField::~PrismaticField()
{
}

PrismaticFunctionSpace& PrismaticField::function_space()
{ 
  return *dynamic_cast<PrismaticFunctionSpace*>(&function_space_);
}


// ------------------------------------------------------------------
// C wrapper interfaces to C++ routines

const char* ecmwf__Field__name (Field* This)
{
  return This->name().c_str();
}

Metadata* ecmwf__Field__metadata (Field* This)
{
  return &This->metadata();
}

void ecmwf__Field__data (Field* This, double* &field_data, int* &field_bounds, int &rank)
{
  field_data = &This->data()[0];
  field_bounds = const_cast<int*>(&(This->bounds()[0]));
  rank = This->bounds().size();
}

// ------------------------------------------------------------------

} // namespace ecmwf

