#include "Field.hpp"
#include "FunctionSpace.hpp"
#include <typeinfo>       // std::bad_cast
#include <iostream>
namespace ecmwf {

Field::Field(const std::string& name, FunctionSpace& function_space) :
    name_(name), function_space_(function_space)
{
}

template <>
std::vector< int >& Field::data<int>()
{
  try {
    return dynamic_cast< FieldT<int>& >(*this).data();
  } 
  catch (std::bad_cast& bc) {
    std::cerr << "bad_cast caught: " << bc.what() << '\n';
    std::cerr << "Could not cast Field " << name() 
              << " with data_type " << data_type() << " to int32" << std::endl;
    throw bc;
  }
}

template <>
std::vector< float >& Field::data<float>()
{
  try {
    return dynamic_cast< FieldT<float>& >(*this).data();
  } 
  catch (std::bad_cast& bc) {
    std::cerr << "bad_cast caught: " << bc.what() << '\n';
    std::cerr << "Could not cast Field " << name() 
              << " with data_type " << data_type() << " to real32" << std::endl;
    throw bc;
  }
}

template <>
std::vector< double >& Field::data<double>()
{
  try {
    return dynamic_cast< FieldT<double>& >(*this).data();
  } 
  catch (std::bad_cast& bc) {
    std::cerr << "bad_cast caught: " << bc.what() << '\n';
    std::cerr << "Could not cast Field " << name() 
              << " with data_type " << data_type() << " to real64" << std::endl;
    throw bc;
  }
}

template<>
std::string FieldT<int>::data_type() const { return "int32"; }
template<>
std::string FieldT<float>::data_type() const { return "real32"; }
template<>
std::string FieldT<double>::data_type() const { return "real64"; }


// ------------------------------------------------------------------
// C wrapper interfaces to C++ routines

const char* ecmwf__Field__name (Field* This)
{
  return This->name().c_str();
}

const char* ecmwf__Field__data_type (Field* This)
{
  return This->data_type().c_str();
}

Metadata* ecmwf__Field__metadata (Field* This)
{
  return &This->metadata();
}

void ecmwf__Field__data_double (Field* This, double* &field_data, int* &field_bounds, int &rank)
{
  field_data = &This->data<double>()[0];
  field_bounds = const_cast<int*>(&(This->bounds()[0]));
  rank = This->bounds().size();
}

void ecmwf__Field__data_float (Field* This, float* &field_data, int* &field_bounds, int &rank)
{
  field_data = &This->data<float>()[0];
  field_bounds = const_cast<int*>(&(This->bounds()[0]));
  rank = This->bounds().size();
}

void ecmwf__Field__data_int (Field* This, int* &field_data, int* &field_bounds, int &rank)
{
  field_data = &This->data<int>()[0];
  field_bounds = const_cast<int*>(&(This->bounds()[0]));
  rank = This->bounds().size();
}

// ------------------------------------------------------------------

} // namespace ecmwf

