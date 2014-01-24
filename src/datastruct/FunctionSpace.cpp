#include "FunctionSpace.hpp"
#include "Field.hpp"
#include <iostream>
namespace ecmwf {

FunctionSpace::FunctionSpace(const std::string& name, const std::string& shape_func, const std::vector<int>& bounds) :
  name_(name), bounds_(bounds)
{ 
  //std::cout << "C++ : FunctionSpace Constructor" << std::endl;
  dof_ = 1;
  size_t bsize = bounds_.size();
  for (size_t i=0; i<bsize; ++i)
  {
    if( bounds_[i] != Field::NB_VARS )
      dof_ *= bounds[i];
  }
}

FunctionSpace::~FunctionSpace() 
{ 
  //std::cout << "C++ : FunctionSpace Destructor ("<<name_<<")" << std::endl;
  index_.clear();
  for( size_t f=0; f<fields_.size(); ++f )
    if( fields_[f] ) delete(fields_[f]);
  fields_.clear();
}

template <>
void FunctionSpace::create_field<double>(const std::string& name, size_t nb_vars)
{
  //std::cout << "C++ : Create field " << name << " with size " << size*nb_nodes_ << std::endl;
  index_[name] = fields_.size();
  fields_.push_back( new FieldT<double>(name,nb_vars,*this) );

  size_t bsize = bounds_.size();
  std::vector< int > bounds(bsize);
  //std::cout << "Allocating field<real64> " << name << " ( ";

  for (size_t i=0; i<bsize; ++i)
  {
    if( bounds_[i] == Field::NB_VARS )
      bounds[i] = nb_vars;
    else
      bounds[i] = bounds_[i];
    //std::cout << bounds[i];
    //if (i<bsize-1) std::cout << " , ";
  }
  //std::cout << " )" << std::endl;
  fields_.back()->allocate(bounds);
}

template <>
void FunctionSpace::create_field<float>(const std::string& name, size_t nb_vars)
{
  //std::cout << "C++ : Create field " << name << " with size " << size*nb_nodes_ << std::endl;
  index_[name] = fields_.size();
  fields_.push_back( new FieldT<float>(name,nb_vars,*this) );

  size_t bsize = bounds_.size();
  std::vector< int > bounds(bsize);
  //std::cout << "Allocating field<real32> " << name << " ( ";

  for (size_t i=0; i<bsize; ++i)
  {
    if( bounds_[i] == Field::NB_VARS )
      bounds[i] = nb_vars;
    else
      bounds[i] = bounds_[i];
    //std::cout << bounds[i];
    //if (i<bsize-1) std::cout << " , ";
  }
  //std::cout << " )" << std::endl;
  fields_.back()->allocate(bounds);
}

template <>
void FunctionSpace::create_field<int>(const std::string& name, size_t nb_vars)
{
  //std::cout << "C++ : Create field " << name << " with size " << size*nb_nodes_ << std::endl;
  index_[name] = fields_.size();
  fields_.push_back( new FieldT<int>(name,nb_vars,*this) );

  size_t bsize = bounds_.size();
  std::vector< int > bounds(bsize);
  //std::cout << "Allocating field<int32> " << name << " ( ";

  for (size_t i=0; i<bsize; ++i)
  {
    if( bounds_[i] == Field::NB_VARS )
      bounds[i] = nb_vars;
    else
      bounds[i] = bounds_[i];
    //std::cout << bounds[i];
    //if (i<bsize-1) std::cout << " , ";
  }
  //std::cout << " )" << std::endl;
  fields_.back()->allocate(bounds);
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

void FunctionSpace::parallelise(const int proc[], const int glb_idx[])
{
  halo_exchange_.setup(proc,glb_idx,bounds_,bounds_.size()-1);
}

// ------------------------------------------------------------------
// C wrapper interfaces to C++ routines

FunctionSpace* ecmwf__FunctionSpace__new (char* name, char* shape_func, int bounds[], int bounds_size) { 
  std::vector<int> bounds_vec(bounds,bounds+bounds_size);
  return new FunctionSpace( std::string(name), std::string(shape_func), bounds_vec ); 
}

void ecmwf__FunctionSpace__create_field_double (FunctionSpace* This, char* name, int nb_vars) {
  This->create_field<double>( std::string(name), nb_vars );
}

void ecmwf__FunctionSpace__create_field_float (FunctionSpace* This, char* name, int nb_vars) {
  This->create_field<float>( std::string(name), nb_vars );
}

void ecmwf__FunctionSpace__create_field_int (FunctionSpace* This, char* name, int nb_vars) {
  This->create_field<int>( std::string(name), nb_vars );
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

void ecmwf__FunctionSpace__parallelise (FunctionSpace* This, int proc[], int glb_idx[]) {
  This->parallelise(proc,glb_idx);
}

void ecmwf__FunctionSpace__halo_exchange_int (FunctionSpace* This, int field_data[], int field_size) {
  This->halo_exchange(field_data,field_size);
}

void ecmwf__FunctionSpace__halo_exchange_float (FunctionSpace* This, float field_data[], int field_size) {
  This->halo_exchange(field_data,field_size);
}

void ecmwf__FunctionSpace__halo_exchange_double (FunctionSpace* This, double field_data[], int field_size) {
  This->halo_exchange(field_data,field_size);
}

HaloExchange const* ecmwf__FunctionSpace__halo_exchange (FunctionSpace* This) {
  return &This->halo_exchange();
}


void ecmwf__FunctionSpace__delete (FunctionSpace* This)  {
  delete This;
}
// ------------------------------------------------------------------

} // namespace ecmwf

