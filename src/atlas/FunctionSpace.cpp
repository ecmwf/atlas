// (C) Copyright 1996-2014 ECMWF.

#include <cassert>
#include <iostream>
#include <sstream>

#include "atlas/FunctionSpace.hpp"
#include "atlas/Field.hpp"

namespace atlas {

FunctionSpace::FunctionSpace(const std::string& name, const std::string& shape_func, const std::vector<int>& bounds) :
  name_(name), bounds_(bounds)
{ 
  //std::cout << "C++ : FunctionSpace Constructor" << std::endl;
  dof_ = 1;
  size_t bsize = bounds_.size();
  for (size_t i=0; i<bsize; ++i)
  {
    if( bounds_[i] != Field::UNDEF_VARS )
      dof_ *= bounds_[i];
  }
  glb_dof_ = dof_;
}

FunctionSpace::~FunctionSpace() 
{ 
//  std::cout << "FunctionSpace Destructor ("<<name_<<")" << std::endl;
  index_.clear();
  for( size_t f=0; f<fields_.size(); ++f )
    if( fields_[f] ) delete(fields_[f]);
  fields_.clear();
}

void FunctionSpace::resize(const std::vector<int>& bounds)
{
  if (bounds.size() != bounds_.size() )
    throw std::runtime_error("Cannot resize functionspace: Bounds sizes don't match.");

  size_t bsize = bounds_.size();
  for (size_t i=0; i<bsize-1; ++i)
  {
    if (bounds[i] != bounds_[i])
      throw std::runtime_error("Only the last bound can be resized for now!");
  }

  bounds_ = bounds;

  dof_ = 1;
  for (size_t i=0; i<bsize; ++i)
  {
    if( bounds_[i] != Field::UNDEF_VARS )
      dof_ *= bounds_[i];
  }

  for( int f=0; f<fields_.size(); ++f)
  {
    size_t bsize = bounds_.size();
    std::vector< int > field_bounds(bsize);
    for (size_t i=0; i<bsize; ++i)
    {
      if( bounds_[i] == Field::UNDEF_VARS )
        field_bounds[i] = fields_[f]->nb_vars();
      else
        field_bounds[i] = bounds_[i];
    }
    fields_[f]->allocate(field_bounds);
  }
}

template <>
FieldT<double>& FunctionSpace::create_field(const std::string& name, size_t nb_vars)
{
    if( has_field(name) )
    {
        std::ostringstream msg; msg << "field with name " << name << "already exists" << std::endl;
        throw std::runtime_error( msg.str() );
    }

  //std::cout << "C++ : Create field " << name << " with size " << size*nb_nodes_ << std::endl;
  index_[name] = fields_.size();
  FieldT<double>* field = new FieldT<double>(name,nb_vars,*this);
  fields_.push_back( field );

  size_t bsize = bounds_.size();
  std::vector< int > bounds(bsize);
  //std::cout << "Allocating field<real64> " << name << " ( ";

  for( size_t i = 0; i < bsize; ++i )
  {
    if( bounds_[i] == Field::UNDEF_VARS )
      bounds[i] = nb_vars;
    else
      bounds[i] = bounds_[i];
    //std::cout << bounds[i];
    //if (i<bsize-1) std::cout << " , ";
  }
  //std::cout << " )" << std::endl;
  field->allocate(bounds);
  return *field;
}

template <>
FieldT<float>& FunctionSpace::create_field(const std::string& name, size_t nb_vars)
{
  //std::cout << "C++ : Create field " << name << " with size " << size*nb_nodes_ << std::endl;
  index_[name] = fields_.size();
  FieldT<float>* field = new FieldT<float>(name,nb_vars,*this);
  fields_.push_back( field );

  size_t bsize = bounds_.size();
  std::vector< int > bounds(bsize);
  //std::cout << "Allocating field<real32> " << name << " ( ";

  for (size_t i=0; i<bsize; ++i)
  {
    if( bounds_[i] == Field::UNDEF_VARS )
      bounds[i] = nb_vars;
    else
      bounds[i] = bounds_[i];
    //std::cout << bounds[i];
    //if (i<bsize-1) std::cout << " , ";
  }
  //std::cout << " )" << std::endl;
  field->allocate(bounds);
  return *field;
}

template <>
FieldT<int>& FunctionSpace::create_field(const std::string& name, size_t nb_vars)
{
  //std::cout << "C++ : Create field " << name << " with size " << size*nb_nodes_ << std::endl;
  index_[name] = fields_.size();
  FieldT<int>* field = new FieldT<int>(name,nb_vars,*this);
  fields_.push_back( field );

  size_t bsize = bounds_.size();
  std::vector< int > bounds(bsize);
  //std::cout << "Allocating field<int32> " << name << " ( ";

  for (size_t i=0; i<bsize; ++i)
  {
    if( bounds_[i] == Field::UNDEF_VARS )
      bounds[i] = nb_vars;
    else
      bounds[i] = bounds_[i];
    //std::cout << bounds[i];
    //if (i<bsize-1) std::cout << " , ";
  }
  //std::cout << " )" << std::endl;
  field->allocate(bounds);
  return *field;
}

void FunctionSpace::remove_field(const std::string& name)
{
  //std::cout << "C++ : Create field " << name << " with size " << size*nb_nodes_ << std::endl;
  delete( fields_[ index_.at(name) ] );
  fields_[ index_.at(name) ] = 0;
  index_.erase(name);
}

Field& FunctionSpace::field( size_t idx )
{
    assert( idx < fields_.size() );
    return *fields_[ idx ];
}

Field& FunctionSpace::field(const std::string& name)
{
    //std::cout << "C++ : Access field " << name << std::endl;
    try {
      return *fields_[ index_.at(name) ];
    }
    catch( std::out_of_range& e ) {
      std::stringstream msg;
      msg << "Could not find field \"" << name << "\" in FunctionSpace \"" << name_ << "\"";
      throw std::out_of_range(msg.str());
    }
}

template<>
  FieldT<double> &FunctionSpace::field(const std::string& name)
{
  try {
    return *dynamic_cast< FieldT<double>* >(fields_[ index_.at(name) ]);
  }
  catch( std::out_of_range& e ) {
    std::stringstream msg;
    msg << "Could not find field \"" << name << "\" in FunctionSpace \"" << name_ << "\"";
    throw std::out_of_range(msg.str());
  }
}

template<>
  FieldT<float> &FunctionSpace::field(const std::string& name)
{
  try {
    return *dynamic_cast< FieldT<float>* >(fields_[ index_.at(name) ]);
  }
  catch( std::out_of_range& e ) {
    std::stringstream msg;
    msg << "Could not find field \"" << name << "\" in FunctionSpace \"" << name_ << "\"";
    throw std::out_of_range(msg.str());
  }
}

template<>
  FieldT<int> &FunctionSpace::field(const std::string& name)
{
  try {
    return *dynamic_cast< FieldT<int>* >(fields_[ index_.at(name) ]);
  }
  catch( std::out_of_range& e ) {
    std::stringstream msg;
    msg << "Could not find field \"" << name << "\" in FunctionSpace \"" << name_ << "\"";
    throw std::out_of_range(msg.str());
  }
}

void FunctionSpace::parallelise(const int proc[], const int glb_idx[], const int master_glb_idx[] )
{
  halo_exchange_.setup(proc,glb_idx,master_glb_idx,bounds_,bounds_.size()-1);
  gather_.setup(proc,glb_idx,master_glb_idx,bounds_,bounds_.size()-1);
  glb_dof_ = gather_.glb_dof();
  for( int b=bounds_.size()-2; b>=0; --b)
  {
    if( bounds_[b] != Field::UNDEF_VARS )
      glb_dof_ *= bounds_[b];
  }
}

void FunctionSpace::parallelise()
{
  FieldT<int>& proc = field<int>("proc");
  FieldT<int>& glb_idx = field<int>("glb_idx");
  FieldT<int>& master_glb_idx = field<int>("master_glb_idx");
  parallelise(proc.data().data(),glb_idx.data().data(),master_glb_idx.data().data());
}

// ------------------------------------------------------------------
// C wrapper interfaces to C++ routines

FunctionSpace* atlas__FunctionSpace__new (char* name, char* shape_func, int bounds[], int bounds_size) { 
  std::vector<int> bounds_vec(bounds,bounds+bounds_size);
  return new FunctionSpace( std::string(name), std::string(shape_func), bounds_vec ); 
}

int atlas__FunctionSpace__dof (FunctionSpace* This) {
  return This->dof();
}

int atlas__FunctionSpace__glb_dof (FunctionSpace* This) {
  return This->glb_dof();
}

void atlas__FunctionSpace__create_field_double (FunctionSpace* This, char* name, int nb_vars) {
  This->create_field<double>( std::string(name), nb_vars );
}

void atlas__FunctionSpace__create_field_float (FunctionSpace* This, char* name, int nb_vars) {
  This->create_field<float>( std::string(name), nb_vars );
}

void atlas__FunctionSpace__create_field_int (FunctionSpace* This, char* name, int nb_vars) {
  This->create_field<int>( std::string(name), nb_vars );
}

void atlas__FunctionSpace__remove_field (FunctionSpace* This, char* name ) {
  This->remove_field( std::string(name) );
}

int atlas__FunctionSpace__has_field (FunctionSpace* This, char* name) {
  return This->has_field( std::string(name) );
}

const char* atlas__FunctionSpace__name (FunctionSpace* This) {
  return This->name().c_str();
}

void atlas__FunctionSpace__bounds (FunctionSpace* This, int* &bounds, int &rank) {
  bounds = const_cast<int*>(&(This->bounds()[0]));
  rank = This->bounds().size();
}

Field* atlas__FunctionSpace__field (FunctionSpace* This, char* name) {
  return &This->field( std::string(name) );
}

void atlas__FunctionSpace__parallelise (FunctionSpace* This, int proc[], int glb_idx[], int master_glb_idx[]) {
  This->parallelise(proc,glb_idx,master_glb_idx);
}

void atlas__FunctionSpace__halo_exchange_int (FunctionSpace* This, int field_data[], int field_size) {
  This->halo_exchange(field_data,field_size);
}

void atlas__FunctionSpace__halo_exchange_float (FunctionSpace* This, float field_data[], int field_size) {
  This->halo_exchange(field_data,field_size);
}

void atlas__FunctionSpace__halo_exchange_double (FunctionSpace* This, double field_data[], int field_size) {
  This->halo_exchange(field_data,field_size);
}

void atlas__FunctionSpace__gather_int (FunctionSpace* This, int field_data[], int field_size, int glbfield_data[], int glbfield_size) {
  This->gather(field_data,field_size, glbfield_data,glbfield_size);
}

void atlas__FunctionSpace__gather_float (FunctionSpace* This, float field_data[], int field_size, float glbfield_data[], int glbfield_size) {
  This->gather(field_data,field_size, glbfield_data,glbfield_size);
}

void atlas__FunctionSpace__gather_double (FunctionSpace* This, double field_data[], int field_size, double glbfield_data[], int glbfield_size) {
  This->gather(field_data,field_size, glbfield_data,glbfield_size);
}

HaloExchange const* atlas__FunctionSpace__halo_exchange (FunctionSpace* This) {
  return &This->halo_exchange();
}


void atlas__FunctionSpace__delete (FunctionSpace* This)  {
  delete This;
}
// ------------------------------------------------------------------

} // namespace atlas

