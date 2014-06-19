/*
 * (C) Copyright 1996-2014 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */



#include <cassert>
#include <iostream>
#include <sstream>

#include "atlas/atlas_defines.h"
#include "atlas/mesh/FunctionSpace.hpp"
#include "atlas/mesh/Field.hpp"
#include "atlas/actions/BuildParallelFields.hpp"
#include "atlas/util/Debug.hpp"

#ifdef HAVE_FORTRAN_NUMBERING
#define REMOTE_IDX_BASE 1
#else
#define REMOTE_IDX_BASE 0
#endif


namespace atlas {

FunctionSpace::FunctionSpace(const std::string& name, const std::string& shape_func, const std::vector<int>& extents) :
  name_(name), extents_(extents)
{ 
  //std::cout << "C++ : FunctionSpace Constructor" << std::endl;
  dof_ = 1;
  size_t extsize = extents_.size();
  bounds_.resize(extsize);
  for (size_t i=0; i<extsize; ++i)
  {
    bounds_[extsize-1-i] = extents_[i];
    if( extents_[i] != Field::UNDEF_VARS )
      dof_ *= extents_[i];
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

void FunctionSpace::resize(const std::vector<int>& extents)
{
  if (extents.size() != extents_.size() )
    throw std::runtime_error("Cannot resize functionspace: extents sizes don't match.");

  size_t extsize = extents_.size();
  
  for (size_t i=1; i<extsize; ++i)
  {
    if (extents[i] != extents_[i])
      throw std::runtime_error("Only the first extent can be resized for now!");
  }

  extents_ = extents;
  bounds_.resize(extsize);
  for (size_t i=0; i<extsize; ++i)
  {
    bounds_[extsize-1-i] = extents_[i];
  }

  dof_ = 1;
  for (size_t i=0; i<extsize; ++i)
  {
    if( extents_[i] != Field::UNDEF_VARS )
      dof_ *= extents_[i];
  }

  for( int f=0; f<fields_.size(); ++f)
  {
    std::vector< int > field_extents(extsize);
    for (size_t i=0; i<extsize; ++i)
    {
      if( extents_[i] == Field::UNDEF_VARS )
        field_extents[i] = fields_[f]->nb_vars();
      else
        field_extents[i] = extents_[i];
    }
    fields_[f]->allocate(field_extents);
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

  // std::cout << "C++ : Create field " << name << " with vars" << nb_vars << std::endl;
  index_[name] = fields_.size();
  FieldT<double>* field = new FieldT<double>(name,nb_vars,*this);
  fields_.push_back( field );

  size_t extsize = extents_.size();
  std::vector< int > field_extents(extsize);
  for (size_t i=0; i<extsize; ++i)
  {
    if( extents_[i] == Field::UNDEF_VARS )
      field_extents[i] = field->nb_vars();
    else
      field_extents[i] = extents_[i];
  }

  field->allocate(field_extents);
  return *field;
}

template <>
FieldT<float>& FunctionSpace::create_field(const std::string& name, size_t nb_vars)
{
  // std::cout << "C++ : Create field " << name << " with size " << size*nb_nodes_ << std::endl;
  index_[name] = fields_.size();
  FieldT<float>* field = new FieldT<float>(name,nb_vars,*this);
  fields_.push_back( field );

  size_t bsize = bounds_.size();
  std::vector< int > bounds(bsize);
  //std::cout << "Allocating field<real32> " << name << " ( ";

  size_t extsize = extents_.size();
  std::vector< int > field_extents(extsize);
  for (size_t i=0; i<extsize; ++i)
  {
    if( extents_[i] == Field::UNDEF_VARS )
      field_extents[i] = field->nb_vars();
    else
      field_extents[i] = extents_[i];
  }

  field->allocate(field_extents);
  return *field;
}

template <>
FieldT<int>& FunctionSpace::create_field(const std::string& name, size_t nb_vars)
{
  // std::cout << "C++ : Create field " << name << " with vars " << nb_vars << std::endl;
  index_[name] = fields_.size();
  FieldT<int>* field = new FieldT<int>(name,nb_vars,*this);
  fields_.push_back( field );

  size_t bsize = bounds_.size();
  std::vector< int > bounds(bsize);
  //std::cout << "Allocating field<int32> " << name << " ( ";

  size_t extsize = extents_.size();
  std::vector< int > field_extents(extsize);
  for (size_t i=0; i<extsize; ++i)
  {
    if( extents_[i] == Field::UNDEF_VARS )
      field_extents[i] = field->nb_vars();
    else
      field_extents[i] = extents_[i];
  }

  field->allocate(field_extents);
  return *field;
}

void FunctionSpace::remove_field(const std::string& name)
{
  //std::cout << "C++ : Create field " << name << " with size " << size*nb_nodes_ << std::endl;
  if( has_field(name) )
  {
    delete( fields_[ index_.at(name) ] );
    fields_[ index_.at(name) ] = 0;
    index_.erase(name);
  }
  else
  {    
    std::stringstream msg;
    msg << "Could not find field \"" << name << "\" in FunctionSpace \"" << name_ << "\"";
    throw std::out_of_range(msg.str());
  }
}

Field& FunctionSpace::field( size_t idx ) const
{
  assert( idx < fields_.size() );
  return *fields_[ idx ];
}

Field& FunctionSpace::field(const std::string& name) const
{
    //std::cout << "C++ : Access field " << name << std::endl;
  if( has_field(name) )
  {
    return *fields_[ index_.at(name) ];
  }
  else
  {
    std::stringstream msg;
    msg << "Could not find field \"" << name << "\" in FunctionSpace \"" << name_ << "\"";
    throw std::out_of_range(msg.str());
  }
}

template<>
  FieldT<double> &FunctionSpace::field(const std::string& name) const
{
  if( has_field(name) )
  {
    return *dynamic_cast< FieldT<double>* >(fields_[ index_.at(name) ]);
  }
  else
  {
    std::stringstream msg;
    msg << "Could not find field \"" << name << "\" in FunctionSpace \"" << name_ << "\"";
    throw std::out_of_range(msg.str());
  }
}

template<>
  FieldT<float> &FunctionSpace::field(const std::string& name) const
{
  if( has_field(name) )
  {
    return *dynamic_cast< FieldT<float>* >(fields_[ index_.at(name) ]);
  }
  else
  {
    std::stringstream msg;
    msg << "Could not find field \"" << name << "\" in FunctionSpace \"" << name_ << "\"";
    throw std::out_of_range(msg.str());
  }
}

template<>
  FieldT<int> &FunctionSpace::field(const std::string& name) const
{
  if( has_field(name) )
  {
    return *dynamic_cast< FieldT<int>* >(fields_[ index_.at(name) ]);
  }
  else
  {
    std::stringstream msg;
    msg << "Could not find field \"" << name << "\" in FunctionSpace \"" << name_ << "\"";
    throw std::out_of_range(msg.str());
  }
}

void FunctionSpace::parallelise(const int part[], const int remote_idx[], int parsize)
{
  halo_exchange_.setup(part,remote_idx,REMOTE_IDX_BASE,parsize);
//  gather_.setup(proc,glb_idx,master_glb_idx,bounds_,bounds_.size()-1);
  glb_dof_ = gather_.glb_dof();
  for( int b=bounds_.size()-2; b>=0; --b)
  {
    if( bounds_[b] != Field::UNDEF_VARS )
      glb_dof_ *= bounds_[b];
  }
}

void FunctionSpace::parallelise()
{
  if( name() == "nodes" )
  {
    FieldT<int>& part       = field<int>("partition");
    FieldT<int>& remote_idx = field<int>("remote_idx");
    parallelise(part.data(),remote_idx.data(),part.size());
  }
  else
  {
    NOTIMP;
  }
}

// ------------------------------------------------------------------
// C wrapper interfaces to C++ routines

FunctionSpace* atlas__FunctionSpace__new (char* name, char* shape_func, int extents[], int extents_size) { 
  std::vector<int> extents_vec(extents,extents+extents_size);
  return new FunctionSpace( std::string(name), std::string(shape_func), extents_vec ); 
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

void atlas__FunctionSpace__boundsf (FunctionSpace* This, int* &bounds, int &rank) {
  bounds = const_cast<int*>(&(This->boundsf()[0]));
  rank = This->boundsf().size();
}

Field* atlas__FunctionSpace__field (FunctionSpace* This, char* name) {
  return &This->field( std::string(name) );
}

void atlas__FunctionSpace__parallelise (FunctionSpace* This) {
  This->parallelise();
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

