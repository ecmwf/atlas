/*
 * (C) Copyright 1996-2014 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include <typeinfo> // std::bad_cast
#include <sstream>
#include <stdexcept>

#include "eckit/exception/Exceptions.h"
#include "eckit/types/Types.h"

#include "atlas/runtime/ErrorHandling.h"
#include "atlas/FunctionSpace.h"
#include "atlas/Field.h"
#include "atlas/field/FieldT.h"
#include "atlas/field/FieldCreator.h"
#include "atlas/field/FieldTCreator.h"

#include "atlas/Grid.h"
#include "atlas/FunctionSpace.h"
#include "atlas/Mesh.h"

using atlas::field::FieldT;

namespace atlas {

template<class T>
inline std::ostream &operator<<(std::ostream &s, const std::vector<T> &v) {
    return eckit::__print_list(s, v);
}

Field* Field::create(const eckit::Parametrisation& params)
{
  std::string creator_factory;
  if( params.get("creator",creator_factory) )
  {
    eckit::ScopedPtr<field::FieldCreator> creator
        (field::FieldCreatorFactory::build(creator_factory,params) );
    return creator->create_field(params);
  }
  else
    throw eckit::Exception("Could not find parameter 'creator' in Parametrisation for call to Field::create()");

  return NULL;
}

Field* Field::create(const ArrayShape& shape, const eckit::Parametrisation& params)
{
  std::string data_type = "real64";
  params.get("data_type",data_type);

  eckit::ScopedPtr<field::FieldTCreator> creator
     (field::FieldTCreatorFactory::build("FieldT<"+data_type+">") );
  return creator->create_field(shape,params);
}

Field::Field(const std::string& name, const size_t nb_vars) :
  name_(name), nb_vars_(nb_vars), grid_(0), function_space_(0)
{
}

Field::Field(const eckit::Parametrisation& params) :
  name_(), nb_vars_(1), grid_(0), function_space_(0)
{
  Grid::Id grid;
  if( params.get("grid",grid) )
    grid_ = &Grid::from_id(grid);

  FunctionSpace::Id function_space;
  if( params.get("function_space",function_space) )
    function_space_ = &FunctionSpace::from_id(function_space);

  params.get("name",name_);
}

Field::~Field()
{
}

namespace {
template< typename DATA_TYPE >
DATA_TYPE* get_field_data( const Field& field )
{
  const FieldT<DATA_TYPE>* fieldT = dynamic_cast< const FieldT<DATA_TYPE>* >(&field);
  if( fieldT == NULL )
  {
    std::stringstream msg;
    msg << "Could not cast Field " << field.name()
        << " with data_type " << field.data_type() << " to "
        << data_type_to_str<DATA_TYPE>();
    throw eckit::BadCast(msg.str(),Here());
  }
  return const_cast<DATA_TYPE*>(fieldT->data());
}
}

template <> const int*    Field::data<int   >() const { return get_field_data<int   >(*this); }
template <>       int*    Field::data<int   >()       { return get_field_data<int   >(*this); }
template <> const long*   Field::data<long  >() const { return get_field_data<long  >(*this); }
template <>       long*   Field::data<long  >()       { return get_field_data<long  >(*this); }
template <> const float*  Field::data<float >() const { return get_field_data<float >(*this); }
template <>       float*  Field::data<float >()       { return get_field_data<float >(*this); }
template <> const double* Field::data<double>() const { return get_field_data<double>(*this); }
template <>       double* Field::data<double>()       { return get_field_data<double>(*this); }

void Field::print(std::ostream& os) const
{
    os << "Field[name=" << name()
       << ",datatype=" << data_type_
       << ",size=" << size()
       << ",shape=" << shape_
       << ",strides=" << strides_
       << "]";
}

std::ostream& operator<<( std::ostream& os, const Field& f)
{
  f.print(os);
  return os;
}



// ===========================================
// TO BE REMOVED


FunctionSpace& Field::function_space() const {
    if( !function_space_ )
        throw eckit::Exception("Field "+name()+" is not associated to any FunctionSpace.");
    return *function_space_;
}

const Grid& Field::grid() const {
    if( function_space_ )
    {
        if( function_space_->mesh().has_grid() )
            return function_space_->mesh().grid();
    }
    if( !grid_ )
        throw eckit::Exception("Field "+name()+" is not associated to any Grid.");
    return *grid_;
}

const Mesh& Field::mesh() const { ASSERT(function_space_); return function_space_->mesh(); }

Mesh& Field::mesh() { ASSERT(function_space_); return function_space_->mesh(); }

void Field::set_function_space(const FunctionSpace& function_space)
{
  function_space_ = const_cast<FunctionSpace*>(&function_space);
}


// END TO BE REMOVED
// ============================================


// ------------------------------------------------------------------
// C wrapper interfaces to C++ routines

Field* atlas__Field__create(eckit::Parametrisation* params)
{
  Field* field = 0;
  ATLAS_ERROR_HANDLING( field = Field::create(*params) );
  ASSERT(field);
  return field;
}

void atlas__Field__delete (Field* This)
{
  delete This;
}

const char* atlas__Field__name (Field* This)
{
  ATLAS_ERROR_HANDLING(
    ASSERT( This != NULL );
    return This->name().c_str();
  );
  return NULL;
}

const char* atlas__Field__data_type (Field* This)
{
  return This->data_type().c_str();
}

int atlas__Field__size (Field* This)
{
  ATLAS_ERROR_HANDLING(
    ASSERT( This != NULL );
    return This->size();
  );
  return 0;
}

double atlas__Field__bytes (Field* This)
{
  ATLAS_ERROR_HANDLING(
    ASSERT( This != NULL );
    return This->bytes();
  );
  return 0;
}

int atlas__Field__nb_vars (Field* This)
{
  return This->nb_vars();
}

Metadata* atlas__Field__metadata (Field* This)
{
  return &This->metadata();
}

FunctionSpace* atlas__Field__function_space (Field* This)
{
  return &This->function_space();
}

void atlas__Field__shapef (Field* This, int* &shape, int &rank)
{
  shape = const_cast<int*>(&This->shapef().front());
  rank = This->shapef().size();
}

void atlas__Field__data_shapef_int (Field* This, int* &field_data, int* &field_bounds, int &rank)
{
  field_data = &This->data<int>()[0];
  field_bounds = const_cast<int*>(&(This->shapef()[0]));
  rank = This->shapef().size();
}

void atlas__Field__data_shapef_long (Field* This, long* &field_data, int* &field_bounds, int &rank)
{
  field_data = &This->data<long>()[0];
  field_bounds = const_cast<int*>(&(This->shapef()[0]));
  rank = This->shapef().size();
}

void atlas__Field__data_shapef_float (Field* This, float* &field_data, int* &field_bounds, int &rank)
{
  field_data = &This->data<float>()[0];
  field_bounds = const_cast<int*>(&(This->shapef()[0]));
  rank = This->shapef().size();
}

void atlas__Field__data_shapef_double (Field* This, double* &field_data, int* &field_bounds, int &rank)
{
  field_data = &This->data<double>()[0];
  field_bounds = const_cast<int*>(&(This->shapef()[0]));
  rank = This->shapef().size();
}

// ------------------------------------------------------------------

} // namespace atlas

