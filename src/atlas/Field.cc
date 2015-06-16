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
#include "eckit/thread/AutoLock.h"
#include "eckit/thread/Mutex.h"

#include "atlas/Field.h"
#include "atlas/FunctionSpace.h"


namespace {

    static eckit::Mutex *local_mutex = 0;
    static std::map<std::string, atlas::FieldCreatorFactory *> *m = 0;
    static pthread_once_t once = PTHREAD_ONCE_INIT;

    static void init() {
        local_mutex = new eckit::Mutex();
        m = new std::map<std::string, atlas::FieldCreatorFactory *>();
    }
}


namespace atlas {


class ArraySpec: public FieldCreator
{
public:
  ArraySpec() {}
  ArraySpec(const eckit::Parametrisation&) {}
  virtual Field* create_field( const eckit::Parametrisation& ) const;
};

Field* ArraySpec::create_field( const eckit::Parametrisation& params ) const
{
  Field* field;

  std::vector<long> shape;
  std::string data_type;
  if( !params.get("shape",shape) )
    throw eckit::Exception("Could not find parameter 'shape' in Parametrisation");
  if( !params.get("data_type",data_type) )
    throw eckit::Exception("Could not find parameter 'data_type' in Parametrisation");
  std::vector<size_t> s(shape.begin(),shape.end());
  if( data_type == "int32" || data_type == "int" )
    field = new FieldT<int>(s);
  else if( data_type == "int64" || data_type == "long" )
    field = new FieldT<long>(s);
  else if( data_type == "real32" || data_type == "float" )
    field = new FieldT<float>(s);
  else if( data_type == "real64" || data_type == "double" )
    field = new FieldT<double>(s);
  else
    throw eckit::Exception("Could not create field. data_type parameter unrecognized: "+data_type);

  filter_params(*field,params);

  return field;
}

namespace {
static FieldCreatorBuilder< ArraySpec > __ArraySpec("ArraySpec");
}




Field::Field(const std::string& name, const size_t nb_vars) :
  name_(name), nb_vars_(nb_vars), function_space_(NULL), grid_(NULL)
{
}

Field::~Field()
{
}

Field* Field::create(const eckit::Parametrisation& params)
{
  std::string creator_factory;
  if( params.get("creator",creator_factory) )
  {
    eckit::ScopedPtr<FieldCreator> creator( FieldCreatorFactory::build(creator_factory,params) );
    return creator->create_field(params);
  }
  else
    throw eckit::Exception("Could not find parameter 'creator' in Parametrisation for call to Field::create()");

  return NULL;
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
        << " with data_type " << field.data_type() << " "
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


template<>
void FieldT<int>::halo_exchange() { function_space().halo_exchange(data_.data(),data_.size()); }
template<>
void FieldT<long>::halo_exchange() { function_space().halo_exchange(data_.data(),data_.size()); }
template<>
void FieldT<float>::halo_exchange() { function_space().halo_exchange(data_.data(),data_.size()); }
template<>
void FieldT<double>::halo_exchange() { function_space().halo_exchange(data_.data(),data_.size()); }

std::ostream& operator<<( std::ostream& os, const Field& f)
{
  f.print(os);
  return os;
}

// ------------------------------------------------------------------

FieldCreator::FieldCreator()
{
}

FieldCreator::~FieldCreator()
{
}

void FieldCreator::filter_params(Field& field, const eckit::Parametrisation& params) const
{
  Grid::Id grid;
  if( params.get("grid",grid) )
    field.set_grid(Grid::from_id(grid));

  FunctionSpace::Id function_space;
  if( params.get("function_space",function_space) )
    field.set_function_space(FunctionSpace::from_id(function_space));

  std::string name;
  if( params.get("name",name) )
    field.set_name(name);
}

FieldCreatorFactory::FieldCreatorFactory(const std::string &name):
    name_(name) {

    pthread_once(&once, init);

    eckit::AutoLock<eckit::Mutex> lock(local_mutex);

    ASSERT(m->find(name) == m->end());
    (*m)[name] = this;
}


FieldCreatorFactory::~FieldCreatorFactory() {
    eckit::AutoLock<eckit::Mutex> lock(local_mutex);
    m->erase(name_);
}


void FieldCreatorFactory::list(std::ostream& out) {
    pthread_once(&once, init);

    eckit::AutoLock<eckit::Mutex> lock(local_mutex);

    const char* sep = "";
    for (std::map<std::string, FieldCreatorFactory *>::const_iterator j = m->begin() ; j != m->end() ; ++j) {
        out << sep << (*j).first;
        sep = ", ";
    }
}


FieldCreator *FieldCreatorFactory::build(const std::string &name) {

    pthread_once(&once, init);

    eckit::AutoLock<eckit::Mutex> lock(local_mutex);
    std::map<std::string, FieldCreatorFactory *>::const_iterator j = m->find(name);

    eckit::Log::debug() << "Looking for FieldCreatorFactory [" << name << "]" << std::endl;

    if (j == m->end()) {
        eckit::Log::error() << "No FieldCreatorFactory for [" << name << "]" << std::endl;
        eckit::Log::error() << "FieldCreatorFactories are:" << std::endl;
        for (j = m->begin() ; j != m->end() ; ++j)
            eckit::Log::error() << "   " << (*j).first << std::endl;
        throw eckit::SeriousBug(std::string("No FieldCreatorFactory called ") + name);
    }

    return (*j).second->make();
}

FieldCreator *FieldCreatorFactory::build(const std::string& name, const eckit::Parametrisation& param) {

    pthread_once(&once, init);

    eckit::AutoLock<eckit::Mutex> lock(local_mutex);

    std::map<std::string, FieldCreatorFactory *>::const_iterator j = m->find(name);

    eckit::Log::debug() << "Looking for FieldCreatorFactory [" << name << "]" << std::endl;

    if (j == m->end()) {
        eckit::Log::error() << "No FieldCreatorFactory for [" << name << "]" << std::endl;
        eckit::Log::error() << "FieldCreatorFactories are:" << std::endl;
        for (j = m->begin() ; j != m->end() ; ++j)
            eckit::Log::error() << "   " << (*j).first << std::endl;
        throw eckit::SeriousBug(std::string("No FieldCreatorFactory called ") + name);
    }

    return (*j).second->make(param);
}



// ------------------------------------------------------------------
// C wrapper interfaces to C++ routines

const char* atlas__Field__name (Field* This)
{
  return This->name().c_str();
}

const char* atlas__Field__data_type (Field* This)
{
  return This->data_type().c_str();
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

