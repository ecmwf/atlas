/*
 * (C) Copyright 1996-2014 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include <sstream>
#include "eckit/exception/Exceptions.h"
#include "eckit/config/Parametrisation.h"
#include "atlas/Field.h"
#include "atlas/field/ArraySpec.h"

namespace atlas {
namespace field {

Field* ArraySpec::create_field( const eckit::Parametrisation& params ) const
{
  Field* field;

  std::vector<long> shape;
  std::string data_type = "real64";
  if( !params.get("shape",shape) )
    throw eckit::Exception("Could not find parameter 'shape' in Parametrisation");
  
  
  std::vector<size_t> s;
  bool fortran=false;
  params.get("fortran",fortran);
  if( fortran ) {
    s.resize(shape.size());
    size_t n=shape.size();
    for( size_t j=0; j<shape.size(); ++j ) s[j] = shape[--n];
  }
  else
    s.assign(shape.begin(),shape.end());
  params.get("data_type",data_type);

  std::string name;
  params.get("name",name);

  std::ostream& out = eckit::Log::debug();
  out << "Creating ArraySpec "<<data_type<<" field: "<<name;
  for( size_t j=0; j<s.size(); ++j ) out << "["<<s[j]<<"]";
  out << '\n';

  if( data_type == "int32" || data_type == "int" )
    field = new FieldT<int>(s,params);
  else if( data_type == "int64" || data_type == "long" )
    field = new FieldT<long>(s,params);
  else if( data_type == "real32" || data_type == "float" )
    field = new FieldT<float>(s,params);
  else if( data_type == "real64" || data_type == "double" )
    field = new FieldT<double>(s,params);
  else
    throw eckit::Exception("Could not create field. data_type parameter unrecognized: "+data_type);

  return field;
}

namespace {
static FieldCreatorBuilder< ArraySpec > __ArraySpec("ArraySpec");
}

// ------------------------------------------------------------------

} // namespace fieldcreator
} // namespace atlas

