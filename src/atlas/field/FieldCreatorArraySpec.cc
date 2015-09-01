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
#include "atlas/field/FieldCreatorArraySpec.h"
#include "atlas/util/DataType.h"

namespace atlas {
namespace field {

Field* FieldCreatorArraySpec::create_field( const eckit::Parametrisation& params ) const
{
  std::vector<long> shape;
  if( !params.get("shape",shape) )
    throw eckit::Exception("Could not find parameter 'shape' in Parametrisation");
  std::vector<size_t> s(shape.begin(),shape.end());

  DataType datatype = DataType::create<double>();
  std::string datatype_str;
  if( params.get("datatype", datatype_str) )
  {
    datatype = DataType(datatype_str);
  }
  else
  {
    DataType::kind_t kind(DataType::kind<double>());
    params.get("kind",kind);
    if( ! DataType::kind_valid(kind) )
    {
      std::stringstream msg;
      msg << "Could not create field. kind parameter unrecognized";
      throw eckit::Exception(msg.str());
    }
    datatype = DataType(kind);
  }

  std::string name;
  params.get("name",name);
  return Field::create(name,datatype,s);
}

namespace {
static FieldCreatorBuilder< FieldCreatorArraySpec > __ArraySpec("ArraySpec");
}

// ------------------------------------------------------------------

} // namespace fieldcreator
} // namespace atlas

