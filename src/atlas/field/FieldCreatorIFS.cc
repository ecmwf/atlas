/*
 * (C) Copyright 1996-2015 ECMWF.
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
#include "atlas/util/ArrayUtil.h"
#include "atlas/util/DataType.h"
#include "atlas/field/FieldCreatorIFS.h"
#include "atlas/Field.h"
#include "atlas/Grid.h"

namespace atlas {
namespace field {

Field* FieldCreatorIFS::create_field( const eckit::Parametrisation& params ) const
{
  size_t ngptot;
  size_t nblk;
  size_t nvar = 1;
  size_t nproma = 1;
  size_t nlev = 1;
  long kind = DataType::kind<double>();

  if( !params.get("ngptot",ngptot) )
    throw eckit::Exception("Could not find parameter 'ngptot' in Parametrisation");
  params.get("nproma",nproma);
  params.get("nlev",nlev);
  params.get("nvar",nvar);

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

  nblk = std::ceil(static_cast<double>(ngptot)/static_cast<double>(nproma));

  ArrayShape s;
  bool fortran (false);
    params.get("fortran",fortran);
  if( fortran ) s = make_shape(nproma,nlev,nvar,nblk);
  else          s = make_shape(nblk,nvar,nlev,nproma);

  std::string name;
  params.get("name",name);
  eckit::Log::debug() << "Creating IFS "<<datatype.str()<<" field: "<<name<<"[nblk="<<nblk<<"][nvar="<<nvar<<"][nlev="<<nlev<<"][nproma="<<nproma<<"]\n";

  return Field::create(name,datatype,s);

}

namespace {
static FieldCreatorBuilder< FieldCreatorIFS > __IFS("IFS");
}

// ------------------------------------------------------------------

} // namespace field
} // namespace atlas

