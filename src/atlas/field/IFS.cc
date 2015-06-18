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
#include "atlas/util/ArrayUtil.h"
#include "atlas/field/IFS.h"
#include "atlas/Field.h"
#include "atlas/Grid.h"

namespace atlas {
namespace field {

Field* IFS::create_field( const eckit::Parametrisation& params ) const
{
  Field* field;
  size_t ngptot;
  size_t nblk;
  size_t nvar = 1;
  size_t nproma = 1;
  size_t nlev = 1;
  long kind = 8; // 8 bytes = double

  if( !params.get("ngptot",ngptot) )
  {
    Grid::Id grid;
    if( params.get("grid",grid ) )
      ngptot = Grid::from_id(grid).npts();
    else
      throw eckit::Exception("Could not find parameter 'ngptot' in Parametrisation");
  }
  params.get("nproma",nproma);
  params.get("nlev",nlev);
  params.get("nvar",nvar);
  params.get("kind",kind);

  std::string data_type;
  if( !params.get("data_type", data_type) )
  {
    // Assume real
    params.get("kind",kind);
    if( kind == 4 ) data_type = "real32";
    else if( kind == 8 ) data_type = "real64";
    else {
      std::stringstream msg;
      msg << "Could not create field. kind parameter unrecognized (expected: 4 or 8; received: " << kind << ") ";
      throw eckit::Exception(msg.str());
    }
  }

  nblk = std::ceil(static_cast<double>(ngptot)/static_cast<double>(nproma));
  ArrayShape s = make_shape(nblk,nvar,nlev,nproma);

  std::string name;
  params.get("name",name);
  eckit::Log::debug() << "Creating IFS "<<data_type<<" field: "<<name<<"[nblk="<<nblk<<"][nvar="<<nvar<<"][nlev="<<nlev<<"][nproma="<<nproma<<"]\n";

  if( data_type == "int32" || data_type == "int" )
    field = new FieldT<int>(s,params);
  else if( data_type == "int64" || data_type == "long" )
    field = new FieldT<long>(s,params);
  else if( data_type == "real32" || data_type == "float" )
    field = new FieldT<float>(s,params);
  else if( data_type == "real64" || data_type == "double" )
    field = new FieldT<double>(s);

  return field;
}

namespace {
static FieldCreatorBuilder< IFS > __IFS("IFS");
}

// ------------------------------------------------------------------

} // namespace field
} // namespace atlas

