/*
 * (C) Copyright 1996-2015 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "atlas/runtime/ErrorHandling.h"
#include "atlas/Parameters.h"
#include "atlas/mesh/ElementType.h"

namespace atlas {
namespace mesh {

//------------------------------------------------------------------------------------------------------

static ElementType* create( const std::string& )
{
  return 0;
}

ElementType::ElementType() {}
ElementType::~ElementType() {}

//-----------------------------------------------------------------------------

extern "C" {

}

}  // namespace mesh
}  // namespace atlas

