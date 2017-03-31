/*
 * (C) Copyright 1996-2017 ECMWF.
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
#include <limits>
#include "eckit/exception/Exceptions.h"
#include "eckit/types/Types.h"
#include "atlas/library/config.h"
#include "atlas/mesh/actions/BuildParallelFields.h"
#include "atlas/field/Field.h"
#include "atlas/functionspace/FunctionSpace.h"
#include "atlas/array/DataType.h"
#include "atlas/runtime/ErrorHandling.h"

namespace atlas {
namespace functionspace {

//-----------------------------------------------------------------------------

FunctionSpace::FunctionSpace() :
  functionspace_( new NoFunctionSpace() ) {
}

FunctionSpace::FunctionSpace( const functionspace_t* functionspace ) :
  functionspace_( functionspace ) {
}

FunctionSpace::FunctionSpace( const FunctionSpace& functionspace ) :
  functionspace_( functionspace.functionspace_ ) {
}

std::string FunctionSpace::name() const {
  return functionspace_->name();
}

FunctionSpace::operator bool() const {
  return functionspace_->operator bool();
}

size_t FunctionSpace::footprint() const {
  return functionspace_->footprint();
}

//-----------------------------------------------------------------------------

// C wrapper interfaces to C++ routines
extern "C" {
 void atlas__FunctionSpace__delete (FunctionSpaceImpl* This)
{
  ATLAS_ERROR_HANDLING(
    ASSERT( This );
    delete This;
    This = 0;
  );
}

const char* atlas__FunctionSpace__name (FunctionSpaceImpl* This) {
  ATLAS_ERROR_HANDLING(
    ASSERT( This );
    return This->name().c_str();
  );
  return 0;
}
}

// ------------------------------------------------------------------

} // namespace functionspace
} // namespace atlas

