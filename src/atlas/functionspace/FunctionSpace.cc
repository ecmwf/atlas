/*
 * (C) Copyright 1996-2016 ECMWF.
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
#include "atlas/atlas_config.h"
#include "atlas/mesh/actions/BuildParallelFields.h"
#include "atlas/field/Field.h"
#include "atlas/functionspace/FunctionSpace.h"
#include "atlas/internals/Bitflags.h"
#include "atlas/array/DataType.h"
#include "atlas/util/runtime/ErrorHandling.h"

namespace atlas {
namespace functionspace {
namespace {

void assert_shared(const eckit::Owned* owned)
{
  if( owned->owners() == 0 )
  {
    throw eckit::SeriousBug("Cannot create shared_ptr from stack allocated or scoped_ptr",Here());
  }
}

}

eckit::SharedPtr<FunctionSpace const> FunctionSpace::shared_from_this() const
{
  assert_shared(this);
  return eckit::SharedPtr<FunctionSpace const>(this);
}

eckit::SharedPtr<FunctionSpace> FunctionSpace::shared_from_this()
{
  assert_shared(this);
  return eckit::SharedPtr<FunctionSpace>(this);
}

eckit::SharedPtr<FunctionSpace const> FunctionSpace::ptr() const
{
  assert_shared(this);
  return eckit::SharedPtr<FunctionSpace const>(this);
}

eckit::SharedPtr<FunctionSpace const> FunctionSpace::cptr() const
{
  ASSERT(owners()!=0);
  return eckit::SharedPtr<FunctionSpace const>(this);
}

eckit::SharedPtr<FunctionSpace> FunctionSpace::ptr()
{
  ASSERT(owners()!=0);
  return eckit::SharedPtr<FunctionSpace>(this);
}


//-----------------------------------------------------------------------------

// C wrapper interfaces to C++ routines
extern "C" {
 void atlas__FunctionSpace__delete (FunctionSpace* This)
{
  ATLAS_ERROR_HANDLING(
    ASSERT( This );
    delete This;
    This = 0;
  );
}

const char* atlas__FunctionSpace__name (FunctionSpace* This) {
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

