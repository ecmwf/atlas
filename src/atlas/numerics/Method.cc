/*
 * (C) Copyright 1996-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */


#include "eckit/exception/Exceptions.h"
#include "atlas/atlas_config.h"
#include "atlas/runtime/ErrorHandling.h"
#include "atlas/numerics/Method.h"

namespace atlas {
namespace numerics {
namespace {

void assert_shared(const eckit::Owned* owned)
{
  if( owned->owners() == 0 )
  {
    throw eckit::SeriousBug("Cannot create shared_ptr from stack allocated or scoped_ptr",Here());
  }
}

}

eckit::SharedPtr<Method const> Method::shared_from_this() const
{
  assert_shared(this);
  return eckit::SharedPtr<Method const>(this);
}

eckit::SharedPtr<Method> Method::shared_from_this()
{
  assert_shared(this);
  return eckit::SharedPtr<Method>(this);
}

eckit::SharedPtr<Method const> Method::ptr() const
{
  assert_shared(this);
  return eckit::SharedPtr<Method const>(this);
}

eckit::SharedPtr<Method const> Method::cptr() const
{
  ASSERT(owners()!=0);
  return eckit::SharedPtr<Method const>(this);
}

eckit::SharedPtr<Method> Method::ptr()
{
  ASSERT(owners()!=0);
  return eckit::SharedPtr<Method>(this);
}

//----------------------------------------------------------------------------------------------------------------------

// C wrapper interfaces to C++ routines
extern "C" {
 void atlas__Method__delete (Method* This)
{
  ATLAS_ERROR_HANDLING(
    ASSERT( This );
    delete This;
    This = 0;
  );
}

const char* atlas__Method__name (Method* This) {
  ATLAS_ERROR_HANDLING(
    ASSERT( This );
    return This->name().c_str();
  );
  return 0;
}
}

// ------------------------------------------------------------------

} // namespace numerics
} // namespace atlas

