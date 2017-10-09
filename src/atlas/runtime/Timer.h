/*
 * (C) Copyright 1996-2017 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#pragma once

#include "atlas/library/config.h"
#include "atlas/runtime/timer/TimerT.h"
#include "atlas/runtime/timer/TimerBarriers.h"
#include "atlas/runtime/timer/TimerLogging.h"
#include "atlas/runtime/timer/Timings.h"
#include "atlas/runtime/timer/TimerNesting.h"

//-----------------------------------------------------------------------------------------------------------

/// Create scoped timer objects
///
/// Example:
///
///     void foo() {
///         ATLAS_TIME();
///         // timer "foo" starts
///
///         /* interesting computations ... */
///
///         ATLAS_TIME_SCOPE("bar") {
///             // timer "bar" starts
///
///             /* interesting computations ... */
///
///             // timer "bar" ends
///         }
///
///         // timer "foo" ends
///     }
///
/// Example 2:
///
///     void foo() {
///         ATLAS_TIME("custom");
///         // timer "custom" starts
///
///         /* interesting computations ... */
///
///         // timer "custom" ends
///     }
///
#define ATLAS_TIME(...)
#define ATLAS_TIME_SCOPE(...)

//-----------------------------------------------------------------------------------------------------------

namespace atlas {

struct TimerTraits {
    using Barriers = runtime::timer::TimerNoBarriers;
    using Logging  = runtime::timer::TimerLogging;
    using Timings  = runtime::timer::Timings;
    using Nesting  = runtime::timer::TimerNesting;
};

class Timer : public runtime::timer::TimerT< TimerTraits > {
    using Base = runtime::timer::TimerT< TimerTraits >;
public:
    using Base::Base;
};

} // namespace atlas


//-----------------------------------------------------------------------------------------------------------

#if ATLAS_HAVE_TIMINGS

#include "atlas/util/detail/BlackMagic.h"

#undef ATLAS_TIME
#undef ATLAS_TIME_SCOPE
#define ATLAS_TIME(...)       ATLAS_TIME_( __ATLAS__NARG(__VA_ARGS__), ##__VA_ARGS__ )
#define ATLAS_TIME_SCOPE(...) ATLAS_TIME_SCOPE_( __ATLAS__NARG(__VA_ARGS__), ##__VA_ARGS__ )

#define ATLAS_TIME_SCOPE_(N, ...) __ATLAS__SPLICE( ATLAS_TIME_SCOPE_, N)(__VA_ARGS__)
#define ATLAS_TIME_SCOPE_0() \
  for( ::atlas::Timer __ATLAS__SPLICE( timer, __LINE__ ) (Here());\
    __ATLAS__SPLICE( timer, __LINE__ ) .running(); \
    __ATLAS__SPLICE( timer, __LINE__ ) .stop() )
#define ATLAS_TIME_SCOPE_1(arg1) \
  for( ::atlas::Timer __ATLAS__SPLICE( timer, __LINE__ ) (Here(),arg1);\
    __ATLAS__SPLICE( timer, __LINE__ ) .running(); \
    __ATLAS__SPLICE( timer, __LINE__ ) .stop() )
#define ATLAS_TIME_SCOPE_2(arg1,arg2) \
  for( ::atlas::Timer __ATLAS__SPLICE( timer, __LINE__ ) (Here(),arg1,arg2);\
    __ATLAS__SPLICE( timer, __LINE__ ) .running(); \
    __ATLAS__SPLICE( timer, __LINE__ ) .stop() )

#define ATLAS_TIME_(N, ...) __ATLAS__SPLICE( ATLAS_TIME_, N)(__VA_ARGS__)
#define ATLAS_TIME_0()           ::atlas::Timer __ATLAS__SPLICE( timer, __LINE__ ) (Here());
#define ATLAS_TIME_1(arg1)       ::atlas::Timer __ATLAS__SPLICE( timer, __LINE__ ) (Here(),arg1);
#define ATLAS_TIME_2(arg1,arg2)  ::atlas::Timer __ATLAS__SPLICE( timer, __LINE__ ) (Here(),arg1,arg2);

#endif

//-----------------------------------------------------------------------------------------------------------
