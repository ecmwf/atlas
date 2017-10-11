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
#include "atlas/runtime/timer/TimerTracing.h"

//-----------------------------------------------------------------------------------------------------------

/// Create scoped timer objects
///
/// Example:
///
///     void foo() {
///         ATLAS_TRACE();
///         // trace "foo" starts
///
///         /* interesting computations ... */
///
///         ATLAS_TRACE_SCOPE("bar") {
///             // trace "bar" starts
///
///             /* interesting computations ... */
///
///             // trace "bar" ends
///         }
///
///         // trace "foo" ends
///     }
///
/// Example 2:
///
///     void foo() {
///         ATLAS_TRACE("custom");
///         // trace "custom" starts
///
///         /* interesting computations ... */
///
///         // trace "custom" ends
///     }
///
#define ATLAS_TRACE(...)
#define ATLAS_TRACE_SCOPE(...)

//-----------------------------------------------------------------------------------------------------------

namespace atlas {

struct TraceTraits {
    using Barriers = runtime::timer::TimerBarriersNone;
    using Tracing  = runtime::timer::TimerTracing;
};

class Trace : public runtime::timer::TimerT< TraceTraits > {
    using Base = runtime::timer::TimerT< TraceTraits >;
public:
    using Base::Base;
};

} // namespace atlas


//-----------------------------------------------------------------------------------------------------------

#if ATLAS_HAVE_TRACE

#include "atlas/util/detail/BlackMagic.h"

#undef ATLAS_TRACE
#undef ATLAS_TRACE_SCOPE

#define ATLAS_TRACE(...)       __ATLAS_TYPE(       ::atlas::Trace, Here(), ##__VA_ARGS__ )
#define ATLAS_TRACE_SCOPE(...) __ATLAS_TYPE_SCOPE( ::atlas::Trace, Here(), ##__VA_ARGS__ )

#endif

//-----------------------------------------------------------------------------------------------------------
