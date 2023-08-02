/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#pragma once

#include "atlas/library/config.h"
#include "atlas/runtime/trace/Barriers.h"
#include "atlas/runtime/trace/Logging.h"
#include "atlas/runtime/trace/TraceT.h"

//-----------------------------------------------------------------------------------------------------------

/// Create scoped trace objects
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
#define ATLAS_TRACE_BARRIERS(enabled)

//-----------------------------------------------------------------------------------------------------------

namespace atlas {

namespace runtime {
namespace trace {
struct TraceTraits {
#if ATLAS_HAVE_TRACE_BARRIERS
    using Barriers = ::atlas::runtime::trace::Barriers;
#else
    using Barriers = ::atlas::runtime::trace::NoBarriers;
#endif
    using Tracing = ::atlas::runtime::trace::Logging;
};

}  // namespace trace
}  // namespace runtime

class Trace : public runtime::trace::TraceT<runtime::trace::TraceTraits> {
    using Base   = runtime::trace::TraceT<runtime::trace::TraceTraits>;
    using Traits = runtime::trace::TraceTraits;

public:
    using Base::Base;
};

}  // namespace atlas

//-----------------------------------------------------------------------------------------------------------

#if ATLAS_HAVE_TRACE

#include "atlas/library/detail/BlackMagic.h"

#undef ATLAS_TRACE
#undef ATLAS_TRACE_SCOPE
#undef ATLAS_TRACE_BARRIERS

#define ATLAS_TRACE(...) __ATLAS_TYPE(::atlas::Trace, Here() __ATLAS_COMMA_ARGS(__VA_ARGS__))
#define ATLAS_TRACE_SCOPE(...) __ATLAS_TYPE_SCOPE(::atlas::Trace, Here() __ATLAS_COMMA_ARGS(__VA_ARGS__))
#define ATLAS_TRACE_BARRIERS(enabled) __ATLAS_TYPE(::atlas::Trace::Barriers, enabled)

#endif

//-----------------------------------------------------------------------------------------------------------
