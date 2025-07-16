/*
 * (C) Copyright 2024- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#pragma once

#include "pluto/pluto_config.h"

#pragma push_macro("STD_MDSPAN_NAMESPACE")
#undef STD_MDSPAN_NAMESPACE

#if defined(__clang__)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunknown-warning-option"
#pragma clang diagnostic ignored "-Wpre-c++2b-compat"
#endif

#if PLUTO_HAVE_MDSPAN
#include <mdspan>
#define STD_MDSPAN_NAMESPACE std
#else
#pragma push_macro("MDSPAN_IMPL_STANDARD_NAMESPACE")
#undef MDSPAN_IMPL_STANDARD_NAMESPACE
#define MDSPAN_IMPL_STANDARD_NAMESPACE pluto::detail::mdspan
#define MDSPAN_USE_PAREN_OPERATOR 0
    // In std::mdspan the call operator doesn't exist, so we don't include it here.
#include "pluto/detail/mdspan/mdspan.hpp"
#define STD_MDSPAN_NAMESPACE MDSPAN_IMPL_STANDARD_NAMESPACE
#endif

#ifdef __clang__
#pragma clang diagnostic pop
#endif

#if defined(__INTEL_COMPILER)
#pragma warning pop
#endif


// ------------------------------------------------------------------------------------------------

namespace pluto {
using ::STD_MDSPAN_NAMESPACE::dynamic_extent;
using ::STD_MDSPAN_NAMESPACE::layout_left;
using ::STD_MDSPAN_NAMESPACE::layout_right;
using ::STD_MDSPAN_NAMESPACE::layout_stride;
using ::STD_MDSPAN_NAMESPACE::default_accessor;
using ::STD_MDSPAN_NAMESPACE::extents;
using ::STD_MDSPAN_NAMESPACE::dextents;
using ::STD_MDSPAN_NAMESPACE::mdspan;
// A C++26 addition:
template< std::size_t Rank, class IndexType = std::size_t >
using dims = dextents<IndexType, Rank>;
} // namespace pluto

// ------------------------------------------------------------------------------------------------

#if !PLUTO_HAVE_MDSPAN
#pragma pop_macro("MDSPAN_IMPL_STANDARD_NAMESPACE")
#endif
#pragma pop_macro("STD_MDSPAN_NAMESPACE")
