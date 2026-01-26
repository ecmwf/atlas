/*
 * (C) Copyright 2025- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#pragma once

#include "pluto/pluto_config.h"

#if defined(__clang__)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunknown-warning-option"
#pragma clang diagnostic ignored "-Wpre-c++2b-compat"
#endif

#if PLUTO_HAVE_MDSPAN
#include <mdspan>

#pragma push_macro("STD_MDSPAN_NAMESPACE")
#undef STD_MDSPAN_NAMESPACE
#define STD_MDSPAN_NAMESPACE std
namespace pluto {
using ::STD_MDSPAN_NAMESPACE::dynamic_extent;
using ::STD_MDSPAN_NAMESPACE::layout_left;
using ::STD_MDSPAN_NAMESPACE::layout_right;
using ::STD_MDSPAN_NAMESPACE::layout_stride;
using ::STD_MDSPAN_NAMESPACE::default_accessor;
using ::STD_MDSPAN_NAMESPACE::extents;
using ::STD_MDSPAN_NAMESPACE::dextents;

#if PLUTO_MDSPAN_USE_PAREN_OPERATOR
#include "pluto/detail/mdspan_paren_operator.h"
#else
using ::STD_MDSPAN_NAMESPACE::mdspan;
#endif
} // namespace pluto
#pragma pop_macro("STD_MDSPAN_NAMESPACE")
#define PLUTO_MDSPAN_USE_BRACKET_OPERATOR 1

#else

#pragma push_macro("MDSPAN_IMPL_STANDARD_NAMESPACE")
#undef MDSPAN_IMPL_STANDARD_NAMESPACE
#define MDSPAN_IMPL_STANDARD_NAMESPACE pluto
#include "pluto/detail/mdspan/mdspan.hpp"
#pragma pop_macro("MDSPAN_IMPL_STANDARD_NAMESPACE")
#define PLUTO_MDSPAN_USE_BRACKET_OPERATOR MDSPAN_USE_BRACKET_OPERATOR
#endif

namespace pluto {
    // A C++26 addition:
    template< std::size_t Rank, class IndexType = std::size_t >
    using dims = dextents<IndexType, Rank>;
}

#ifdef __clang__
#pragma clang diagnostic pop
#endif

#if defined(__INTEL_COMPILER)
#pragma warning pop
#endif

// ------------------------------------------------------------------------------------------------
