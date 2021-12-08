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

#if ATLAS_HAVE_GRIDTOOLS_STORAGE
#include "atlas/array/gridtools/GridToolsArrayView.h"
#else
#include "atlas/array/native/NativeArrayView.h"
#endif

namespace atlas {
namespace array {

#define EXPLICIT_TEMPLATE_DECLARATION_TYPE_RANK(TYPE, RANK) \
    extern template class ArrayView<TYPE, RANK>;            \
    extern template class ArrayView<const TYPE, RANK>;

#define EXPLICIT_TEMPLATE_DECLARATION(RANK)               \
    EXPLICIT_TEMPLATE_DECLARATION_TYPE_RANK(int, RANK);   \
    EXPLICIT_TEMPLATE_DECLARATION_TYPE_RANK(long, RANK);  \
    EXPLICIT_TEMPLATE_DECLARATION_TYPE_RANK(float, RANK); \
    EXPLICIT_TEMPLATE_DECLARATION_TYPE_RANK(double, RANK);

// The Fujitsu compiler complains about missing references with the below
// so we just skip it in that case
#if !defined(__FUJITSU)

// For each RANK in [1..9]
EXPLICIT_TEMPLATE_DECLARATION(1)
EXPLICIT_TEMPLATE_DECLARATION(2)
EXPLICIT_TEMPLATE_DECLARATION(3)
EXPLICIT_TEMPLATE_DECLARATION(4)
EXPLICIT_TEMPLATE_DECLARATION(5)
EXPLICIT_TEMPLATE_DECLARATION(6)
EXPLICIT_TEMPLATE_DECLARATION(7)
EXPLICIT_TEMPLATE_DECLARATION(8)
EXPLICIT_TEMPLATE_DECLARATION(9)

#endif

#undef EXPLICIT_TEMPLATE_DECLARATION_TYPE_RANK
#undef EXPLICIT_TEMPLATE_DECLARATION

}  // namespace array
}  // namespace atlas
