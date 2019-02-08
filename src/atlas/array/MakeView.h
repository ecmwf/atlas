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

#include "atlas/array/ArrayView.h"
#include "atlas/array/IndexView.h"
#include "atlas/array/LocalView.h"
#include "atlas/array_fwd.h"

namespace atlas {
namespace array {

extern template IndexView<idx_t, 1> make_indexview<idx_t, 1>( const Array& );

#define EXPLICIT_TEMPLATE_INSTANTIATION( Rank )                                                                        \
    extern template ArrayView<int, Rank, Intent::ReadOnly> make_view<int, Rank, Intent::ReadOnly>( const Array& );     \
    extern template ArrayView<int, Rank, Intent::ReadWrite> make_view<int, Rank, Intent::ReadWrite>( const Array& );   \
    extern template ArrayView<long, Rank, Intent::ReadOnly> make_view<long, Rank, Intent::ReadOnly>( const Array& );   \
    extern template ArrayView<long, Rank, Intent::ReadWrite> make_view<long, Rank, Intent::ReadWrite>( const Array& ); \
    extern template ArrayView<float, Rank, Intent::ReadOnly> make_view<float, Rank, Intent::ReadOnly>( const Array& ); \
    extern template ArrayView<float, Rank, Intent::ReadWrite> make_view<float, Rank, Intent::ReadWrite>(               \
        const Array& );                                                                                                \
    extern template ArrayView<double, Rank, Intent::ReadOnly> make_view<double, Rank, Intent::ReadOnly>(               \
        const Array& );                                                                                                \
    extern template ArrayView<double, Rank, Intent::ReadWrite> make_view<double, Rank, Intent::ReadWrite>(             \
        const Array& );                                                                                                \
                                                                                                                       \
    extern template LocalView<int, Rank, Intent::ReadOnly> make_view<int, Rank, Intent::ReadOnly>(                     \
        const int data[], const ArrayShape& );                                                                         \
    extern template LocalView<int, Rank, Intent::ReadWrite> make_view<int, Rank, Intent::ReadWrite>(                   \
        const int data[], const ArrayShape& );                                                                         \
    extern template LocalView<long, Rank, Intent::ReadOnly> make_view<long, Rank, Intent::ReadOnly>(                   \
        const long data[], const ArrayShape& );                                                                        \
    extern template LocalView<long, Rank, Intent::ReadWrite> make_view<long, Rank, Intent::ReadWrite>(                 \
        const long data[], const ArrayShape& );                                                                        \
    extern template LocalView<float, Rank, Intent::ReadOnly> make_view<float, Rank, Intent::ReadOnly>(                 \
        const float data[], const ArrayShape& );                                                                       \
    extern template LocalView<float, Rank, Intent::ReadWrite> make_view<float, Rank, Intent::ReadWrite>(               \
        const float data[], const ArrayShape& );                                                                       \
    extern template LocalView<double, Rank, Intent::ReadOnly> make_view<double, Rank, Intent::ReadOnly>(               \
        const double data[], const ArrayShape& );                                                                      \
    extern template LocalView<double, Rank, Intent::ReadWrite> make_view<double, Rank, Intent::ReadWrite>(             \
        const double data[], const ArrayShape& );                                                                      \
                                                                                                                       \
    extern template LocalView<int, Rank, Intent::ReadOnly> make_view<int, Rank, Intent::ReadOnly>( const int data[],   \
                                                                                                   size_t );           \
    extern template LocalView<int, Rank, Intent::ReadWrite> make_view<int, Rank, Intent::ReadWrite>( const int data[], \
                                                                                                     size_t );         \
    extern template LocalView<long, Rank, Intent::ReadOnly> make_view<long, Rank, Intent::ReadOnly>(                   \
        const long data[], size_t );                                                                                   \
    extern template LocalView<long, Rank, Intent::ReadWrite> make_view<long, Rank, Intent::ReadWrite>(                 \
        const long data[], size_t );                                                                                   \
    extern template LocalView<float, Rank, Intent::ReadOnly> make_view<float, Rank, Intent::ReadOnly>(                 \
        const float data[], size_t );                                                                                  \
    extern template LocalView<float, Rank, Intent::ReadWrite> make_view<float, Rank, Intent::ReadWrite>(               \
        const float data[], size_t );                                                                                  \
    extern template LocalView<double, Rank, Intent::ReadOnly> make_view<double, Rank, Intent::ReadOnly>(               \
        const double data[], size_t );                                                                                 \
    extern template LocalView<double, Rank, Intent::ReadWrite> make_view<double, Rank, Intent::ReadWrite>(             \
        const double data[], size_t );

// For each NDims in [1..9]
EXPLICIT_TEMPLATE_INSTANTIATION( 1 )
EXPLICIT_TEMPLATE_INSTANTIATION( 2 )
EXPLICIT_TEMPLATE_INSTANTIATION( 3 )
EXPLICIT_TEMPLATE_INSTANTIATION( 4 )
EXPLICIT_TEMPLATE_INSTANTIATION( 5 )
EXPLICIT_TEMPLATE_INSTANTIATION( 6 )
EXPLICIT_TEMPLATE_INSTANTIATION( 7 )
EXPLICIT_TEMPLATE_INSTANTIATION( 8 )
EXPLICIT_TEMPLATE_INSTANTIATION( 9 )

#undef EXPLICIT_TEMPLATE_INSTANTIATION

}  // namespace array
}  // namespace atlas
