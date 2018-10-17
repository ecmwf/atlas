#pragma once

#include "atlas/array/ArrayView.h"
#include "atlas/array/IndexView.h"
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
        const Array& );

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
