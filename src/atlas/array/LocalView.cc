/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "LocalView.h"

#include <iostream>

#include "atlas/array/helpers/ArrayAssigner.h"
#include "atlas/runtime/Exception.h"

//------------------------------------------------------------------------------------------------------

namespace atlas {
namespace array {

//------------------------------------------------------------------------------------------------------

template <typename Value, int Rank, Intent AccessMode>
void LocalView<Value, Rank, AccessMode>::assign( const value_type& value ) {
    helpers::array_assigner<Value, Rank>::apply( *this, value );
}

//------------------------------------------------------------------------------------------------------

template <typename Value, int Rank, Intent AccessMode>
void LocalView<Value, Rank, AccessMode>::dump( std::ostream& os ) const {
    ATLAS_ASSERT( contiguous(), "Cannot dump non-contiguous view" );
    const value_type* data_ = data();
    os << "size: " << size() << " , values: ";
    os << "[ ";
    for ( idx_t j = 0; j < size(); ++j ) {
        os << data_[j] << " ";
    }
    os << "]";
}

//------------------------------------------------------------------------------------------------------

template <typename Value, unsigned int Rank, Intent AccessMode>
LocalView<Value, Rank, AccessMode> make_view( const Value data[], const ArrayShape& shape ) {
    return LocalView<Value, Rank, AccessMode>( data, shape );
}

//------------------------------------------------------------------------------------------------------

template <typename Value, unsigned int Rank, Intent AccessMode>
LocalView<Value, Rank, AccessMode> make_view( const Value data[], size_t size ) {
    return LocalView<Value, Rank, AccessMode>( data, ArrayShape{idx_t( size )} );
}

//------------------------------------------------------------------------------------------------------

}  // namespace array
}  // namespace atlas

//-----------------------------------------------------------------------
// Explicit instantiation
namespace atlas {
namespace array {
#define EXPLICIT_TEMPLATE_INSTANTIATION( Rank )                                                                        \
    template class LocalView<int, Rank, Intent::ReadOnly>;                                                             \
    template class LocalView<int, Rank, Intent::ReadWrite>;                                                            \
    template class LocalView<long, Rank, Intent::ReadOnly>;                                                            \
    template class LocalView<long, Rank, Intent::ReadWrite>;                                                           \
    template class LocalView<long unsigned, Rank, Intent::ReadOnly>;                                                   \
    template class LocalView<long unsigned, Rank, Intent::ReadWrite>;                                                  \
    template class LocalView<float, Rank, Intent::ReadOnly>;                                                           \
    template class LocalView<float, Rank, Intent::ReadWrite>;                                                          \
    template class LocalView<double, Rank, Intent::ReadOnly>;                                                          \
    template class LocalView<double, Rank, Intent::ReadWrite>;                                                         \
    template LocalView<int, Rank, Intent::ReadOnly> make_view<int, Rank, Intent::ReadOnly>( const int data[],          \
                                                                                            const ArrayShape& );       \
    template LocalView<int, Rank, Intent::ReadWrite> make_view<int, Rank, Intent::ReadWrite>( const int data[],        \
                                                                                              const ArrayShape& );     \
    template LocalView<long, Rank, Intent::ReadOnly> make_view<long, Rank, Intent::ReadOnly>( const long data[],       \
                                                                                              const ArrayShape& );     \
    template LocalView<long, Rank, Intent::ReadWrite> make_view<long, Rank, Intent::ReadWrite>( const long data[],     \
                                                                                                const ArrayShape& );   \
    template LocalView<float, Rank, Intent::ReadOnly> make_view<float, Rank, Intent::ReadOnly>( const float data[],    \
                                                                                                const ArrayShape& );   \
    template LocalView<float, Rank, Intent::ReadWrite> make_view<float, Rank, Intent::ReadWrite>( const float data[],  \
                                                                                                  const ArrayShape& ); \
    template LocalView<double, Rank, Intent::ReadOnly> make_view<double, Rank, Intent::ReadOnly>( const double data[], \
                                                                                                  const ArrayShape& ); \
    template LocalView<double, Rank, Intent::ReadWrite> make_view<double, Rank, Intent::ReadWrite>(                    \
        const double data[], const ArrayShape& );                                                                      \
                                                                                                                       \
    template LocalView<int, Rank, Intent::ReadOnly> make_view<int, Rank, Intent::ReadOnly>( const int data[],          \
                                                                                            size_t );                  \
    template LocalView<int, Rank, Intent::ReadWrite> make_view<int, Rank, Intent::ReadWrite>( const int data[],        \
                                                                                              size_t );                \
    template LocalView<long, Rank, Intent::ReadOnly> make_view<long, Rank, Intent::ReadOnly>( const long data[],       \
                                                                                              size_t );                \
    template LocalView<long, Rank, Intent::ReadWrite> make_view<long, Rank, Intent::ReadWrite>( const long data[],     \
                                                                                                size_t );              \
    template LocalView<float, Rank, Intent::ReadOnly> make_view<float, Rank, Intent::ReadOnly>( const float data[],    \
                                                                                                size_t );              \
    template LocalView<float, Rank, Intent::ReadWrite> make_view<float, Rank, Intent::ReadWrite>( const float data[],  \
                                                                                                  size_t );            \
    template LocalView<double, Rank, Intent::ReadOnly> make_view<double, Rank, Intent::ReadOnly>( const double data[], \
                                                                                                  size_t );            \
    template LocalView<double, Rank, Intent::ReadWrite> make_view<double, Rank, Intent::ReadWrite>(                    \
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
