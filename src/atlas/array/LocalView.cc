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

#define ENABLE_IF_NON_CONST \
    template <bool EnableBool, typename std::enable_if<( !std::is_const<Value>::value && EnableBool ), int>::type*>


template <typename Value, int Rank>
ENABLE_IF_NON_CONST void LocalView<Value, Rank>::assign( const Value& value ) {
    helpers::array_assigner<Value, Rank>::apply( *this, value );
}

//------------------------------------------------------------------------------------------------------

template <typename Value, int Rank>
void LocalView<Value, Rank>::dump( std::ostream& os ) const {
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


#define ENABLE_IF_NOT_CONST typename std::enable_if<!std::is_const<Value>::value, Value>::type*
#define ENABLE_IF_CONST typename std::enable_if<std::is_const<Value>::value, Value>::type*

template <class Value, int Rank, ENABLE_IF_NOT_CONST>
LocalView<Value, Rank> make_view( Value data[], const ArrayShape& shape ) {
    return LocalView<Value, Rank>( data, shape );
}

template <class Value, int Rank, ENABLE_IF_NOT_CONST>
LocalView<const Value, Rank> make_view( const Value data[], const ArrayShape& shape ) {
    return LocalView<const Value, Rank>( data, shape );
}

template <class Value, int Rank, ENABLE_IF_CONST>
LocalView<Value, Rank> make_view( Value data[], const ArrayShape& shape ) {
    return LocalView<Value, Rank>( data, shape );
}

template <class Value, int Rank, ENABLE_IF_CONST>
LocalView<Value, Rank> make_view( typename std::remove_const<Value>::type data[], const ArrayShape& shape ) {
    return LocalView<Value, Rank>( data, shape );
}

//------------------------------------------------------------------------------------------------------

template <typename Value, unsigned int Rank, ENABLE_IF_NOT_CONST>
LocalView<Value, Rank> make_view( Value data[], size_t size ) {
    return LocalView<Value, Rank>( data, ArrayShape{idx_t( size )} );
}

template <typename Value, unsigned int Rank, ENABLE_IF_NOT_CONST>
LocalView<const Value, Rank> make_view( const Value data[], size_t size ) {
    return LocalView<const Value, Rank>( data, ArrayShape{idx_t( size )} );
}

template <typename Value, unsigned int Rank, ENABLE_IF_CONST>
LocalView<Value, Rank> make_view( Value data[], size_t size ) {
    return LocalView<Value, Rank>( data, ArrayShape{idx_t( size )} );
}

template <typename Value, unsigned int Rank, ENABLE_IF_CONST>
LocalView<Value, Rank> make_view( typename std::remove_const<Value>::type data[], size_t size ) {
    return LocalView<Value, Rank>( data, ArrayShape{idx_t( size )} );
}

//------------------------------------------------------------------------------------------------------

}  // namespace array
}  // namespace atlas

//-----------------------------------------------------------------------
// Explicit instantiation
namespace atlas {
namespace array {
#undef EXPLICIT_TEMPLATE_INSTANTIATION

#define EXPLICIT_TEMPLATE_INSTANTIATION_TYPE_RANK( TYPE, RANK )                                                        \
    template LocalView<TYPE, RANK> make_view<TYPE, RANK, nullptr>( TYPE data[], const ArrayShape& );                   \
    template LocalView<const TYPE, RANK> make_view<const TYPE, RANK, nullptr>( TYPE data[], const ArrayShape& );       \
    template LocalView<const TYPE, RANK> make_view<TYPE, RANK, nullptr>( const TYPE data[], const ArrayShape& );       \
    template LocalView<const TYPE, RANK> make_view<const TYPE, RANK, nullptr>( const TYPE data[], const ArrayShape& ); \
                                                                                                                       \
    template class LocalView<TYPE, RANK>;                                                                              \
    template class LocalView<const TYPE, RANK>;                                                                        \
    template void LocalView<TYPE, RANK>::assign<true, nullptr>( const TYPE& );

#define EXPLICIT_TEMPLATE_INSTANTIATION_TYPE_RANK1( TYPE )                                      \
    template LocalView<TYPE, 1> make_view<TYPE, 1, nullptr>( TYPE data[], size_t );             \
    template LocalView<const TYPE, 1> make_view<const TYPE, 1, nullptr>( TYPE data[], size_t ); \
    template LocalView<const TYPE, 1> make_view<TYPE, 1, nullptr>( const TYPE data[], size_t ); \
    template LocalView<const TYPE, 1> make_view<const TYPE, 1, nullptr>( const TYPE data[], size_t );

#define EXPLICIT_TEMPLATE_INSTANTIATION( RANK )              \
    EXPLICIT_TEMPLATE_INSTANTIATION_TYPE_RANK( int, RANK )   \
    EXPLICIT_TEMPLATE_INSTANTIATION_TYPE_RANK( long, RANK )  \
    EXPLICIT_TEMPLATE_INSTANTIATION_TYPE_RANK( float, RANK ) \
    EXPLICIT_TEMPLATE_INSTANTIATION_TYPE_RANK( double, RANK )


// For each RANK in [1..9]
EXPLICIT_TEMPLATE_INSTANTIATION( 1 )
EXPLICIT_TEMPLATE_INSTANTIATION( 2 )
EXPLICIT_TEMPLATE_INSTANTIATION( 3 )
EXPLICIT_TEMPLATE_INSTANTIATION( 4 )
EXPLICIT_TEMPLATE_INSTANTIATION( 5 )
EXPLICIT_TEMPLATE_INSTANTIATION( 6 )
EXPLICIT_TEMPLATE_INSTANTIATION( 7 )
EXPLICIT_TEMPLATE_INSTANTIATION( 8 )
EXPLICIT_TEMPLATE_INSTANTIATION( 9 )
EXPLICIT_TEMPLATE_INSTANTIATION_TYPE_RANK1( int )
EXPLICIT_TEMPLATE_INSTANTIATION_TYPE_RANK1( long )
EXPLICIT_TEMPLATE_INSTANTIATION_TYPE_RANK1( float )
EXPLICIT_TEMPLATE_INSTANTIATION_TYPE_RANK1( double )


#undef EXPLICIT_TEMPLATE_INSTANTIATION

}  // namespace array
}  // namespace atlas
