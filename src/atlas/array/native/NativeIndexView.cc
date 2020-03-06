/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <iostream>

#include "atlas/array/native/NativeIndexView.h"

//------------------------------------------------------------------------------------------------------

namespace atlas {
namespace array {

//------------------------------------------------------------------------------------------------------

template <typename Value, int Rank>
IndexView<Value, Rank>::IndexView( Value* data, const idx_t shape[1] ) : data_( data ) {
    strides_[0] = 1;
    shape_[0]   = shape[0];
}

template <typename Value, int Rank>
IndexView<Value, Rank>::IndexView( Value* data, const idx_t shape[1], const idx_t strides[1] ) :
    data_( const_cast<Value*>( data ) ) {
    strides_[0] = strides[0];
    shape_[0]   = shape[0];
}

template <typename Value, int Rank>
void IndexView<Value, Rank>::dump( std::ostream& os ) const {
    os << "size: " << size() << " , values: ";
    os << "[ ";
    for ( idx_t j = 0; j < size(); ++j ) {
        os << ( *this )( j ) << " ";
    }
    os << "]" << std::endl;
}

//------------------------------------------------------------------------------------------------------
// Explicit template instatiation

#define EXPLICIT_TEMPLATE_INSTANTIATION_TYPE_RANK( TYPE, RANK ) \
    template class IndexView<TYPE, RANK>;                       \
    template class IndexView<const TYPE, RANK>;

#define EXPLICIT_TEMPLATE_INSTANTIATION( RANK )            \
    EXPLICIT_TEMPLATE_INSTANTIATION_TYPE_RANK( int, RANK ) \
    EXPLICIT_TEMPLATE_INSTANTIATION_TYPE_RANK( long, RANK )

EXPLICIT_TEMPLATE_INSTANTIATION( 1 )
EXPLICIT_TEMPLATE_INSTANTIATION( 2 )

#undef EXPLICIT_TEMPLATE_INSTANTIATION_TYPE_RANK
#undef EXPLICIT_TEMPLATE_INSTANTIATION

//------------------------------------------------------------------------------------------------------

}  // namespace array
}  // namespace atlas
