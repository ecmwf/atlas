/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <initializer_list>
#include <ostream>

#include "atlas/array/ArrayView.h"
#include "atlas/array/helpers/ArrayAssigner.h"
#include "atlas/array/helpers/ArrayWriter.h"

//------------------------------------------------------------------------------------------------------

namespace atlas {
namespace array {

//------------------------------------------------------------------------------------------------------

template <typename Value, int Rank, Intent AccessMode>
void ArrayView<Value, Rank, AccessMode>::assign( const value_type& value ) {
    helpers::array_assigner<Value, Rank>::apply( *this, value );
}

//------------------------------------------------------------------------------------------------------

template <typename Value, int Rank, Intent AccessMode>
void ArrayView<Value, Rank, AccessMode>::assign( const std::initializer_list<value_type>& list ) {
    helpers::array_assigner<Value, Rank>::apply( *this, list );
}

//------------------------------------------------------------------------------------------------------

template <typename Value, int Rank, Intent AccessMode>
void ArrayView<Value, Rank, AccessMode>::dump( std::ostream& os ) const {
    os << "size: " << size() << " , values: ";
    os << "[ ";
    helpers::array_writer::apply( *this, os );
    os << " ]";
}

//------------------------------------------------------------------------------------------------------

}  // namespace array
}  // namespace atlas

//-----------------------------------------------------------------------
// Explicit instantiation
namespace atlas {
namespace array {
#define EXPLICIT_TEMPLATE_INSTANTIATION( Rank )                       \
    template class ArrayView<int, Rank, Intent::ReadOnly>;            \
    template class ArrayView<int, Rank, Intent::ReadWrite>;           \
    template class ArrayView<long, Rank, Intent::ReadOnly>;           \
    template class ArrayView<long, Rank, Intent::ReadWrite>;          \
    template class ArrayView<long unsigned, Rank, Intent::ReadOnly>;  \
    template class ArrayView<long unsigned, Rank, Intent::ReadWrite>; \
    template class ArrayView<float, Rank, Intent::ReadOnly>;          \
    template class ArrayView<float, Rank, Intent::ReadWrite>;         \
    template class ArrayView<double, Rank, Intent::ReadOnly>;         \
    template class ArrayView<double, Rank, Intent::ReadWrite>;

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
