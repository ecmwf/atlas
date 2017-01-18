/*
 * (C) Copyright 1996-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "atlas/array/LocalView.h"
#include "eckit/exception/Exceptions.h"

//------------------------------------------------------------------------------------------------------

namespace atlas {
namespace array {

//------------------------------------------------------------------------------------------------------

template <typename Value, int Rank>
void LocalView<Value,Rank>::assign(const value_type& value) {
    ASSERT( contiguous() );
    value_type* raw_data = data();
    for( size_t j=0; j<size_; ++j ) {
        raw_data[j] = value;
    }
}

//------------------------------------------------------------------------------------------------------

template <typename Value, int Rank>
void LocalView<Value,Rank>::dump(std::ostream& os) const {
ASSERT( contiguous() );
const value_type* data_ = data();
os << "size: " << size() << " , values: ";
os << "[ ";
for( size_t j=0; j<size(); ++ j )
  os << data_[j] << " ";
os << "]";
}


template <int Dim>
static char array_dim();

//------------------------------------------------------------------------------------------------------

} // namespace array
} // namespace atlas


//-----------------------------------------------------------------------
// Explicit instantiation
namespace atlas {
namespace array {
#define EXPLICIT_TEMPLATE_INSTANTIATION(Rank) \
template class LocalView<int,Rank>;\
template class LocalView<long,Rank>;\
template class LocalView<long unsigned,Rank>;\
template class LocalView<float,Rank>;\
template class LocalView<double,Rank>;\

// For each NDims in [1..9]
EXPLICIT_TEMPLATE_INSTANTIATION(1)
EXPLICIT_TEMPLATE_INSTANTIATION(2)
EXPLICIT_TEMPLATE_INSTANTIATION(3)
EXPLICIT_TEMPLATE_INSTANTIATION(4)
EXPLICIT_TEMPLATE_INSTANTIATION(5)
EXPLICIT_TEMPLATE_INSTANTIATION(6)
EXPLICIT_TEMPLATE_INSTANTIATION(7)
EXPLICIT_TEMPLATE_INSTANTIATION(8)
EXPLICIT_TEMPLATE_INSTANTIATION(9)

#undef EXPLICIT_TEMPLATE_INSTANTIATION
}
}
