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
#include "atlas/array/helpers/ArrayCopier.h"
#include "atlas/array/helpers/ArrayWriter.h"

//------------------------------------------------------------------------------------------------------

namespace atlas {
namespace array {

//------------------------------------------------------------------------------------------------------

#undef ENABLE_IF_NON_CONST
#define ENABLE_IF_NON_CONST \
    template <bool EnableBool, typename std::enable_if<(!std::is_const<Value>::value && EnableBool), int>::type*>


template <typename Value, int Rank>
ENABLE_IF_NON_CONST void ArrayView<Value, Rank>::assign(const value_type& value) {
    helpers::array_assigner<Value, Rank>::apply(*this, value);
}

//------------------------------------------------------------------------------------------------------

template <typename Value, int Rank>
ENABLE_IF_NON_CONST void ArrayView<Value, Rank>::assign(const std::initializer_list<value_type>& list) {
    helpers::array_assigner<Value, Rank>::apply(*this, list);
}

//------------------------------------------------------------------------------------------------------


template <typename Value, int Rank>
ENABLE_IF_NON_CONST void ArrayView<Value, Rank>::assign(const ArrayView& other) {
    helpers::array_copier<Value, Rank>::apply(other, *this);
}

//------------------------------------------------------------------------------------------------------

template <typename Value, int Rank>
ENABLE_IF_NON_CONST void ArrayView<Value, Rank>::assign(const LocalView<Value, Rank>& other) {
    helpers::array_copier<Value, Rank>::apply(other, *this);
}

//------------------------------------------------------------------------------------------------------

template <typename Value, int Rank>
void ArrayView<Value, Rank>::dump(std::ostream& os) const {
    os << "size: " << size() << " , values: ";
    os << "[ ";
    helpers::array_writer::apply(*this, os);
    os << " ]";
}

//------------------------------------------------------------------------------------------------------

}  // namespace array
}  // namespace atlas

//-----------------------------------------------------------------------
// Explicit instantiation
namespace atlas {
namespace array {
#define EXPLICIT_TEMPLATE_INSTANTIATION(Rank)                                                           \
    template class ArrayView<int, Rank>;                                                                \
    template class ArrayView<const int, Rank>;                                                          \
    template class ArrayView<long, Rank>;                                                               \
    template class ArrayView<const long, Rank>;                                                         \
    template class ArrayView<unsigned int, Rank>;                                                       \
    template class ArrayView<const unsigned int, Rank>;                                                 \
    template class ArrayView<long unsigned, Rank>;                                                      \
    template class ArrayView<const long unsigned, Rank>;                                                \
    template class ArrayView<float, Rank>;                                                              \
    template class ArrayView<const float, Rank>;                                                        \
    template class ArrayView<double, Rank>;                                                             \
    template class ArrayView<const double, Rank>;                                                       \
    template void ArrayView<int, Rank>::assign<true, nullptr>(int const&);                              \
    template void ArrayView<long, Rank>::assign<true, nullptr>(long const&);                            \
    template void ArrayView<float, Rank>::assign<true, nullptr>(float const&);                          \
    template void ArrayView<double, Rank>::assign<true, nullptr>(double const&);                        \
    template void ArrayView<int, Rank>::assign<true, nullptr>(std::initializer_list<int> const&);       \
    template void ArrayView<long, Rank>::assign<true, nullptr>(std::initializer_list<long> const&);     \
    template void ArrayView<float, Rank>::assign<true, nullptr>(std::initializer_list<float> const&);   \
    template void ArrayView<double, Rank>::assign<true, nullptr>(std::initializer_list<double> const&); \
    template void ArrayView<int, Rank>::assign<true, nullptr>(ArrayView<int, Rank> const&);             \
    template void ArrayView<long, Rank>::assign<true, nullptr>(ArrayView<long, Rank> const&);           \
    template void ArrayView<float, Rank>::assign<true, nullptr>(ArrayView<float, Rank> const&);         \
    template void ArrayView<double, Rank>::assign<true, nullptr>(ArrayView<double, Rank> const&);       \
    template void ArrayView<int, Rank>::assign<true, nullptr>(LocalView<int, Rank> const&);             \
    template void ArrayView<long, Rank>::assign<true, nullptr>(LocalView<long, Rank> const&);           \
    template void ArrayView<float, Rank>::assign<true, nullptr>(LocalView<float, Rank> const&);         \
    template void ArrayView<double, Rank>::assign<true, nullptr>(LocalView<double, Rank> const&);

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

}  // namespace array
}  // namespace atlas
