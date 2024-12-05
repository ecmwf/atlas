/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <sstream>

#include "atlas/array.h"
#include "atlas/array/ArrayView.h"
#include "atlas/array/IndexView.h"
#include "atlas/library/config.h"
#include "atlas/runtime/Exception.h"

namespace atlas {
namespace array {

namespace {
template <typename Value, int Rank>
inline static void check_metadata(const Array& array) {
    if (array.rank() != Rank) {
        std::stringstream err;
        err << "Number of dimensions do not match: template argument " << Rank << " expected to be " << array.rank();
        throw_Exception(err.str(), Here());
    }
    if (array.datatype() != array::DataType::create<Value>()) {
        std::stringstream err;
        err << "Data Type does not match: template argument expected to be " << array.datatype().str();
        throw_Exception(err.str(), Here());
    }
}
}  // namespace

//------------------------------------------------------------------------------

template <typename Value, int Rank>
ArrayView<Value, Rank> make_host_view(Array& array) {
    return ArrayView<Value, Rank>(array.host_data<Value>(), array.shape(), array.strides());
}

template <typename Value, int Rank>
ArrayView<const Value, Rank> make_host_view(const Array& array) {
    return ArrayView<const Value, Rank>(array.host_data<const Value>(), array.shape(), array.strides());
}

template <typename Value, int Rank>
ArrayView<Value, Rank> make_device_view(Array& array) {
#if ATLAS_HAVE_GPU
    ATLAS_ASSERT(array.deviceAllocated(),"make_device_view: Array not allocated on device");
    return ArrayView<Value, Rank>((array.device_data<Value>()), array.shape(), array.device_strides());
#else
    return make_host_view<Value, Rank>(array);
#endif
}

template <typename Value, int Rank>
ArrayView<const Value, Rank> make_device_view(const Array& array) {
#if ATLAS_HAVE_GPU
    ATLAS_ASSERT(array.deviceAllocated(),"make_device_view: Array not allocated on device");
    return ArrayView<const Value, Rank>(array.device_data<const Value>(), array.shape(), array.device_strides());
#else
    return make_host_view<const Value, Rank>(array);
#endif
}

template <typename Value, int Rank>
IndexView<Value, Rank> make_host_indexview(Array& array) {
    return IndexView<Value, Rank>((Value*)(array.storage()), array.shape().data());
}

template <typename Value, int Rank>
IndexView<const Value, Rank> make_host_indexview(const Array& array) {
    return IndexView<const Value, Rank>((const Value*)(array.storage()), array.shape().data());
}

template <typename Value, int Rank>
IndexView<Value, Rank> make_indexview(Array& array) {
    check_metadata<Value, Rank>(array);
    return make_host_indexview<Value, Rank>(array);
}

template <typename Value, int Rank>
IndexView<const Value, Rank> make_indexview(const Array& array) {
    check_metadata<Value, Rank>(array);
    return make_host_indexview<const Value, Rank>(array);
}

template <typename Value, int Rank>
ArrayView<Value, Rank> make_view(Array& array) {
    check_metadata<Value, Rank>(array);
    return make_host_view<Value, Rank>(array);
}

template <typename Value, int Rank>
ArrayView<const Value, Rank> make_view(const Array& array) {
    check_metadata<Value, Rank>(array);
    return make_host_view<const Value, Rank>(array);
}

// --------------------------------------------------------------------------------------------

}  // namespace array
}  // namespace atlas

//-----------------------------------------------------------------------
// Explicit instantiation
namespace atlas {
namespace array {

#define EXPLICIT_TEMPLATE_INSTANTIATION_TYPE_RANK(TYPE, RANK)                            \
    template ArrayView<TYPE, RANK> make_view<TYPE, RANK>(Array&);                        \
    template ArrayView<const TYPE, RANK> make_view<const TYPE, RANK>(Array&);            \
    template ArrayView<const TYPE, RANK> make_view<TYPE, RANK>(const Array&);            \
    template ArrayView<const TYPE, RANK> make_view<const TYPE, RANK>(const Array&);      \
                                                                                         \
    template ArrayView<TYPE, RANK> make_host_view<TYPE, RANK>(Array&);                   \
    template ArrayView<const TYPE, RANK> make_host_view<const TYPE, RANK>(Array&);       \
    template ArrayView<const TYPE, RANK> make_host_view<TYPE, RANK>(const Array&);       \
    template ArrayView<const TYPE, RANK> make_host_view<const TYPE, RANK>(const Array&); \
                                                                                         \
    template ArrayView<TYPE, RANK> make_device_view<TYPE, RANK>(Array&);                 \
    template ArrayView<const TYPE, RANK> make_device_view<const TYPE, RANK>(Array&);     \
    template ArrayView<const TYPE, RANK> make_device_view<TYPE, RANK>(const Array&);     \
    template ArrayView<const TYPE, RANK> make_device_view<const TYPE, RANK>(const Array&);

#define EXPLICIT_TEMPLATE_INSTATIATION(RANK)                \
    EXPLICIT_TEMPLATE_INSTANTIATION_TYPE_RANK(int, RANK)    \
    EXPLICIT_TEMPLATE_INSTANTIATION_TYPE_RANK(long, RANK)   \
    EXPLICIT_TEMPLATE_INSTANTIATION_TYPE_RANK(float, RANK)  \
    EXPLICIT_TEMPLATE_INSTANTIATION_TYPE_RANK(double, RANK) \
    EXPLICIT_TEMPLATE_INSTANTIATION_TYPE_RANK(unsigned int, RANK) \
    EXPLICIT_TEMPLATE_INSTANTIATION_TYPE_RANK(unsigned long, RANK)


EXPLICIT_TEMPLATE_INSTATIATION(1)
EXPLICIT_TEMPLATE_INSTATIATION(2)
EXPLICIT_TEMPLATE_INSTATIATION(3)
EXPLICIT_TEMPLATE_INSTATIATION(4)
EXPLICIT_TEMPLATE_INSTATIATION(5)
EXPLICIT_TEMPLATE_INSTATIATION(6)
EXPLICIT_TEMPLATE_INSTATIATION(7)
EXPLICIT_TEMPLATE_INSTATIATION(8)
EXPLICIT_TEMPLATE_INSTATIATION(9)

#undef EXPLICIT_TEMPLATE_INSTATIATION_TYPE_RANK
#undef EXPLICIT_TEMPLATE_INSTATIATION

#define EXPLICIT_TEMPLATE_INSTANTIATION_INDEXVIEW_TYPE_RANK(TYPE, RANK)            \
    template IndexView<TYPE, RANK> make_indexview<TYPE, RANK>(Array&);             \
    template IndexView<const TYPE, RANK> make_indexview<const TYPE, RANK>(Array&); \
    template IndexView<const TYPE, RANK> make_indexview<TYPE, RANK>(const Array&); \
    template IndexView<const TYPE, RANK> make_indexview<const TYPE, RANK>(const Array&);

#define EXPLICIT_TEMPLATE_INSTANTIATION_INDEXVIEW(RANK)            \
    EXPLICIT_TEMPLATE_INSTANTIATION_INDEXVIEW_TYPE_RANK(int, RANK) \
    EXPLICIT_TEMPLATE_INSTANTIATION_INDEXVIEW_TYPE_RANK(long, RANK)

EXPLICIT_TEMPLATE_INSTANTIATION_INDEXVIEW(1)
EXPLICIT_TEMPLATE_INSTANTIATION_INDEXVIEW(2)

#undef EXPLICIT_TEMPLATE_INSTANTIATION_INDEXVIEW
#undef EXPLICIT_TEMPLATE_INSTANTIATION_INDEXVIEW

}  // namespace array
}  // namespace atlas
