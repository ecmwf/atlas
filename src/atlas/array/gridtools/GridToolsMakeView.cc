/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "atlas/array/gridtools/GridToolsMakeView.h"

#include <sstream>
#include <vector>

#include "atlas/array/Array.h"
#include "atlas/array/ArrayView.h"
#include "atlas/array/ArrayViewDefs.h"
#include "atlas/array/IndexView.h"
#include "atlas/array/gridtools/GridToolsTraits.h"
#include "atlas/library/config.h"
#include "atlas/runtime/Exception.h"

//------------------------------------------------------------------------------

namespace atlas {
namespace array {

namespace {
template <typename Value, unsigned Rank>
static void check_metadata(const Array& array) {
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

namespace gridtools {


template <typename Value, unsigned int Rank>
typename gt_view<Value, Rank>::type make_gt_host_view(const Array& array) {
    using value_t         = typename std::remove_const<Value>::type;
    using storage_info_ty = storage_traits::storage_info_t<0, Rank>;
    using data_store_t    = storage_traits::data_store_t<value_t, storage_info_ty>;

    data_store_t* ds = reinterpret_cast<data_store_t*>(const_cast<void*>(array.storage()));
    return ::gridtools::make_host_view<get_access_mode<Value>()>(*ds);
}

template <typename Value, unsigned int Rank>
typename gt_view<Value, Rank>::type make_gt_device_view(const Array& array) {
    using value_t         = typename std::remove_const<Value>::type;
    using storage_info_ty = storage_traits::storage_info_t<0, Rank>;
    using data_store_t    = storage_traits::data_store_t<value_t, storage_info_ty>;

    data_store_t* ds = reinterpret_cast<data_store_t*>(const_cast<void*>(array.storage()));
#if ATLAS_GRIDTOOLS_STORAGE_BACKEND_CUDA
    return ::gridtools::make_device_view<get_access_mode<Value>()>(*ds);
#else
    return ::gridtools::make_host_view<get_access_mode<Value>()>(*ds);
#endif
}

}  // namespace gridtools


constexpr bool host_view   = false;
constexpr bool device_view = true;


template <typename Value, int Rank>
ArrayView<Value, Rank> make_host_view(Array& array) {
    check_metadata<Value, Rank>(array);
    return ArrayView<Value, Rank>(array, host_view);
}

template <typename Value, int Rank>
ArrayView<const Value, Rank> make_host_view(const Array& array) {
    check_metadata<Value, Rank>(array);
    return ArrayView<const Value, Rank>(array, host_view);
}

template <typename Value, int Rank>
ArrayView<Value, Rank> make_device_view(Array& array) {
    check_metadata<Value, Rank>(array);
    return ArrayView<Value, Rank>(array, device_view);
}

template <typename Value, int Rank>
ArrayView<const Value, Rank> make_device_view(const Array& array) {
    check_metadata<Value, Rank>(array);
    return ArrayView<const Value, Rank>(array, device_view);
}

template <typename Value, int Rank>
ArrayView<Value, Rank> make_view(Array& array) {
    return make_host_view<Value, Rank>(array);
}

template <typename Value, int Rank>
ArrayView<const Value, Rank> make_view(const Array& array) {
    return make_host_view<Value, Rank>(array);
}

//template <typename Value, int Rank>
//IndexView<Value, Rank> make_host_indexview( Array& array ) {
//    return IndexView<Value, Rank>( (Value*)( array.storage() ), array.shape().data() );
//}

//template <typename Value, int Rank>
//IndexView<const Value, Rank> make_host_indexview( const Array& array ) {
//    return IndexView<const Value, Rank>( (const Value*)( array.storage() ), array.shape().data() );
//}

//template <typename Value, int Rank>
//IndexView<Value, Rank> make_indexview( Array& array ) {
//    check_metadata<Value, Rank>( array );
//    return make_host_indexview<Value, Rank>( array );
//}

//template <typename Value, int Rank>
//IndexView<const Value, Rank> make_indexview( const Array& array ) {
//    check_metadata<Value, Rank>( array );
//    return make_host_indexview<Value, Rank>( array );
//}


template <typename Value, int Rank>
IndexView<Value, Rank> make_host_indexview(Array& array) {
    // using value_t         = typename std::remove_const<Value>::type;
    // using storage_info_ty = gridtools::storage_traits::storage_info_t<0, Rank>;
    // using data_store_t    = gridtools::storage_traits::data_store_t<value_t, storage_info_ty>;

    // data_store_t* ds = reinterpret_cast<data_store_t*>(const_cast<void*>(array.storage()));

    // return IndexView<Value, Rank>(::gridtools::make_host_view<gridtools::get_access_mode<Value>()>(*ds));
    check_metadata<Value, Rank>(array);
    constexpr bool device_view = false;
    return IndexView<Value, Rank>(array, device_view);
}

template <typename Value, int Rank>
IndexView<const Value, Rank> make_host_indexview(const Array& array) {
    // using value_t         = typename std::remove_const<Value>::type;
    // using storage_info_ty = gridtools::storage_traits::storage_info_t<0, Rank>;
    // using data_store_t    = gridtools::storage_traits::data_store_t<value_t, storage_info_ty>;

    // data_store_t* ds = reinterpret_cast<data_store_t*>(const_cast<void*>(array.storage()));

    // return IndexView<const Value, Rank>(::gridtools::make_host_view<gridtools::get_access_mode<const Value>()>(*ds));
    check_metadata<Value, Rank>(array);
    constexpr bool device_view = false;
    return IndexView<const Value, Rank>(array, device_view);
}


// --------------------------------------------------------------------------------------------

template <typename Value, int Rank>
IndexView<Value, Rank> make_indexview(Array& array) {
    // check_metadata<Value, Rank>(array);
    // constexpr bool device_view = false;
    // return IndexView<cValue, Rank>(array, device_view);
    return make_host_indexview<Value,Rank>(array);

}

template <typename Value, int Rank>
IndexView<const Value, Rank> make_indexview(const Array& array) {
    return make_host_indexview<const Value, Rank>(array);
}

}  // namespace array
}  // namespace atlas

//-----------------------------------------------------------------------
// Explicit instantiation
namespace atlas {
namespace array {

//#define EXPLICIT_TEMPLATE_INSTANTIATION( RANK )                                                                        \
//    namespace gridtools {                                                                                              \
//    template typename gt_view<long unsigned, RANK, Intent::ReadOnly>::type                                             \
//    make_gt_host_view<long unsigned, RANK, Intent::ReadOnly>( const Array& array );                                    \
//    template typename gt_view<long unsigned, RANK, Intent::ReadWrite>::type                                            \
//    make_gt_host_view<long unsigned, RANK, Intent::ReadWrite>( const Array& array );                                   \
//    template typename gt_view<float, RANK, Intent::ReadOnly>::type make_gt_host_view<float, RANK, Intent::ReadOnly>(   \
//        const Array& array );                                                                                          \
//    template typename gt_view<float, RANK, Intent::ReadWrite>::type make_gt_host_view<float, RANK, Intent::ReadWrite>( \
//        const Array& array );                                                                                          \
//    template typename gt_view<double, RANK, Intent::ReadOnly>::type make_gt_host_view<double, RANK, Intent::ReadOnly>( \
//        const Array& array );                                                                                          \
//    template typename gt_view<double, RANK, Intent::ReadWrite>::type                                                   \
//    make_gt_host_view<double, RANK, Intent::ReadWrite>( const Array& array );                                          \
//                                                                                                                       \
//    template typename gt_view<int, RANK, Intent::ReadOnly>::type make_gt_device_view<int, RANK, Intent::ReadOnly>(     \
//        const Array& array );                                                                                          \
//    template typename gt_view<int, RANK, Intent::ReadWrite>::type make_gt_device_view<int, RANK, Intent::ReadWrite>(   \
//        const Array& array );                                                                                          \
//    template typename gt_view<long, RANK, Intent::ReadOnly>::type make_gt_device_view<long, RANK, Intent::ReadOnly>(   \
//        const Array& array );                                                                                          \
//    template typename gt_view<long, RANK, Intent::ReadWrite>::type make_gt_device_view<long, RANK, Intent::ReadWrite>( \
//        const Array& array );                                                                                          \
//    template typename gt_view<long unsigned, RANK, Intent::ReadOnly>::type                                             \
//    make_gt_device_view<long unsigned, RANK, Intent::ReadOnly>( const Array& array );                                  \
//    template typename gt_view<long unsigned, RANK, Intent::ReadWrite>::type                                            \
//    make_gt_device_view<long unsigned, RANK, Intent::ReadWrite>( const Array& array );                                 \
//    template typename gt_view<float, RANK, Intent::ReadOnly>::type make_gt_device_view<float, RANK, Intent::ReadOnly>( \
//        const Array& array );                                                                                          \
//    template typename gt_view<float, RANK, Intent::ReadWrite>::type                                                    \
//    make_gt_device_view<float, RANK, Intent::ReadWrite>( const Array& array );                                         \
//    template typename gt_view<double, RANK, Intent::ReadOnly>::type                                                    \
//    make_gt_device_view<double, RANK, Intent::ReadOnly>( const Array& array );                                         \
//    template typename gt_view<double, RANK, Intent::ReadWrite>::type                                                   \
//    make_gt_device_view<double, RANK, Intent::ReadWrite>( const Array& array );                                        \
//    }

//// For each Rank in [1..9]
//EXPLICIT_TEMPLATE_INSTANTIATION( 1 )
//EXPLICIT_TEMPLATE_INSTANTIATION( 2 )
//EXPLICIT_TEMPLATE_INSTANTIATION( 3 )
//EXPLICIT_TEMPLATE_INSTANTIATION( 4 )
//EXPLICIT_TEMPLATE_INSTANTIATION( 5 )
//EXPLICIT_TEMPLATE_INSTANTIATION( 6 )
//EXPLICIT_TEMPLATE_INSTANTIATION( 7 )
//EXPLICIT_TEMPLATE_INSTANTIATION( 8 )
//EXPLICIT_TEMPLATE_INSTANTIATION( 9 )

//#undef EXPLICIT_TEMPLATE_INSTANTIATION


#define EXPLICIT_TEMPLATE_INSTANTIATION_TYPE_RANK(TYPE, RANK)                                                    \
    template ArrayView<TYPE, RANK> make_view<TYPE, RANK>(Array&);                                                \
    template ArrayView<const TYPE, RANK> make_view<const TYPE, RANK>(Array&);                                    \
    template ArrayView<const TYPE, RANK> make_view<TYPE, RANK>(const Array&);                                    \
    template ArrayView<const TYPE, RANK> make_view<const TYPE, RANK>(const Array&);                              \
                                                                                                                 \
    template ArrayView<TYPE, RANK> make_host_view<TYPE, RANK>(Array&);                                           \
    template ArrayView<const TYPE, RANK> make_host_view<const TYPE, RANK>(Array&);                               \
    template ArrayView<const TYPE, RANK> make_host_view<TYPE, RANK>(const Array&);                               \
    template ArrayView<const TYPE, RANK> make_host_view<const TYPE, RANK>(const Array&);                         \
                                                                                                                 \
    template ArrayView<TYPE, RANK> make_device_view<TYPE, RANK>(Array&);                                         \
    template ArrayView<const TYPE, RANK> make_device_view<const TYPE, RANK>(Array&);                             \
    template ArrayView<const TYPE, RANK> make_device_view<TYPE, RANK>(const Array&);                             \
    template ArrayView<const TYPE, RANK> make_device_view<const TYPE, RANK>(const Array&);                       \
                                                                                                                 \
    namespace gridtools {                                                                                        \
    template typename gt_view<TYPE, RANK>::type make_gt_host_view<TYPE, RANK>(const Array& array);               \
    template typename gt_view<const TYPE, RANK>::type make_gt_host_view<const TYPE, RANK>(const Array& array);   \
    template typename gt_view<TYPE, RANK>::type make_gt_device_view<TYPE, RANK>(const Array& array);             \
    template typename gt_view<const TYPE, RANK>::type make_gt_device_view<const TYPE, RANK>(const Array& array); \
    }


#define EXPLICIT_TEMPLATE_INSTATIATION(RANK)                \
    EXPLICIT_TEMPLATE_INSTANTIATION_TYPE_RANK(int, RANK)    \
    EXPLICIT_TEMPLATE_INSTANTIATION_TYPE_RANK(long, RANK)   \
    EXPLICIT_TEMPLATE_INSTANTIATION_TYPE_RANK(float, RANK)  \
    EXPLICIT_TEMPLATE_INSTANTIATION_TYPE_RANK(double, RANK) \
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

#define EXPLICIT_TEMPLATE_INSTANTIATION_INDEXVIEW_TYPE_RANK(TYPE, RANK)                       \
    template IndexView<TYPE, RANK> make_host_indexview<TYPE, RANK>(Array&);                   \
    template IndexView<const TYPE, RANK> make_host_indexview<const TYPE, RANK>(Array&);       \
    template IndexView<const TYPE, RANK> make_host_indexview<TYPE, RANK>(const Array&);       \
    template IndexView<const TYPE, RANK> make_host_indexview<const TYPE, RANK>(const Array&); \
                                                                                              \
    template IndexView<TYPE, RANK> make_indexview<TYPE, RANK>(Array&);                        \
    template IndexView<const TYPE, RANK> make_indexview<const TYPE, RANK>(Array&);            \
    template IndexView<const TYPE, RANK> make_indexview<TYPE, RANK>(const Array&);            \
    template IndexView<const TYPE, RANK> make_indexview<const TYPE, RANK>(const Array&);

#define EXPLICIT_TEMPLATE_INSTANTIATION_INDEXVIEW(RANK)            \
    EXPLICIT_TEMPLATE_INSTANTIATION_INDEXVIEW_TYPE_RANK(int, RANK) \
    EXPLICIT_TEMPLATE_INSTANTIATION_INDEXVIEW_TYPE_RANK(long, RANK)

EXPLICIT_TEMPLATE_INSTANTIATION_INDEXVIEW(1)
EXPLICIT_TEMPLATE_INSTANTIATION_INDEXVIEW(2)

#undef EXPLICIT_TEMPLATE_INSTANTIATION_INDEXVIEW
#undef EXPLICIT_TEMPLATE_INSTANTIATION_INDEXVIEW_TYPE_RANK


}  // namespace array
}  // namespace atlas
