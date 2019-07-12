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
static void check_metadata( const Array& array ) {
    if ( array.rank() != Rank ) {
        std::stringstream err;
        err << "Number of dimensions do not match: template argument " << Rank << " expected to be " << array.rank();
        throw_Exception( err.str(), Here() );
    }
    if ( array.datatype() != array::DataType::create<Value>() ) {
        std::stringstream err;
        err << "Data Type does not match: template argument expected to be " << array.datatype().str();
        throw_Exception( err.str(), Here() );
    }
}
}  // namespace

namespace gridtools {


template <typename Value, unsigned int Rank, Intent AccessMode>
typename gt_view<Value, Rank, AccessMode>::type make_gt_host_view( const Array& array ) {
    using storage_info_ty = storage_traits::storage_info_t<0, Rank>;
    using data_store_t    = storage_traits::data_store_t<Value, storage_info_ty>;

    data_store_t* ds = reinterpret_cast<data_store_t*>( const_cast<void*>( array.storage() ) );
    return ::gridtools::make_host_view<get_access_mode( AccessMode )>( *ds );
}

template <typename Value, unsigned int Rank, Intent AccessMode>
typename gt_view<Value, Rank, AccessMode>::type make_gt_device_view( const Array& array ) {
    using storage_info_ty = storage_traits::storage_info_t<0, Rank>;
    using data_store_t    = storage_traits::data_store_t<Value, storage_info_ty>;

    data_store_t* ds = reinterpret_cast<data_store_t*>( const_cast<void*>( array.storage() ) );
#if ATLAS_GRIDTOOLS_STORAGE_BACKEND_CUDA
    return ::gridtools::make_device_view<get_access_mode( AccessMode )>( *ds );
#else
    return ::gridtools::make_host_view<get_access_mode( AccessMode )>( *ds );
#endif
}

}  // namespace gridtools

template <typename Value, unsigned int Rank, Intent AccessMode>
ArrayView<Value, Rank, AccessMode> make_host_view( const Array& array ) {
    check_metadata<Value, Rank>( array );
    return ArrayView<Value, Rank, AccessMode>( gridtools::make_gt_host_view<Value, Rank, AccessMode>( array ), array );
}

template <typename Value, unsigned int Rank, Intent AccessMode>
ArrayView<Value, Rank, AccessMode> make_device_view( const Array& array ) {
    check_metadata<Value, Rank>( array );
    return ArrayView<Value, Rank, AccessMode>( gridtools::make_gt_device_view<Value, Rank, AccessMode>( array ),
                                               array );
}

template <typename Value, unsigned int Rank, Intent AccessMode>
IndexView<Value, Rank> make_host_indexview( const Array& array ) {
    typedef gridtools::storage_traits::storage_info_t<0, Rank> storage_info_ty;
    typedef gridtools::storage_traits::data_store_t<Value, storage_info_ty> data_store_t;

    data_store_t* ds = reinterpret_cast<data_store_t*>( const_cast<void*>( array.storage() ) );

    return IndexView<Value, Rank>( ::gridtools::make_host_view<::gridtools::access_mode::ReadWrite>( *ds ) );
}

// --------------------------------------------------------------------------------------------

template <typename Value, unsigned int Rank, Intent AccessMode>
IndexView<Value, Rank> make_indexview( const Array& array ) {
    check_metadata<Value, Rank>( array );
    return make_host_indexview<Value, Rank>( array );
}

template <typename Value, unsigned int Rank, Intent AccessMode>
ArrayView<Value, Rank, AccessMode> make_view( const Array& array ) {
    check_metadata<Value, Rank>( array );

    return make_host_view<Value, Rank, AccessMode>( array );
}

}  // namespace array
}  // namespace atlas

//-----------------------------------------------------------------------
// Explicit instantiation
namespace atlas {
namespace array {

#define EXPLICIT_TEMPLATE_INSTANTIATION( RANK )                                                                        \
    template ArrayView<int, RANK, Intent::ReadOnly> make_view<int, RANK, Intent::ReadOnly>( const Array& );            \
    template ArrayView<int, RANK, Intent::ReadWrite> make_view<int, RANK, Intent::ReadWrite>( const Array& );          \
    template ArrayView<long, RANK, Intent::ReadOnly> make_view<long, RANK, Intent::ReadOnly>( const Array& );          \
    template ArrayView<long, RANK, Intent::ReadWrite> make_view<long, RANK, Intent::ReadWrite>( const Array& );        \
    template ArrayView<long unsigned, RANK, Intent::ReadOnly> make_view<long unsigned, RANK, Intent::ReadOnly>(        \
        const Array& );                                                                                                \
    template ArrayView<long unsigned, RANK, Intent::ReadWrite> make_view<long unsigned, RANK, Intent::ReadWrite>(      \
        const Array& );                                                                                                \
    template ArrayView<float, RANK, Intent::ReadOnly> make_view<float, RANK, Intent::ReadOnly>( const Array& );        \
    template ArrayView<float, RANK, Intent::ReadWrite> make_view<float, RANK, Intent::ReadWrite>( const Array& );      \
    template ArrayView<double, RANK, Intent::ReadOnly> make_view<double, RANK, Intent::ReadOnly>( const Array& );      \
    template ArrayView<double, RANK, Intent::ReadWrite> make_view<double, RANK, Intent::ReadWrite>( const Array& );    \
                                                                                                                       \
    template ArrayView<int, RANK, Intent::ReadOnly> make_host_view<int, RANK, Intent::ReadOnly>( const Array& );       \
    template ArrayView<int, RANK, Intent::ReadWrite> make_host_view<int, RANK, Intent::ReadWrite>( const Array& );     \
    template ArrayView<long, RANK, Intent::ReadOnly> make_host_view<long, RANK, Intent::ReadOnly>( const Array& );     \
    template ArrayView<long, RANK, Intent::ReadWrite> make_host_view<long, RANK, Intent::ReadWrite>( const Array& );   \
    template ArrayView<long unsigned, RANK, Intent::ReadOnly> make_host_view<long unsigned, RANK, Intent::ReadOnly>(   \
        const Array& );                                                                                                \
    template ArrayView<long unsigned, RANK, Intent::ReadWrite> make_host_view<long unsigned, RANK, Intent::ReadWrite>( \
        const Array& );                                                                                                \
    template ArrayView<float, RANK, Intent::ReadOnly> make_host_view<float, RANK, Intent::ReadOnly>( const Array& );   \
    template ArrayView<float, RANK, Intent::ReadWrite> make_host_view<float, RANK, Intent::ReadWrite>( const Array& ); \
    template ArrayView<double, RANK, Intent::ReadOnly> make_host_view<double, RANK, Intent::ReadOnly>( const Array& ); \
    template ArrayView<double, RANK, Intent::ReadWrite> make_host_view<double, RANK, Intent::ReadWrite>(               \
        const Array& );                                                                                                \
                                                                                                                       \
    template ArrayView<int, RANK, Intent::ReadOnly> make_device_view<int, RANK, Intent::ReadOnly>( const Array& );     \
    template ArrayView<int, RANK, Intent::ReadWrite> make_device_view<int, RANK, Intent::ReadWrite>( const Array& );   \
    template ArrayView<long, RANK, Intent::ReadOnly> make_device_view<long, RANK, Intent::ReadOnly>( const Array& );   \
    template ArrayView<long, RANK, Intent::ReadWrite> make_device_view<long, RANK, Intent::ReadWrite>( const Array& ); \
    template ArrayView<long unsigned, RANK, Intent::ReadOnly> make_device_view<long unsigned, RANK, Intent::ReadOnly>( \
        const Array& );                                                                                                \
    template ArrayView<long unsigned, RANK, Intent::ReadWrite>                                                         \
    make_device_view<long unsigned, RANK, Intent::ReadWrite>( const Array& );                                          \
    template ArrayView<float, RANK, Intent::ReadOnly> make_device_view<float, RANK, Intent::ReadOnly>( const Array& ); \
    template ArrayView<float, RANK, Intent::ReadWrite> make_device_view<float, RANK, Intent::ReadWrite>(               \
        const Array& );                                                                                                \
    template ArrayView<double, RANK, Intent::ReadOnly> make_device_view<double, RANK, Intent::ReadOnly>(               \
        const Array& );                                                                                                \
    template ArrayView<double, RANK, Intent::ReadWrite> make_device_view<double, RANK, Intent::ReadWrite>(             \
        const Array& );                                                                                                \
                                                                                                                       \
    template IndexView<int, RANK> make_indexview<int, RANK, Intent::ReadOnly>( const Array& );                         \
    template IndexView<int, RANK> make_indexview<int, RANK, Intent::ReadWrite>( const Array& );                        \
                                                                                                                       \
    template IndexView<int, RANK> make_host_indexview<int, RANK, Intent::ReadOnly>( const Array& );                    \
    template IndexView<int, RANK> make_host_indexview<int, RANK, Intent::ReadWrite>( const Array& );                   \
                                                                                                                       \
    template IndexView<long, RANK> make_indexview<long, RANK, Intent::ReadOnly>( const Array& );                       \
    template IndexView<long, RANK> make_indexview<long, RANK, Intent::ReadWrite>( const Array& );                      \
                                                                                                                       \
    template IndexView<long, RANK> make_host_indexview<long, RANK, Intent::ReadOnly>( const Array& );                  \
    template IndexView<long, RANK> make_host_indexview<long, RANK, Intent::ReadWrite>( const Array& );                 \
    namespace gridtools {                                                                                              \
    template typename gt_view<int, RANK, Intent::ReadOnly>::type make_gt_host_view<int, RANK, Intent::ReadOnly>(       \
        const Array& array );                                                                                          \
    template typename gt_view<int, RANK, Intent::ReadWrite>::type make_gt_host_view<int, RANK, Intent::ReadWrite>(     \
        const Array& array );                                                                                          \
    template typename gt_view<long, RANK, Intent::ReadOnly>::type make_gt_host_view<long, RANK, Intent::ReadOnly>(     \
        const Array& array );                                                                                          \
    template typename gt_view<long, RANK, Intent::ReadWrite>::type make_gt_host_view<long, RANK, Intent::ReadWrite>(   \
        const Array& array );                                                                                          \
    template typename gt_view<long unsigned, RANK, Intent::ReadOnly>::type                                             \
    make_gt_host_view<long unsigned, RANK, Intent::ReadOnly>( const Array& array );                                    \
    template typename gt_view<long unsigned, RANK, Intent::ReadWrite>::type                                            \
    make_gt_host_view<long unsigned, RANK, Intent::ReadWrite>( const Array& array );                                   \
    template typename gt_view<float, RANK, Intent::ReadOnly>::type make_gt_host_view<float, RANK, Intent::ReadOnly>(   \
        const Array& array );                                                                                          \
    template typename gt_view<float, RANK, Intent::ReadWrite>::type make_gt_host_view<float, RANK, Intent::ReadWrite>( \
        const Array& array );                                                                                          \
    template typename gt_view<double, RANK, Intent::ReadOnly>::type make_gt_host_view<double, RANK, Intent::ReadOnly>( \
        const Array& array );                                                                                          \
    template typename gt_view<double, RANK, Intent::ReadWrite>::type                                                   \
    make_gt_host_view<double, RANK, Intent::ReadWrite>( const Array& array );                                          \
                                                                                                                       \
    template typename gt_view<int, RANK, Intent::ReadOnly>::type make_gt_device_view<int, RANK, Intent::ReadOnly>(     \
        const Array& array );                                                                                          \
    template typename gt_view<int, RANK, Intent::ReadWrite>::type make_gt_device_view<int, RANK, Intent::ReadWrite>(   \
        const Array& array );                                                                                          \
    template typename gt_view<long, RANK, Intent::ReadOnly>::type make_gt_device_view<long, RANK, Intent::ReadOnly>(   \
        const Array& array );                                                                                          \
    template typename gt_view<long, RANK, Intent::ReadWrite>::type make_gt_device_view<long, RANK, Intent::ReadWrite>( \
        const Array& array );                                                                                          \
    template typename gt_view<long unsigned, RANK, Intent::ReadOnly>::type                                             \
    make_gt_device_view<long unsigned, RANK, Intent::ReadOnly>( const Array& array );                                  \
    template typename gt_view<long unsigned, RANK, Intent::ReadWrite>::type                                            \
    make_gt_device_view<long unsigned, RANK, Intent::ReadWrite>( const Array& array );                                 \
    template typename gt_view<float, RANK, Intent::ReadOnly>::type make_gt_device_view<float, RANK, Intent::ReadOnly>( \
        const Array& array );                                                                                          \
    template typename gt_view<float, RANK, Intent::ReadWrite>::type                                                    \
    make_gt_device_view<float, RANK, Intent::ReadWrite>( const Array& array );                                         \
    template typename gt_view<double, RANK, Intent::ReadOnly>::type                                                    \
    make_gt_device_view<double, RANK, Intent::ReadOnly>( const Array& array );                                         \
    template typename gt_view<double, RANK, Intent::ReadWrite>::type                                                   \
    make_gt_device_view<double, RANK, Intent::ReadWrite>( const Array& array );                                        \
    }

// For each Rank in [1..9]
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
