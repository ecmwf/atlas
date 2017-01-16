#pragma once

#include "atlas/internals/atlas_defines.h"
#include "storage-facility.hpp"


namespace atlas {
namespace array {
namespace gridtools {

//------------------------------------------------------------------------------

#if ATLAS_GRIDTOOLS_STORAGE_BACKEND_CUDA
using storage_traits = ::gridtools::storage_traits< ::gridtools::enumtype::Cuda >;
#elif ATLAS_GRIDTOOLS_STORAGE_BACKEND_HOST
using storage_traits = ::gridtools::storage_traits< ::gridtools::enumtype::Host >;
#else
#error ATLAS_GRIDTOOLS_STORAGE_BACKEND_<HOST,CUDA> not set
#endif

//------------------------------------------------------------------------------

} // namespace gridtools

//------------------------------------------------------------------------------

template <typename Value, unsigned int NDims, bool ReadOnly = false>
using data_view_tt = ::gridtools::data_view<
        gridtools::storage_traits::data_store_t<
          Value,
          gridtools::storage_traits::storage_info_t<0, NDims> >,
          ReadOnly>;

//------------------------------------------------------------------------------

} // namespace array
} // namespace atlas
