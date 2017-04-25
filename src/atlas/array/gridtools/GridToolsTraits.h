#pragma once

#include "atlas/library/config.h"

//------------------------------------------------------------------------------
#if ATLAS_GRIDTOOLS_STORAGE_BACKEND_CUDA
#define ENABLE_GPU
#endif
#include "storage-facility.hpp"
#ifdef ENABLE_GPU
#undef ENABLE_GPU
#endif
//------------------------------------------------------------------------------

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

template <typename Value, unsigned int Rank, bool ReadOnly = false>
using data_view_tt = ::gridtools::data_view<
        gridtools::storage_traits::data_store_t<
          Value,
          gridtools::storage_traits::storage_info_t<0, Rank> >,
          ReadOnly>;

//------------------------------------------------------------------------------

} // namespace gridtools
} // namespace array
} // namespace atlas
