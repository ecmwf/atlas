#pragma once

#include "atlas/library/config.h"

//------------------------------------------------------------------------------
#if ATLAS_GRIDTOOLS_STORAGE_BACKEND_CUDA
#ifndef _USE_GPU_
#define _USE_GPU_
#endif
#endif
#include "common/generic_metafunctions/all_integrals.hpp"
#include "storage/storage-facility.hpp"
#include "atlas/array/ArrayViewDefs.h"

#ifdef _USE_GPU_
#undef _USE_GPU_
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

template <typename Value, unsigned int Rank, ::gridtools::access_mode AccessMode = ::gridtools::access_mode::ReadWrite >
using data_view_tt = ::gridtools::data_view<
        gridtools::storage_traits::data_store_t<
          Value,
          gridtools::storage_traits::storage_info_t<0, Rank> >,
          AccessMode>;

inline constexpr ::gridtools::access_mode get_access_mode(Intent kind) {
   return (kind == Intent::ReadOnly) ? ::gridtools::access_mode::ReadOnly : ::gridtools::access_mode::ReadWrite;
}

//------------------------------------------------------------------------------

} // namespace gridtools
} // namespace array
} // namespace atlas
