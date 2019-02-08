/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#pragma once

#include "gridtools/common/generic_metafunctions/accumulate.hpp"
#include "gridtools/common/generic_metafunctions/is_all_integrals.hpp"
#include "gridtools/storage/storage-facility.hpp"

#include "atlas/array/ArrayViewDefs.h"
#include "atlas/library/config.h"

//------------------------------------------------------------------------------

namespace atlas {
namespace array {
namespace gridtools {

//------------------------------------------------------------------------------

#if ATLAS_GRIDTOOLS_STORAGE_BACKEND_CUDA
using storage_traits = ::gridtools::storage_traits<::gridtools::target::cuda>;
#elif ATLAS_GRIDTOOLS_STORAGE_BACKEND_HOST
using storage_traits = ::gridtools::storage_traits<::gridtools::target::x86>;
#else
#error ATLAS_GRIDTOOLS_STORAGE_BACKEND_<HOST,CUDA> not set
#endif

//------------------------------------------------------------------------------

template <typename Value, unsigned int Rank, ::gridtools::access_mode AccessMode = ::gridtools::access_mode::ReadWrite>
using data_view_tt = ::gridtools::data_view<
    gridtools::storage_traits::data_store_t<Value, gridtools::storage_traits::storage_info_t<0, Rank>>, AccessMode>;

inline constexpr ::gridtools::access_mode get_access_mode( Intent kind ) {
    return ( kind == Intent::ReadOnly ) ? ::gridtools::access_mode::ReadOnly : ::gridtools::access_mode::ReadWrite;
}

//------------------------------------------------------------------------------

}  // namespace gridtools
}  // namespace array
}  // namespace atlas
