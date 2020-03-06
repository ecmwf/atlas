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

#include <type_traits>
#include <utility>

#include "gridtools/common/generic_metafunctions/accumulate.hpp"
#include "gridtools/storage/storage_facility.hpp"

#include "atlas/array/ArrayViewDefs.h"
#include "atlas/library/config.h"

//------------------------------------------------------------------------------

namespace atlas {
namespace array {
namespace gridtools {

//------------------------------------------------------------------------------

#if ATLAS_GRIDTOOLS_STORAGE_BACKEND_CUDA
using backend_t      = ::gridtools::backend::cuda;
using storage_traits = ::gridtools::storage_traits<backend_t>;
#elif ATLAS_GRIDTOOLS_STORAGE_BACKEND_HOST
using backend_t      = ::gridtools::backend::x86;
using storage_traits = ::gridtools::storage_traits<backend_t>;
#else
#error ATLAS_GRIDTOOLS_STORAGE_BACKEND_<HOST,CUDA> not set
#endif

//------------------------------------------------------------------------------

template <typename Value, unsigned int Rank, ::gridtools::access_mode AccessMode = ::gridtools::access_mode::read_write>
using data_view_tt =
    ::gridtools::data_view<gridtools::storage_traits::data_store_t<typename std::remove_const<Value>::type,
                                                                   gridtools::storage_traits::storage_info_t<0, Rank>>,
                           AccessMode>;

template <typename Value, typename std::enable_if<std::is_const<Value>::value, Value>::type* = nullptr>
inline constexpr ::gridtools::access_mode get_access_mode() {
    return ::gridtools::access_mode::read_only;
}
template <typename Value, typename std::enable_if<!std::is_const<Value>::value, Value>::type* = nullptr>
inline constexpr ::gridtools::access_mode get_access_mode() {
    return ::gridtools::access_mode::read_write;
}

//------------------------------------------------------------------------------

}  // namespace gridtools
}  // namespace array
}  // namespace atlas
