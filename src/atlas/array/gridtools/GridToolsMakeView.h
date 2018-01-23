/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */
#pragma once

#include "atlas/array_fwd.h"
#include "atlas/array/gridtools/GridToolsTraits.h"

namespace atlas {
namespace array {
namespace gridtools {

template <typename Value, unsigned int Rank, Intent AccessMode = Intent::ReadWrite>
data_view_tt<Value, Rank, get_access_mode(AccessMode) >
make_gt_host_view(const Array& array);

template <typename Value, unsigned int Rank, Intent AccessMode = Intent::ReadWrite>
data_view_tt<Value, Rank, get_access_mode(AccessMode) >
make_gt_device_view(const Array& array);

} // namespace gridtools
} // namespace array
} // namespace atlas
