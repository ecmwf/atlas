/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

// TL63

#include "N.h"

namespace atlas {
namespace grid {
namespace detail {
namespace pl {
namespace classic_gaussian {

DEFINE_POINTS_PER_LATITUDE(32, LIST(20, 27, 36, 40, 45, 50, 60, 64, 72, 75, 80, 90, 90, 96, 100, 108, 108, 120, 120,
                                    120, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128))

}  // namespace classic_gaussian
}  // namespace pl
}  // namespace detail
}  // namespace grid
}  // namespace atlas
