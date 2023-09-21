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

#include "atlas/grid/detail/distribution/DistributionFunction.h"

namespace atlas {
namespace grid {
namespace detail {
namespace distribution {

class SerialDistribution : public DistributionFunctionT<SerialDistribution> {
public:
    SerialDistribution(const Grid& grid);
    SerialDistribution(const Grid& grid, int part);

    ATLAS_ALWAYS_INLINE int function(gidx_t gidx) const { return part_; }

private:
    int part_{0};
};

}  // namespace distribution
}  // namespace detail
}  // namespace grid
}  // namespace atlas
