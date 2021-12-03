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

#include <string>

#include "atlas/grid/detail/distribution/DistributionFunction.h"

namespace atlas {
namespace grid {
namespace detail {
namespace distribution {

template <typename Int>
class BandsDistribution : public DistributionFunctionT<BandsDistribution<Int>> {
private:
    Int blocksize_;
    Int nb_blocks_;
    Int nb_partitions_Int_;

public:
    BandsDistribution(const Grid& grid, idx_t nb_partitions, const std::string& type, size_t blocksize = 1);

    ATLAS_ALWAYS_INLINE int function(gidx_t index) const {
        Int iblock = static_cast<Int>(index / blocksize_);
        return (iblock * nb_partitions_Int_) / nb_blocks_;
    }

    static bool detectOverflow(size_t gridsize, size_t nb_partitions, size_t blocksize);
};


}  // namespace distribution
}  // namespace detail
}  // namespace grid
}  // namespace atlas
