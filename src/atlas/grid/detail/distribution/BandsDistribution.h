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


class BandsDistribution : public DistributionFunctionT<BandsDistribution> {
private:
    size_t blocksize_;
    idx_t nb_blocks_;

public:
    BandsDistribution( const Grid& grid, idx_t nb_partitions, const std::string& type, size_t blocksize = 1 );

    int function( gidx_t gidx ) const {
        idx_t iblock = gidx / blocksize_;
        return ( iblock * nb_partitions_ ) / nb_blocks_;
    }
};


}  // namespace distribution
}  // namespace detail
}  // namespace grid
}  // namespace atlas
