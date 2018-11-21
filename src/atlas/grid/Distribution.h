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

#include <vector>

#include "eckit/memory/SharedPtr.h"

#include "atlas/grid/detail/distribution/DistributionImpl.h"
#include "atlas/library/config.h"

namespace atlas {
class Grid;
namespace grid {
class Partitioner;
class DistributionImpl;
class Partitioner;
}  // namespace grid
}  // namespace atlas

namespace atlas {
namespace grid {

class Distribution {
    friend class Partitioner;

public:
    using Implementation = DistributionImpl;

public:
    Distribution();
    Distribution( const Implementation* );
    Distribution( const Distribution& );

    Distribution( const Grid& );

    Distribution( const Grid&, const Partitioner& );

    Distribution( idx_t npts, int partition[], int part0 = 0 );

    ~Distribution();

    int partition( const gidx_t gidx ) const;

    const std::vector<int>& partition() const;

    idx_t nb_partitions() const;

    operator const std::vector<int>&() const;

    const int* data() const;

    const std::vector<idx_t>& nb_pts() const;

    idx_t max_pts() const;
    idx_t min_pts() const;

    const std::string& type() const;

    friend std::ostream& operator<<( std::ostream& os, const Distribution& distribution );

    const Implementation* get() const;

private:
    eckit::SharedPtr<const Implementation> impl_;
};

}  // namespace grid
}  // namespace atlas
