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

#include "atlas/grid/detail/distribution/DistributionImpl.h"
#include "atlas/library/config.h"
#include "atlas/util/Config.h"
#include "atlas/util/ObjectHandle.h"
#include "atlas/util/vector.h"

namespace atlas {
class Grid;
namespace grid {
class Partitioner;
}  // namespace grid
}  // namespace atlas

namespace atlas {
namespace grid {

class Distribution : DOXYGEN_HIDE(public util::ObjectHandle<DistributionImpl>) {
    friend class Partitioner;

public:
    using Config      = DistributionImpl::Config;
    using partition_t = atlas::vector<int>;

    using Handle::Handle;
    Distribution() = default;

    /// @brief Create a serial distribution
    Distribution(const Grid&);

    /// @brief Create a distribution specified by a configuration
    Distribution(const Grid&, const Config&);

    /// @brief Create a distribution using a given partitioner
    Distribution(const Grid&, const Partitioner&);

    /// @brief Create a distribution by given array, and make internal copy
    Distribution(int nb_partitions, idx_t npts, int partition[], int part0 = 0);

    /// @brief Create a distribution by given array, and take ownership (move)
    Distribution(int nb_partitions, partition_t&& partition);

    ~Distribution();

    ATLAS_ALWAYS_INLINE int partition(gidx_t index) const { return get()->partition(index); }

    template <typename PartitionContainer>
    void partition(gidx_t begin, gidx_t end, PartitionContainer& partitions) const {
        ATLAS_ASSERT(end - begin <= partitions.size());
        return get()->partition(begin, end, partitions.data());
    }

    size_t footprint() const { return get()->footprint(); }

    ATLAS_ALWAYS_INLINE idx_t nb_partitions() const { return get()->nb_partitions(); }

    ATLAS_ALWAYS_INLINE gidx_t size() const { return get()->size(); }

    const std::vector<idx_t>& nb_pts() const;

    idx_t max_pts() const;

    idx_t min_pts() const;

    const std::string& type() const;

    friend std::ostream& operator<<(std::ostream& os, const Distribution& distribution);

    std::string hash() const;

    void hash(eckit::Hash&) const;
};

}  // namespace grid
}  // namespace atlas
