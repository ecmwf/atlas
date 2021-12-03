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
#include <vector>


#include "atlas/grid/detail/distribution/DistributionImpl.h"
#include "atlas/util/vector.h"


namespace atlas {
namespace grid {
class Partitioner;

namespace detail {
namespace distribution {

class DistributionArray : public DistributionImpl {
public:
    using partition_t = atlas::vector<int>;

    DistributionArray() = default;

    DistributionArray(const Grid&, const eckit::Parametrisation&);

    DistributionArray(const Grid&, const Partitioner&);

    DistributionArray(int nb_partitions, idx_t npts, int partition[], int part0 = 0);

    DistributionArray(int nb_partitions, partition_t&& partition);

    virtual ~DistributionArray();

    int partition(const gidx_t gidx) const override { return part_[gidx]; }

    const partition_t& partition() const { return part_; }

    idx_t nb_partitions() const override { return nb_partitions_; }

    operator const partition_t&() const { return part_; }

    const int* data() const { return part_.data(); }

    const std::vector<idx_t>& nb_pts() const override { return nb_pts_; }

    idx_t max_pts() const override { return max_pts_; }
    idx_t min_pts() const override { return min_pts_; }

    const std::string& type() const override { return type_; }

    void print(std::ostream&) const override;

    size_t footprint() const override { return nb_pts_.size() * sizeof(nb_pts_[0]) + part_.size() * sizeof(part_[0]); }

    bool functional() const override { return false; }

    gidx_t size() const override { return part_.size(); }

    void hash(eckit::Hash&) const override;

    void partition(gidx_t begin, gidx_t end, int partitions[]) const override {
        size_t i = 0;
        for (gidx_t n = begin; n < end; ++n, ++i) {
            partitions[i] = part_[n];
        }
    }

protected:
    idx_t nb_partitions_ = 0;

    partition_t part_;
    std::vector<idx_t> nb_pts_;
    idx_t max_pts_;
    idx_t min_pts_;
    std::string type_;
};

}  // namespace distribution
}  // namespace detail
}  // namespace grid
}  // namespace atlas
