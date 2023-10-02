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

#include "atlas/grid/detail/partitioner/Partitioner.h"

#include "atlas/grid/Distribution.h"
#include "atlas/grid/Grid.h"

namespace atlas {
namespace grid {
namespace detail {
namespace partitioner {

class BandsPartitioner : public Partitioner {
private:
    int blocksize_;
    static size_t extract_blocksize(const eckit::Parametrisation& config) {
        size_t blocksize{1};
        config.get("blocksize", blocksize);
        return blocksize;
    }

    size_t blocksize(const Grid& grid) const;

protected:
    static constexpr int BLOCKSIZE_NX{-1};

public:
    BandsPartitioner(const eckit::Parametrisation& config = util::NoConfig());
    BandsPartitioner(int N, const eckit::Parametrisation& config): BandsPartitioner(N, extract_blocksize(config), config) {}
    BandsPartitioner(int N, int blocksize, const eckit::Parametrisation& config = util::NoConfig());

    std::string type() const override { return static_type(); }
    static std::string static_type() { return "bands"; }

    Distribution partition(const Grid& grid) const override;

    void partition(const Grid& grid, int part[]) const override;
};

}  // namespace partitioner
}  // namespace detail
}  // namespace grid
}  // namespace atlas
