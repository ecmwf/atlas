/*
 * (C) Copyright 1996-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */


#pragma once

#include <vector>
#include "eckit/memory/Owned.h"
#include "eckit/memory/SharedPtr.h"
#include "atlas/internals/atlas_config.h"

namespace atlas {
namespace grid {
class Grid;
}
}
namespace atlas {
namespace grid {
class Partitioner;
}
}

namespace atlas {
namespace grid {

class GridDistribution: public eckit::Owned {

public:

     typedef eckit::SharedPtr<GridDistribution> Ptr;

public:

    GridDistribution(const Grid&);

    GridDistribution(const Partitioner&);

    GridDistribution(size_t npts, int partition[], int part0 = 0);

    virtual ~GridDistribution() {}

    int partition(const gidx_t gidx) const {
        return part_[gidx];
    }

    const std::vector<int>& partition() const {
        return part_;
    }

    size_t nb_partitions() const {
        return nb_partitions_;
    }

    operator const std::vector<int>&() const {
        return part_;
    }

    const int* data() const {
        return part_.data();
    }

    const std::vector<int>& nb_pts() const {
        return nb_pts_;
    }

    size_t max_pts() const {
        return max_pts_;
    }
    size_t min_pts() const {
        return min_pts_;
    }

private:

    size_t nb_partitions_;
    std::vector<int> part_;
    std::vector<int> nb_pts_;
    size_t max_pts_;
    size_t min_pts_;
};

extern "C" {
    GridDistribution* atlas__GridDistribution__new(int npts, int part[], int part0);
    void atlas__GridDistribution__delete(GridDistribution* This);
}

} // namespace grid
} // namespace atlas
