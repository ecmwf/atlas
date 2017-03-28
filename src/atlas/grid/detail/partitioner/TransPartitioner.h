/*
 * (C) Copyright 1996-2017 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#pragma once

#include "atlas/grid/detail/partitioner/Partitioner.h"

namespace atlas {
namespace grid {
class StructuredGrid;
class Distribution;
}
}

namespace atlas {
namespace trans {
class Trans;
}
}


namespace atlas {
namespace grid {
namespace detail {
namespace partitioner {

class TransPartitioner: public Partitioner {

public:

    /// @brief Constructor
    /// This constructor allocates a new Trans, but without the computations
    /// of the spectral coefficients (LDGRIDONLY=TRUE)
    TransPartitioner();

    TransPartitioner(const size_t nb_partitions );

    virtual ~TransPartitioner();

    virtual void partition(const Grid&, int part[]) const;

    int nb_bands() const;

    int nb_regions(int b) const;

private:

    trans::Trans* new_trans(const Grid&) const;

private:

    size_t nbands_;
    std::vector<size_t> nregions_;
};

} // namespace partitioner
} // namespace detail
} // namespace grid
} // namespace atlas
