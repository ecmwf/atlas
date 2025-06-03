/*
 * (C) Crown copyright 2025 Met Office.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include "atlas/grid/Distribution.h"
#include "atlas/grid/Grid.h"
#include "atlas/grid/detail/partitioner/Partitioner.h"

namespace atlas {
namespace grid {
namespace detail {
namespace partitioner {

class GlobalIndexPartitioner : public Partitioner {
public:
    GlobalIndexPartitioner();
    GlobalIndexPartitioner(int N);
    GlobalIndexPartitioner(int N, const eckit::Parametrisation& config);
    GlobalIndexPartitioner(const eckit::Parametrisation& config);

    using Partitioner::partition;
    void partition(const Grid& grid, int part[]) const override;
    Distribution partition(const Grid& grid) const override;

    std::string type() const override { return "global_index"; }
};

}  // namespace partitioner
}  // namespace detail
}  // namespace grid
}  // namespace atlas
