/*
 * (C) Copyright 2023 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#pragma once

#include <vector>

#include "atlas/grid/detail/partitioner/EqualRegionsPartitioner.h"
#include "atlas/grid/detail/partitioner/Partitioner.h"

namespace atlas {
namespace grid {
namespace detail {
namespace partitioner {

class EqualAreaPartitioner : public Partitioner {
public:
    EqualAreaPartitioner();
    EqualAreaPartitioner(int N);
    EqualAreaPartitioner(int N, const eckit::Parametrisation& config);
    EqualAreaPartitioner(const eckit::Parametrisation& config);

    using Partitioner::partition;
    void partition(const Grid&, int part[]) const override;

    std::string type() const override { return "equal_area"; }

private:
    EqualRegionsPartitioner partitioner_;
};


}  // namespace partitioner
}  // namespace detail
}  // namespace grid
}  // namespace atlas
