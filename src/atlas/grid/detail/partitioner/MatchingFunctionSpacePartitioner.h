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

#include "atlas/functionspace/FunctionSpace.h"
#include "atlas/grid/detail/partitioner/Partitioner.h"

namespace atlas {
namespace grid {
namespace detail {
namespace partitioner {

class MatchingFunctionSpacePartitioner : public Partitioner {
public:
    MatchingFunctionSpacePartitioner();

    // MatchingFunctionSpacePartitioner(const idx_t nb_partitions);

    MatchingFunctionSpacePartitioner(const FunctionSpace&, const eckit::Parametrisation&);

    virtual ~MatchingFunctionSpacePartitioner() override {}

protected:
    const FunctionSpace partitioned_;
};

}  // namespace partitioner
}  // namespace detail
}  // namespace grid
}  // namespace atlas
