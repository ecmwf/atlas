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

#include <iosfwd>
#include <vector>

#include "atlas/library/config.h"
#include "atlas/util/Object.h"

//----------------------------------------------------------------------------------------------------------------------

namespace atlas {
namespace mesh {
namespace detail {

class MeshImpl;

//----------------------------------------------------------------------------------------------------------------------

class PartitionGraph : public util::Object {
public:
    using Neighbours = std::vector<idx_t>;

public:
    PartitionGraph();
    PartitionGraph(idx_t values[], idx_t rows, idx_t displs[], idx_t counts[]);
    size_t footprint() const;
    idx_t size() const;
    Neighbours nearestNeighbours(const idx_t partition) const;
    idx_t maximumNearestNeighbours() const;
    operator bool() const;

private:
    void print(std::ostream&) const;
    friend std::ostream& operator<<(std::ostream& s, const PartitionGraph& p);

private:
    std::vector<idx_t> counts_;
    std::vector<idx_t> displs_;
    std::vector<idx_t> values_;
    idx_t maximum_nearest_neighbours_;
};

PartitionGraph* build_partition_graph(const MeshImpl& mesh);

//----------------------------------------------------------------------------------------------------------------------

}  // namespace detail
}  // namespace mesh
}  // namespace atlas
