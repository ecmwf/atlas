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
#include <array>

#include "atlas/grid/detail/partitioner/Partitioner.h"

namespace atlas {
namespace grid {
namespace detail {
namespace partitioner {

class CheckerboardPartitioner : public Partitioner {
public:
    CheckerboardPartitioner();

    // CheckerboardPartitioner(int N);  // N is the number of parts (aka MPI tasks)
    CheckerboardPartitioner(int N, const eckit::Parametrisation&);

    CheckerboardPartitioner(int N, int nbands);
    CheckerboardPartitioner(int N, int nbands, bool checkerboard);

    CheckerboardPartitioner(const eckit::Parametrisation&);

    // Node struct that holds the x and y indices (for global, it's longitude and
    // latitude in millidegrees (integers))
    // This structure is used in sorting algorithms, and uses less memory than
    // if x and y were in double precision.
    struct NodeInt {
        int x, y;
        int n;
    };

    virtual std::string type() const { return "checkerboard"; }

    std::array<int,2> checkerboardDimensions(const Grid&);

private:
    struct Checkerboard {
        idx_t nbands;  // number of bands
        idx_t nx, ny;  // grid dimensions
    };

    Checkerboard checkerboard(const Grid&) const;

    // Doesn't matter if nodes[] is in degrees or radians, as a sorting
    // algorithm is used internally
    void partition(const Checkerboard& cb, int nb_nodes, NodeInt nodes[], int part[]) const;

    using Partitioner::partition;
    virtual void partition(const Grid&, int part[]) const;

    void check() const;

private:
    idx_t nbands_      = 0;  // number of bands from configuration
    bool regular_      = false;
    bool split_x_      = true;
    bool split_y_      = true;
    bool checkerboard_ = true;  // exact (true) or approximate (false) checkerboard
};

}  // namespace partitioner
}  // namespace detail
}  // namespace grid
}  // namespace atlas
