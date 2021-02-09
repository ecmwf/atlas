
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

#include "atlas/grid/detail/partitioner/Partitioner.h"

namespace atlas {
namespace grid {
namespace detail {
namespace partitioner {

class ZonalBoardPartitioner : public Partitioner {
public:
    ZonalBoardPartitioner();

    ZonalBoardPartitioner( int N );  // N is the number of parts (aka MPI tasks)
    ZonalBoardPartitioner( int N, const eckit::Parametrisation& );

    ZonalBoardPartitioner( int N, int nbands );
    ZonalBoardPartitioner( int N, int nbands, bool zonalboard );

    // Node struct that holds the x and y indices (for global, it's longitude and
    // latitude in millidegrees (integers))
    // This structure is used in sorting algorithms, and uses less memory than
    // if x and y were in double precision.
    struct NodeInt {
        int x, y;
        int n;
    };

    virtual std::string type() const { return "zonalboard"; }

private:
    struct Zonalboard {
        idx_t nbands;  // number of bands
        idx_t nx, ny;  // grid dimensions
    };

    Zonalboard zonalboard( const Grid& ) const;

    // Doesn't matter if nodes[] is in degrees or radians, as a sorting
    // algorithm is used internally
    void partition( const Zonalboard& cb, int nb_nodes, NodeInt nodes[], int part[] ) const;

    virtual void partition( const Grid&, int part[] ) const;

    void check() const;

private:
    idx_t nbands_;       // number of bands from configuration
    bool zonalboard_;  // exact (true) or approximate (false) zonalboard
};

}  // namespace partitioner
}  // namespace detail
}  // namespace grid
}  // namespace atlas
