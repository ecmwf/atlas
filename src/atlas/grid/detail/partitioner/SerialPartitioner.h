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

#include "atlas/grid/Distribution.h"
#include "atlas/grid/Grid.h"
#include "atlas/grid/detail/distribution/SerialDistribution.h"
#include "atlas/grid/detail/partitioner/Partitioner.h"

namespace atlas {
namespace grid {
namespace detail {
namespace partitioner {

class SerialPartitioner : public Partitioner {
public:
    SerialPartitioner() : Partitioner( 1 ) {}
    SerialPartitioner( int N, const eckit::Parametrisation& config ) : Partitioner( 1 ) {}

    SerialPartitioner( int N ) : Partitioner( N ) {}

    std::string type() const override { return static_type(); }
    static std::string static_type() { return "serial"; }

    virtual Distribution partition( const Grid& grid ) const override {
        return Distribution{new distribution::SerialDistribution{grid}};
    }

    virtual void partition( const Grid& grid, int part[] ) const {
        gidx_t gridsize = grid.size();
        for ( gidx_t n = 0; n < gridsize; ++n ) {
            part[n] = 0.;
        }
    }
};


}  // namespace partitioner
}  // namespace detail
}  // namespace grid
}  // namespace atlas
