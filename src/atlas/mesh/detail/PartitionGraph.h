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

#include "eckit/memory/Owned.h"
#include "eckit/memory/SharedPtr.h"

//----------------------------------------------------------------------------------------------------------------------

namespace atlas {
namespace mesh {
namespace detail {

class MeshImpl;

//----------------------------------------------------------------------------------------------------------------------

class PartitionGraph : public eckit::Owned {
public:
    using Neighbours = std::vector<size_t>;

public:
    PartitionGraph();
    PartitionGraph( size_t values[], size_t rows, size_t displs[], size_t counts[] );
    size_t footprint() const;
    size_t size() const;
    Neighbours nearestNeighbours( const size_t partition ) const;
    size_t maximumNearestNeighbours() const;
    operator bool() const;

private:
    void print( std::ostream& ) const;
    friend std::ostream& operator<<( std::ostream& s, const PartitionGraph& p );

private:
    std::vector<size_t> counts_;
    std::vector<size_t> displs_;
    std::vector<size_t> values_;
    size_t maximum_nearest_neighbours_;
};

PartitionGraph* build_partition_graph( const MeshImpl& mesh );

//----------------------------------------------------------------------------------------------------------------------

}  // namespace detail
}  // namespace mesh
}  // namespace atlas
