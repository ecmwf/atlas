/*
 * (C) Copyright 1996-2017 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "atlas/mesh/detail/Polygon.h"

#include <algorithm>
#include <iostream>
#include <limits>
#include "eckit/exception/Exceptions.h"
#include "eckit/types/FloatCompare.h"
#include "atlas/util/CoordinateEnums.h"

namespace atlas {
namespace mesh {
namespace detail {

//------------------------------------------------------------------------------------------------------


Polygon::Polygon() {
}


Polygon::Polygon(const Polygon::edge_set_t& edges) {
    // get external edges by attempting to remove reversed edges, if any
    edge_set_t extEdges;
    for (const edge_t& e : edges) {
        if (!extEdges.erase(e.reverse())) {
            extEdges.insert(e);
        }
    }
    ASSERT(extEdges.size() >= 2);

    // set one polygon cycle, by picking next edge with first node same as second node of last edge
    clear();
    reserve(extEdges.size() + 1);

    push_back(extEdges.begin()->first);
    for (edge_set_t::iterator e = extEdges.begin(); e != extEdges.end() && e->first == back();
         e = extEdges.lower_bound(edge_t(back(), std::numeric_limits< idx_t >::min()))) {
        push_back(e->second);
        extEdges.erase(*e);
    }
    ASSERT(front() == back());

    // exhaust remaining edges, should represent additional cycles, if any
    while (!extEdges.empty()) {
        operator+=(Polygon(extEdges));
    }
}


Polygon::operator bool() const {
    return !Polygon::empty();
}


Polygon& Polygon::operator+=(const Polygon& other) {

    if (empty()) {
        return operator=(other);
    }

    // polygon can have multiple cycles, but must be connected graphs
    // Note: a 'cycle' is handled by repeating the indices, excluding (repeated) last index
    ASSERT(other.front() == other.back());
    const difference_type N = difference_type(other.size()) - 1;

    container_t cycle(2 * size_t(N));
    std::copy(other.begin(), other.begin() + N, cycle.begin());
    std::copy(other.begin(), other.begin() + N, cycle.begin() + N);

    for (const_iterator c = cycle.begin(); c != cycle.begin() + N; ++c) {
        iterator here = std::find(begin(), end(), *c);
        if (here != end()) {
            insert(here, c, c + N);
            return *this;
        }
    }

    throw eckit::AssertionFailed("Polygon: could not merge polygons, they are not connected", Here());
}


void Polygon::print(std::ostream& s) const {
    char z = '{';
    for (auto n : static_cast<const container_t&>(*this)) {
        s << z << n;
        z = ',';
    }
    s << '}';
}


//------------------------------------------------------------------------------------------------------

}  // detail
}  // mesh
}  // atlas

