/*
 * (C) Copyright 1996-2017 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

/// @author Pedro Maciel
/// @author Willem Deconinck
/// @date September 2017

#pragma once

#include <vector>
#include "eckit/memory/Owned.h"
#include "atlas/library/config.h"
#include "atlas/mesh/detail/MeshImpl.h"
#include "atlas/mesh/detail/Polygon.h"
#include "atlas/util/Point.h"

namespace atlas {
namespace mesh {

/**
 * @brief Polygon class that holds the boundary of a mesh partition
 */
class PartitionPolygon : public detail::Polygon, public eckit::Owned {

public: // methods

    //-- Constructors

    /// @brief Construct "size" Polygon
    PartitionPolygon(const detail::MeshImpl& mesh, size_t halo);

    //-- Accessors

    size_t halo() const {
        return halo_;
    }

    /// @brief Return the memory footprint of the Polygon
    size_t footprint() const;

    void outputPythonScript(const eckit::PathName&) const;

    /*
     * Point-in-partition test based on winding number for a point in a polygon
     * @note reference <a href="http://geomalgorithms.com/a03-_inclusion.html">Inclusion of a Point in a Polygon</a>
     * @param[in] points vertex points of a polygon (closed, where poly.front() == poly.back())
     * @param[in] P given point
     * @return if point is in partition
     */
    bool containsPointInLonLatGeometry(const std::vector<PointLonLat>&, const PointLonLat&) const;

    /*
     * Point-in-partition test based on winding number for a point in a polygon
     * @note reference <a href="http://geomalgorithms.com/a03-_inclusion.html">Inclusion of a Point in a Polygon</a>
     * @param[in] points vertex points of a polygon (closed, where poly.front() == poly.back())
     * @param[in] P given point
     * @return if point is in partition
     */
    bool containsPointInSphericalGeometry(const std::vector<PointLonLat>&, const PointLonLat&) const;

private:

    void print(std::ostream&) const;

    friend std::ostream& operator<<(std::ostream& s, const PartitionPolygon& p) {
        p.print(s);
        return s;
    }

private:

    const detail::MeshImpl& mesh_;
    size_t halo_;

};

//------------------------------------------------------------------------------------------------------

} // namespace mesh
} // namespace atlas
