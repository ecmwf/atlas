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

#include "atlas/domain/Domain.h"
#include "atlas/functionspace/FunctionSpace.h"
#include "atlas/library/config.h"
#include "atlas/projection/Projection.h"
#include "atlas/util/Config.h"
#include "atlas/util/Point.h"
#include "atlas/util/Polygon.h"

namespace atlas {
namespace grid {

// -------------------------------------------------------------------

/**
 * @brief StructuredPartitionPolygon class that holds the boundary of a structured grid partition
 */
class StructuredPartitionPolygon : public util::PartitionPolygon {
public:  // methods
    //-- Constructors

    /// @brief Construct "size" polygon
    StructuredPartitionPolygon(const functionspace::FunctionSpaceImpl& fs, idx_t halo);

    //-- Accessors

    idx_t halo() const override { return halo_; }

    /// @brief Return the memory footprint of the Polygon
    size_t footprint() const override;

    void outputPythonScript(const eckit::PathName&, const eckit::Configuration& = util::NoConfig()) const override;

    PointsXY xy() const override;
    PointsLonLat lonlat() const override;

    const RectangleXY& inscribedDomain() const override { return inscribed_domain_; }

    void allGather(util::PartitionPolygons&) const override;

private:
    // void print( std::ostream& ) const;

    // friend std::ostream& operator<<( std::ostream& s, const StructuredPartitionPolygon& p ) {
    //     p.print( s );
    //     return s;
    // }

    //util::Polygon::edge_set_t compute_edges( std::vector<Point2>&, std::vector<Point2>&  );


private:
    PointsXY points_;
    PointsXY inner_bounding_box_;
    RectangleXY inscribed_domain_;
    const functionspace::FunctionSpaceImpl& fs_;
    idx_t halo_;
};

// -------------------------------------------------------------------


}  // namespace grid
}  // namespace atlas
