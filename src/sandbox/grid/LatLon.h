/*
 * (C) Copyright 1996-2015 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

/// @author Peter Bispham
/// @author Tiago Quintino
/// @date Oct 2013

#ifndef atlas_grid_LatLon_H
#define atlas_grid_LatLon_H

#include <cstddef>
#include <vector>

#include "atlas/grid/Grid.h"
#include "atlas/grid/GridFactory.h"

//-----------------------------------------------------------------------------

namespace atlas {
namespace grid {

//-----------------------------------------------------------------------------

class LatLon : public Grid {
   REGISTER(LatLon);

public: // methods

    LatLon() : nlat_(0), nlon_(0) {}

    LatLon( size_t nlat, size_t nlon, const BoundBox& bb );

    virtual ~LatLon();

    virtual std::string hash() const;

    virtual BoundBox boundingBox() const;

    virtual size_t nPoints() const { return points_.size(); }

    virtual void coordinates( Grid::Coords & ) const;

    virtual std::string gridType() const { return std::string("lonlat") ;}

    virtual GridSpec* spec() const;

    virtual void constructFrom(const GridSpec& );

    virtual bool compare(const Grid&) const;

    /// @deprecated will be removed soon as it exposes the inner storage of the coordinates
    virtual const std::vector<Point>& coordinates() const { return points_; }

protected:

    size_t nlat_;                     ///< number of latitude  increments - ODD number for coindidence with 0,0 on Earth
    size_t nlon_;                     ///< number of longitude increments - can be any size as no requirement for

    std::vector< Point > points_;     ///< storage of coordinate points

    BoundBox bound_box_;              ///< bounding box for the domain
};

//-----------------------------------------------------------------------------

} // namespace grid
} // namespace eckit

#endif
