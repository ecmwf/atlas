/*
 * (C) Copyright 1996-2013 ECMWF.
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

#ifndef eckit_grid_LatLon_H
#define eckit_grid_LatLon_H

#include <cstddef>

#include "eckit/types/Coord.h"
#include "eckit/grid/Grid.h"

//-----------------------------------------------------------------------------

namespace eckit {
namespace grid {

//-----------------------------------------------------------------------------

class LatLon : public Grid {

public: // methods

    LatLon( size_t nlat, size_t nlon, const BoundBox2D& bb );

    virtual ~LatLon();

     virtual const std::vector< Point2D >& coordinates() const { return points_; }
     virtual BoundBox2D boundingBox() const;

protected:

    size_t nlat_;                       ///< number of latitude  increments - ODD number for coindidence with 0,0 on Earth 
    size_t nlon_;                       ///< number of longitude increments - can be any size as no requirement for 

    std::vector< Point2D > points_;     ///< storage of coordinate points

    BoundBox2D bound_box_;              ///< bounding box for the domain

};

//-----------------------------------------------------------------------------

} // namespace grid
} // namespace eckit

#endif
