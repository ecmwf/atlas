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

#ifndef atlas_grid_Gaussian_H
#define atlas_grid_Gaussian_H

#include <cstddef>
#include <vector>

#include "eckit/types/Coord.h"
#include "atlas/grid/Grid.h"

//-----------------------------------------------------------------------------

namespace atlas {
namespace grid {

//-----------------------------------------------------------------------------

class Gaussian : public Grid {

public: // methods

    Gaussian( size_t resolution, const BoundBox2D& bb );

    virtual ~Gaussian();

    virtual const std::vector<Point2D>& coordinates() const { return coordinates_; }
    virtual BoundBox2D boundingBox() const;

protected:

    size_t resolution_;                 ///< number of longitude increments - can be any size as no requirement for 

    std::vector< Point2D > coordinates_;     ///< storage of coordinate points

    BoundBox2D bound_box_;              ///< bounding box for the domain

private:

};

//-----------------------------------------------------------------------------

} // namespace grid
} // namespace eckit

#endif
