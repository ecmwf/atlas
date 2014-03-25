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

#ifndef eckit_grid_Grid_H
#define eckit_grid_Grid_H

#include <cstddef>
#include <vector>
#include <cmath>

#include "eckit/memory/NonCopyable.h"
#include "eckit/exception/Exceptions.h"

//-----------------------------------------------------------------------------

namespace eckit {
namespace grid {

//-----------------------------------------------------------------------------

class Point2D {
public:

    Point2D( double lat, double lon ) : lat_(lat),lon_(lon) {}
    Point2D( ) : lat_(0.0), lon_(0.0) {}

    double lat_;
    double lon_;

    /// Calculates distance between two points
    static double distance(const Point2D& a, const Point2D& b)
    {
        return std::sqrt(distance2(a, b));
    }
    
    /// Calculates distance squared between two points
    static double distance2(const Point2D& a, const Point2D& b)
    {
        const double dlat = a.lat_ - b.lat_;
        const double dlon = a.lon_ - b.lon_;
        return dlat*dlat + dlon*dlon;
    }

    /// Tests whether two points can be considered equal
    /// @todo the epsilon could be imported from config file
    static bool equal(const Point2D& a, const Point2D& b, const double epsilon = 1.0e-8 )
    {
		/// @todo take epsilon from some general config
		
        return ((std::fabs(a.lat_ - b.lat_) < epsilon ) &&
                (std::fabs(a.lon_ - b.lon_) < epsilon ) );
                
    }
};

//-----------------------------------------------------------------------------

struct BoundBox2D
{
    BoundBox2D( const Point2D& bottom_left, const Point2D& top_right ) :
        bottom_left_(bottom_left),
        top_right_(top_right)
    {
        ASSERT( bottom_left_.lat_ < top_right_.lat_ );
        ASSERT( bottom_left_.lon_ < top_right_.lon_ );
    }

    Point2D bottom_left_;
    Point2D top_right_;
};

//-----------------------------------------------------------------------------

class Grid : private eckit::NonCopyable {

public: // methods

    Grid();

    virtual ~Grid();

    virtual const std::vector<Point2D>& coordinates() const = 0;

    virtual BoundBox2D boundingBox() const = 0;

    virtual std::string hash() const;

protected:

};

//-----------------------------------------------------------------------------

} // namespace grid
} // namespace eckit

#endif
