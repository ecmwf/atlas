/*
 * (C) Copyright 1996-2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "eckit/log/Log.h"

#include "eckit/grid/LatLon.h"

//-----------------------------------------------------------------------------

namespace eckit {
namespace grid {

//-----------------------------------------------------------------------------

LatLon::LatLon( size_t nlat, size_t nlon, const BoundBox2D& bb) :
    nlat_(nlat),
    nlon_(nlon),
    bound_box_(bb)
{
    ASSERT( nlat > 0 );
    ASSERT( nlon > 0 );

    Log::info() << "Build a LatLon with " << nlat << " x " <<  nlon << std::endl;

    points_.reserve( (nlat_ + 1) * (nlon_ + 1) );

    double dlat = ( bb.top_right_.lat_ - bb.bottom_left_.lat_ ) / nlat ;
    double dlon = ( bb.top_right_.lon_ - bb.bottom_left_.lon_ ) / nlon ;

    double plat = bb.bottom_left_.lat_;
    for( size_t i = 0; i <= nlat_; ++i )
    {
        double plon = bb.bottom_left_.lon_;
        for( size_t j = 0; j <= nlon_; ++j )
        {
            points_.push_back( Point2D(  plat,plon ) );
            plon += dlon;
        }
        plat += dlat;
    }
}

LatLon::~LatLon()
{
    Log::info() << "Destroy a LatLon" << std::endl;
}

BoundBox2D LatLon::boundingBox() const
{
    return bound_box_;
}

//-----------------------------------------------------------------------------

} // namespace grid
} // namespace eckit
