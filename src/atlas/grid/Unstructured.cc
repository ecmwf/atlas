/*
 * (C) Copyright 1996-2014 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include <limits>

#include "eckit/log/Log.h"

#include "atlas/grid/Unstructured.h"

using namespace eckit;

namespace atlas {
namespace grid {

//-----------------------------------------------------------------------------

Unstructured::Unstructured( std::vector< Point >* pts, const std::string& hash ) :
    points_(pts),
    hash_(hash),
    the_grid_spec_("unstructured","U")

{
    const std::vector<Point>& p = *points_;
    const size_t npts = p.size();

    double lat_min = std::numeric_limits<double>::max();
    double lat_max = std::numeric_limits<double>::min();
    double lon_min = lat_min;
    double lon_max = lat_max;

    for( size_t n = 0; n < npts; ++n )
    {
        lat_min = std::min( lat_min, p[n].lat() );
        lat_max = std::max( lat_max, p[n].lat() );
        lon_min = std::min( lon_min, p[n].lon() );
        lon_max = std::max( lon_max, p[n].lon() );
    }

    bound_box_ = BoundBox( Point(lat_min,lon_min), Point(lat_max,lon_max) );
}

Unstructured::~Unstructured()
{
}

std::string Unstructured::hash() const
{
    return hash_;
}

Grid::BoundBox Unstructured::boundingBox() const
{
    return bound_box_;
}

size_t Unstructured::nPoints() const
{
    return points_->size();
}

void Unstructured::coordinates( Grid::Coords& r ) const
{
    ASSERT( r.size() == points_->size() );

    const std::vector<Point>& p = *points_;
    const size_t npts = p.size();

    for( size_t i = 0; i < npts; ++i )
    {
        r.lat(i) = p[i].lat();
        r.lon(i) = p[i].lon();
    }
}

//-----------------------------------------------------------------------------

} // namespace grid
} // namespace eckit
