/*
 * (C) Copyright 1996-2014 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "eckit/log/Log.h"
#include "atlas/grid/GridSpec.h"

#include "atlas/grid/LatLon.h"

using namespace eckit;

namespace atlas {
namespace grid {

//-----------------------------------------------------------------------------

LatLon::LatLon( size_t nlat, size_t nlon, const BoundBox& bb) :
    nlat_(nlat),
    nlon_(nlon),
    bound_box_(bb)
{
    ASSERT( nlat > 0 );
    ASSERT( nlon > 0 );

//    Log::info() << "Build a LatLon with " << nlat << " x " <<  nlon << std::endl;

    points_.reserve( (nlat_ + 1) * (nlon_ + 1) );

    double dlat = ( bb.top_right_.lat() - bb.bottom_left_.lat() ) / nlat ;
    double dlon = ( bb.top_right_.lon() - bb.bottom_left_.lon() ) / nlon ;

    double plat = bb.bottom_left_.lat();
    for( size_t i = 0; i <= nlat_; ++i )
    {
        double plon = bb.bottom_left_.lon();
        for( size_t j = 0; j <= nlon_; ++j )
        {
            points_.push_back( Point( plat, plon ) );
            plon += dlon;
        }
        plat += dlat;
    }
}

LatLon::~LatLon()
{
    Log::info() << "Destroy a LatLon" << std::endl;
}

std::string LatLon::hash() const
{
    NOTIMP;
}

Grid::BoundBox LatLon::boundingBox() const
{
    return bound_box_;
}

void LatLon::coordinates( Grid::Coords& r ) const
{
    ASSERT( r.size() == points_.size() );

    for( size_t i = 0; i < points_.size(); ++i )
    {
        r.lat(i) = points_[i].lat();
        r.lon(i) = points_[i].lon();
    }
}

GridSpec* LatLon::spec() const
{
   GridSpec* grid_spec = new GridSpec(gridType(),"LL");

   grid_spec->set("nlat",eckit::Value(nlat_));
   grid_spec->set("nlon",eckit::Value(nlon_));

   grid_spec->set_bounding_box(bound_box_);
   grid_spec->set_points(points_);

   return grid_spec;
}

void LatLon::constructFrom(const GridSpec& grid_spec)
{
   nlat_ = grid_spec.get("nlat");
   nlon_ = grid_spec.get("nlon");
   grid_spec.get_bounding_box(bound_box_);
   grid_spec.get_points(points_);
}

//-----------------------------------------------------------------------------

} // namespace grid
} // namespace eckit
