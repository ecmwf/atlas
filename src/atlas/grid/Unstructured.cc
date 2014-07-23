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
#include "eckit/memory/Builder.h"

#include "atlas/grid/GridSpec.h"
#include "atlas/grid/Unstructured.h"

using namespace eckit;

namespace atlas {
namespace grid {

//-----------------------------------------------------------------------------

ConcreteBuilderT1<Grid,Unstructured> unstructured_builder;

//void Unstructured::constructFrom(const GridSpec& grid_spec)
//{
//    if (grid_spec.has("hash")) hash_ = (std::string)grid_spec.get("hash");
//    grid_spec.get_bounding_box(bound_box_);
//       std::vector< Grid::Point >* pts = new std::vector< Grid::Point >(0);
//       grid_spec.get_points(*pts);
//       points_.reset(pts);
//}

Unstructured::Unstructured(const Params& p)
{
	NOTIMP;
}

Unstructured::Unstructured( std::vector< Point >* pts, const std::string& hash ) :
    points_(pts),
    hash_(hash)
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

void Unstructured::coordinates(std::vector<double>& r ) const
{
	NOTIMP;
}

void Unstructured::coordinates(std::vector<Grid::Point>&) const
{
	NOTIMP;
}

GridSpec* Unstructured::spec() const
{
	NOTIMP;

	GridSpec* grid_spec = new GridSpec(gridType(),"U");

	grid_spec->set("hash",eckit::Value(hash_));
	grid_spec->set_bounding_box(bound_box_);
	//   grid_spec->set_points(coordinates());

	return grid_spec;
}

bool Unstructured::same(const Grid& grid) const
{
   if (gridType() != grid.gridType()) return false;

   if ( static_cast<const Unstructured&>(grid).hash_ != hash_) return false;
   if ( static_cast<const Unstructured&>(grid).bound_box_ != bound_box_) return false;
   if ( *static_cast<const Unstructured&>(grid).points_ != *points_) return false;

   return true;
}

//-----------------------------------------------------------------------------

} // namespace grid
} // namespace eckit
