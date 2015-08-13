/*
 * (C) Copyright 1996-2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include <cmath>

#include "eckit/log/Log.h"
#include "eckit/memory/Builder.h"
#include "eckit/value/Value.h"
#include "eckit/geometry/PolarStereoGraphicProj.h"
#include "eckit/types/FloatCompare.h"

#include "atlas/grids/PolarStereoGraphic.h"

using namespace eckit;
using namespace eckit::geometry;
using namespace std;

namespace atlas {
namespace grids {

//------------------------------------------------------------------------------------------------------

register_BuilderT1(Grid, PolarStereoGraphic, PolarStereoGraphic::grid_type_str());

PolarStereoGraphic::PolarStereoGraphic(size_t Nx, size_t Ny,
                     size_t Dx, size_t Dy,
                     double longitudeOfFirstGridPoint,
                     double latitudeOfFirstGridPoint,
                     double orientationOfTheGrid,
                     bool southPoleOnProjectionPlane,
                     double radius,
                     bool earth_is_oblate,
                     double semi_major,
                     double semi_minor):
    npts_xaxis_(Nx),
    npts_yaxis_(Ny),
    x_grid_length_(Dx),
    y_grid_length_(Dy),
    orientationOfTheGrid_(orientationOfTheGrid),
    radius_(radius),
    semi_major_(semi_major),
    semi_minor_(semi_minor),
    first_grid_pt_(longitudeOfFirstGridPoint, latitudeOfFirstGridPoint),
    southPoleOnProjectionPlane_(southPoleOnProjectionPlane),
    earth_is_oblate_(earth_is_oblate)

{

    ASSERT( npts_xaxis_ > 0);
    ASSERT( npts_yaxis_ > 0);
    ASSERT( x_grid_length_ > 0);
    ASSERT( y_grid_length_ > 0);
    if (earth_is_oblate_) {
        ASSERT( semi_major_ > semi_minor_);
    }

}

PolarStereoGraphic::PolarStereoGraphic( const eckit::Parametrisation &p )
    : npts_xaxis_(0),
      npts_yaxis_(0),
      x_grid_length_(0),
      y_grid_length_(0),
      orientationOfTheGrid_(0),
      radius_(6371229),
      semi_major_(6378137),
      semi_minor_(6356752.3),
      southPoleOnProjectionPlane_(false),
      earth_is_oblate_(false)

{

    NOTIMP;
}

PolarStereoGraphic::~PolarStereoGraphic() {
}


eckit::Properties PolarStereoGraphic::spec() const {
    NOTIMP;
}

void PolarStereoGraphic::print(ostream& os) const
{
    os << "PolarStereoGraphic("
          << "npts_xaxis:" << npts_xaxis_
          << "npts_yaxis:" << npts_yaxis_
          << "x_grid_length:" << x_grid_length_
          << "y_grid_length:" << y_grid_length_
          << "first_grid_pt:" << first_grid_pt_
          << "orientationOfTheGrid:" << orientationOfTheGrid_
          << "southPoleOnProjectionPlane:" << southPoleOnProjectionPlane_
          << "earth_is_oblate:" << earth_is_oblate_
          << "semi_major:" << semi_major_
          << "semi_minor:" << semi_minor_
          << "radius:" << radius_
          << ")";
}

std::string PolarStereoGraphic::shortName() const {

    if (shortName_.empty()) {
        std::ostringstream s;
        s << grid_type_str() << "." << npts_xaxis_ << "x" << npts_yaxis_;
        shortName_ = s.str();
    }

    return shortName_;
}

void PolarStereoGraphic::hash(eckit::MD5 &md5) const {

    md5.add(grid_type_str());

    md5.add(npts_xaxis_);
    md5.add(npts_yaxis_);

    md5.add(x_grid_length_);
    md5.add(y_grid_length_);

    md5.add(first_grid_pt_.lat());
    md5.add(first_grid_pt_.lon());

    md5.add(orientationOfTheGrid_);
    md5.add(southPoleOnProjectionPlane_);

    md5.add(earth_is_oblate_);
    md5.add(semi_major_);
    md5.add(semi_minor_);
    md5.add(radius_);

}

size_t PolarStereoGraphic::npts() const {
    return npts_xaxis_ * npts_yaxis_;
}

//Enable to test against grib iterator
//#define GRIB_COMPAT 1

void PolarStereoGraphic::lonlat(std::vector<Grid::Point> &points) const {
    points.resize(npts());

    PolarStereoGraphicProj ps(southPoleOnProjectionPlane_, earth_is_oblate_, orientationOfTheGrid_);
    if (earth_is_oblate_) {
        ps.set_radius(semi_major_);
        ps.set_eccentricity(sqrt( 1.0 - (semi_minor_ * semi_minor_) / (semi_major_ * semi_major_)));
    } else {
        ps.set_radius(radius_);
    }

    // Points go from North,West --> South,east
    BoundBox bbox = boundingBox();
    const double north = bbox.north();
    const double west  = bbox.west();

    Point2 first_pt_on_plane = ps.map_to_plane(eckit::geometry::LLPoint2(north, west));
    double x = first_pt_on_plane[0];
    double y = first_pt_on_plane[1];
    size_t k = 0;
    for (long j = 0; j < npts_yaxis_; j++) {
        x = first_pt_on_plane[0];
        for (long i = 0; i < npts_xaxis_; i++) {

            points[k++] = ps.map_to_spherical(x, y);
            x += x_grid_length_;
        }
        y -= y_grid_length_;
    }

}

BoundBox PolarStereoGraphic::boundingBox() const {
    // Map first_grid_pt_ to the plane.
    PolarStereoGraphicProj ps(southPoleOnProjectionPlane_, earth_is_oblate_, orientationOfTheGrid_);
    if (earth_is_oblate_) {
        ps.set_radius(semi_major_);
        ps.set_eccentricity(sqrt( 1.0 - (semi_minor_ * semi_minor_) / (semi_major_ * semi_major_)));
    } else {
        ps.set_radius(radius_);
    }

    Point2 first_pt_on_plane = ps.map_to_plane(first_grid_pt_);

    // *Depending* on the scanning mode, find the last point on the plane
    // Note: without the scanning mode we can not find the correct bounding box.
    double last_point_x = first_pt_on_plane[0] + npts_xaxis_ * x_grid_length_;
    double last_point_y = first_pt_on_plane[1] + npts_yaxis_ * y_grid_length_;

    // Transform this back into spherical co-ordinates
    Point last_pt_in_sp = ps.map_to_spherical(last_point_x, last_point_y);

    double top = std::max(first_grid_pt_.lat(), last_pt_in_sp.lat());
    double bottom = std::min(first_grid_pt_.lat(), last_pt_in_sp.lat());
    double right = std::max(first_grid_pt_.lon(), last_pt_in_sp.lon());
    double left = std::min(first_grid_pt_.lon(), last_pt_in_sp.lon());

    return BoundBox(top, bottom, right, left);
}

string PolarStereoGraphic::gridType() const {
    return PolarStereoGraphic::grid_type_str();
}

//-----------------------------------------------------------------------------


} // namespace grids
} // namespace atlas
