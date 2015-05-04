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
#include "atlas/GridSpec.h"
#include "atlas/grids/Unstructured.h"
#include "atlas/util/ArrayView.h"

using eckit::MD5;

namespace atlas {
namespace grids {

register_BuilderT1(Grid,Unstructured,Unstructured::grid_type_str());

//void Unstructured::constructFrom(const GridSpec& grid_spec)
//{
//    if (grid_spec.has("hash")) hash_ = (std::string)grid_spec.get("hash");
//    grid_spec.get_bounding_box(bound_box_);
//       std::vector< Grid::Point >* pts = new std::vector< Grid::Point >(0);
//       grid_spec.get_points(*pts);
//       points_.reset(pts);
//}


Unstructured::Unstructured(const eckit::Params& p)
{
  NOTIMP;
}


Unstructured::Unstructured( std::vector< Point >* pts ) :
  points_(pts)
{
  const std::vector<Point>& p = *points_;
  const size_t npts = p.size();

  double lat_min = std::numeric_limits<double>::max();
  double lat_max = std::numeric_limits<double>::min();
  double lon_min = lat_min;
  double lon_max = lat_max;

  for (size_t n=0; n<npts; ++n) {
    lat_min = std::min( lat_min, p[n].lat() );
    lat_max = std::max( lat_max, p[n].lat() );
    lon_min = std::min( lon_min, p[n].lon() );
    lon_max = std::max( lon_max, p[n].lon() );
  }

  bound_box_ = BoundBox( Point(lon_min,lat_min), Point(lon_max,lat_max) );
}


Unstructured::~Unstructured()
{
}


Grid::uid_t Unstructured::shortName() const
{
  if( shortName_.empty() )
  {
    std::ostringstream s;
    s <<  "Unst." << hash().substr(0, 7) << eckit::StrStream::ends;
    shortName_ = s.str();
  }
  return shortName_;
}

Grid::uid_t Unstructured::unique_id() const {

  if (uid_.empty()) {
    std::ostringstream s;
    s <<  "Unst." << hash() << eckit::StrStream::ends;
    uid_ = s.str();
  }

  return uid_;
}

MD5::digest_t Unstructured::hash() const {

  if (hash_.empty()) {

    ASSERT(points_);
    const std::vector< Point >& pts = *points_;
    eckit::MD5 md5;
    md5.add(&pts[0], sizeof(Point)*pts.size());
    hash_ = md5.digest();
  }

  return hash_;
}

BoundBox Unstructured::bounding_box() const
{
  return bound_box_;
}


size_t Unstructured::npts() const
{
  return points_->size();
}


void Unstructured::lonlat(double crds[]) const
{
  for (size_t i=0, c=0; i<npts(); ++i) {
    crds[c++] = (*points_)[i].lon();
    crds[c++] = (*points_)[i].lat();
  }
}

void Unstructured::lonlat(std::vector< Grid::Point >& crds) const
{
  crds.resize(npts());
  for (size_t i=0; i<npts(); ++i)
    crds[i].assign(
          (*points_)[i].lon(),
          (*points_)[i].lat() );
}

void Unstructured::lonlat(std::vector< double >& v) const
{
  Grid::lonlat(v);
}


GridSpec Unstructured::spec() const
{
  if(cachedGridSpec_)
    return *cachedGridSpec_;

  cachedGridSpec_.reset( new GridSpec(grid_type()) );

  cachedGridSpec_->set_bounding_box(bound_box_);

  std::vector<double> coords;
  lonlat(coords);

  cachedGridSpec_->set( "lonlat", eckit::makeVectorValue<double>(coords) );

  return *cachedGridSpec_;
}

} // namespace grids
} // namespace atlas
