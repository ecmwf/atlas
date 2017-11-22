/*
 * (C) Copyright 1996-2017 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "atlas/functionspace/PointCloud.h"
#include "atlas/array.h"

namespace atlas {
namespace functionspace {

namespace detail {

PointCloud::PointCloud(const std::vector<PointXY>& points) {
  lonlat_ = Field("lonlat", array::make_datatype<double>(), array::make_shape(points.size(),2));
  auto lonlat = array::make_view<double,2>(lonlat_);
  for( size_t j=0; j<points.size(); ++j ) {
    lonlat(j,0) = points[j].x();
    lonlat(j,1) = points[j].y();
  }
}

PointCloud::PointCloud(const Field& lonlat) :
  lonlat_(lonlat) {
}

PointCloud::PointCloud(const Field& lonlat, const Field& ghost) :
  lonlat_(lonlat),
  ghost_(ghost) {
}

const Field& PointCloud::ghost() const {
 if (not ghost_) {
   ghost_ = Field( "ghost", array::make_datatype<int>(), array::make_shape(size()) );
   array::make_view<int,1>(ghost_).assign(0);
 }
 return ghost_;
}


Field PointCloud::createField(const eckit::Configuration& options) const {
  NOTIMP;
  return Field();
}

Field PointCloud::createField(
    const Field& other,
    const eckit::Configuration& config ) const
{
  return createField(
    option::datatype ( other.datatype()  ) |
    config );
}

std::string PointCloud::distribution() const {
  return std::string("serial");
}

}

PointCloud::PointCloud( const FunctionSpace& functionspace ) :
  FunctionSpace(functionspace),
  functionspace_( dynamic_cast< const detail::PointCloud* >( get() ) ) {
}

PointCloud::PointCloud( const Field& points ) :
  FunctionSpace( new detail::PointCloud(points) ),
  functionspace_( dynamic_cast< const detail::PointCloud* >( get() ) ) {
}

PointCloud::PointCloud( const std::vector<PointXY>& points ) :
  FunctionSpace( new detail::PointCloud(points) ),
  functionspace_( dynamic_cast< const detail::PointCloud* >( get() ) ) {
}


}
}
