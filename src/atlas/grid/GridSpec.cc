/*
 * (C) Copyright 1996-2014 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "eckit/exception/Exceptions.h"
#include "atlas/grid/GridSpec.h"

//------------------------------------------------------------------------------------------------------

namespace atlas {
namespace grid {

//------------------------------------------------------------------------------------------------------

GridSpec::GridSpec(const std::string& grid_type)
{
   set("grid_type",Properties::property_t(grid_type));
}

GridSpec::~GridSpec(){}

std::string GridSpec::grid_type() const
{
   Properties::property_t val = get("grid_type");
   if (val.isNil()) {
      throw eckit::SeriousBug("GridSpec with no grid type specified", Here());
   }
   return val.as<std::string>();
}

void GridSpec::uid(const std::string& uid)
{
   set("uid",Properties::property_t(uid));
}

std::string GridSpec::uid() const
{
   Properties::property_t val = get("shortName");
   if (val.isNil()) {
      throw eckit::SeriousBug("GridSpec with no short name specified", Here());
   }
   return val.as<std::string>();
}

void GridSpec::print_simple(std::ostream& s) const
{
   s << "GridSpec[ ";
   s << "gridType:" << grid_type();
   s << ", shortName:" << uid();

   Grid::BoundBox bbox;
   get_bounding_box( bbox ) ;
   s << ", bbox:(" << bbox.bottom_left().lat() << "," << bbox.bottom_left().lon() << "," << bbox.top_right().lat() << "," << bbox.top_right().lon() << ")";

   s << " ]";
}

void GridSpec::print( std::ostream& s) const
{
   s << "GridSpec[ ";
   Properties::print(s) ;
   s << " ]";
}

void GridSpec::set_points(const std::vector<Grid::Point>& points)
{
   points_ = points;

//   std::vector<eckit::Value> points_lat; points_lat.reserve(points.size());
//   std::vector<eckit::Value> points_lon; points_lon.reserve(points.size());
//   for(size_t i = 0; i < points.size(); ++i) {
//      points_lat.push_back(eckit::Value(points[i].lat()));
//      points_lon.push_back(eckit::Value(points[i].lon()));
//   }
//   set("points_lat",eckit::Value(points_lat) );
//   set("points_lon",eckit::Value(points_lon) );
}

void GridSpec::set_latitudes(const std::vector<double>& latitudes_vec)
{
   std::vector<eckit::Value> latitudes; latitudes.reserve(latitudes_vec.size());
   for(size_t i = 0; i < latitudes_vec.size(); ++i) {
      latitudes.push_back(eckit::Value(latitudes_vec[i]));
   }
   set("latitudes",eckit::Value(latitudes) );
}

void GridSpec::set_rgspec(const std::vector<long>& rgSpec_vec)
{
   std::vector<eckit::Value> rgSpec; rgSpec.reserve(rgSpec_vec.size());
   for(size_t i = 0; i < rgSpec_vec.size(); ++i) {
      rgSpec.push_back(eckit::Value(rgSpec_vec[i]));
   }
   set("rgSpec",eckit::Value(rgSpec) );
}

void GridSpec::set_bounding_box(const Grid::BoundBox& bbox )
{
   set("area_s", bbox.bottom_left().lat());
   set("area_w", bbox.bottom_left().lon());
   set("area_n", bbox.top_right().lat());
   set("area_e", bbox.top_right().lon());
}

void GridSpec::get_points(std::vector<Grid::Point>& points) const
{
   points = points_;

//   if (has("points_lat")) {
//      eckit::ValueList vec_lat = get("points_lat");
//      eckit::ValueList vec_lon = get("points_lon");
//      ASSERT(vec_lat.size() ==  vec_lon.size());
//      points.reserve(vec_lat.size());
//      for(size_t i=0; i < vec_lat.size();++i) {
//         points.push_back(Grid::Point(vec_lat[i],vec_lon[i]));
//      }
//   }
}

void GridSpec::get_latitudes(std::vector<double>& latitudes) const
{
   if (has("latitudes")) {
      eckit::ValueList vec = get("latitudes");
      latitudes.reserve(vec.size());
      for(size_t i=0; i < vec.size();++i) {
         latitudes.push_back(vec[i]);
      }
   }
}

void GridSpec::get_rgspec(std::vector<long>& rgSpec) const
{
   if (has("rgSpec")) {
      eckit::ValueList vec = get("rgSpec");
      rgSpec.reserve(vec.size());
      for(size_t i=0; i < vec.size();++i) {
         rgSpec.push_back(vec[i]);
      }
   }
}

void GridSpec::get_bounding_box(Grid::BoundBox& bbox) const
{
   if (has("area_s")) {
      double area_s = get("area_s");
      double area_w = get("area_w");
      double area_n = get("area_n");
      double area_e = get("area_e");
      bbox.bottom_left().assign(area_s,area_w);
      bbox.top_right().assign(area_n,area_e);
      ASSERT( bbox.validate() );
   }
}

//------------------------------------------------------------------------------------------------------

} // namespace grid
} // namespace atlas

