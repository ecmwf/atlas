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

GridSpec::GridSpec(const std::string& the_grid_type, const std::string& the_short_name)
{
   set("gridType",Properties::property_t(the_grid_type));
   set("shortName",Properties::property_t(the_short_name));
}

GridSpec::GridSpec(const std::string& the_grid_type)
{
   set("gridType",Properties::property_t(the_grid_type));
}

GridSpec::~GridSpec(){}


std::string GridSpec::grid_type() const
{
   Properties::property_t val = get("gridType");
   if (val.isNil()) {
      throw eckit::SeriousBug("GridSpec with no grid type specified", Here());
   }

   std::string the_grid_type = val;
   return the_grid_type;
}

void GridSpec::set_short_name(const std::string& the_short_name)
{
   set("shortName",Properties::property_t(the_short_name));
}

std::string GridSpec::short_name() const
{
   Properties::property_t val = get("shortName");
   if (val.isNil()) {
      throw eckit::SeriousBug("GridSpec with no short name specified", Here());
   }
   std::string the_short_name = val;
   return the_short_name;
}


void GridSpec::print( std::ostream& s) const
{
   s << "GridSpec[ ";
   Properties::print(s) ;
   s << " ]";
}

void GridSpec::set_points(const std::vector<Grid::Point>& points)
{
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
   set("bottom_left_lat",eckit::Value(bbox.bottom_left_.lat()));
   set("bottom_left_lon",eckit::Value(bbox.bottom_left_.lon()));
   set("top_right_lat",eckit::Value(bbox.top_right_.lat()));
   set("top_right_lon",eckit::Value(bbox.top_right_.lon()));
}

void GridSpec::get_points(std::vector<Grid::Point>& points) const
{
   eckit::ValueList vec_lat = get("points_lat");
   eckit::ValueList vec_lon = get("points_lon");
   ASSERT(vec_lat.size() ==  vec_lon.size());
   points.reserve(vec_lat.size());
   for(size_t i=0; i < vec_lat.size();++i) {
      points.push_back(Grid::Point(vec_lat[i],vec_lon[i]));
   }
}

void GridSpec::get_latitudes(std::vector<double>& latitudes) const
{
   eckit::ValueList vec = get("latitudes");
   latitudes.reserve(vec.size());
   for(size_t i=0; i < vec.size();++i) {
      latitudes.push_back(vec[i]);
   }
}

void GridSpec::get_rgspec(std::vector<long>& rgSpec) const
{
   eckit::ValueList vec = get("rgSpec");
   rgSpec.reserve(vec.size());
   for(size_t i=0; i < vec.size();++i) {
      rgSpec.push_back(vec[i]);
   }
}

void GridSpec::get_bounding_box(Grid::BoundBox& bbox ) const
{
   double bottom_left_lat = get("bottom_left_lat");
   double bottom_left_lon = get("bottom_left_lon");
   double top_right_lat = get("top_right_lat");
   double top_right_lon = get("top_right_lon");
   bbox.bottom_left_.assign(bottom_left_lat,bottom_left_lon);
   bbox.top_right_.assign(top_right_lat,top_right_lon);
}

//------------------------------------------------------------------------------------------------------

} // namespace grid
} // namespace atlas

