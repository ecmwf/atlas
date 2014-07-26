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

GridSpec::~GridSpec()
{
}

std::string GridSpec::grid_type() const
{
   Properties::property_t val = get("grid_type");
   if( val.isNil() )
      throw eckit::SeriousBug("GridSpec with no grid type specified", Here());

   return val.as<std::string>();
}

void GridSpec::uid(const std::string& uid)
{
   set("uid",Properties::property_t(uid));
}

std::string GridSpec::uid() const
{
   Properties::property_t val = get("uid");
   if(val.isNil())
      throw eckit::SeriousBug("GridSpec with no short name specified", Here());

   return val.as<std::string>();
}

void GridSpec::print( std::ostream& s) const
{
   s << "GridSpec[ ";
   Properties::print(s) ;
   s << " ]";
}

void GridSpec::set_rgspec(const std::vector<long>& rgSpec_vec)
{
   std::vector<eckit::Value> rgSpec; rgSpec.reserve(rgSpec_vec.size());
   for(size_t i = 0; i < rgSpec_vec.size(); ++i)
   {
      rgSpec.push_back(eckit::Value(rgSpec_vec[i]));
   }
   set("rgSpec",eckit::Value(rgSpec) );
}

void GridSpec::set_bounding_box(const Grid::BoundBox& bbox )
{
   set("grid_bbox_s", bbox.bottom_left().lat());
   set("grid_bbox_w", bbox.bottom_left().lon());
   set("grib_bbox_n", bbox.top_right().lat());
   set("grid_bbox_e", bbox.top_right().lon());
}

void GridSpec::get_rgspec(std::vector<long>& rgSpec) const
{
   if(has("rgSpec"))
   {
      eckit::ValueList vec = get("rgSpec");
      rgSpec.reserve(vec.size());
      for(size_t i=0; i < vec.size();++i)
      {
         rgSpec.push_back(vec[i]);
      }
   }
}

void GridSpec::get_bounding_box(Grid::BoundBox& bbox) const
{
   if(has("grid_bbox_s"))
   {
      double area_s = get("grid_bbox_s");
      double area_w = get("grid_bbox_w");
      double area_n = get("grib_bbox_n");
      double area_e = get("grid_bbox_e");
      bbox.bottom_left().assign(area_s,area_w);
      bbox.top_right().assign(area_n,area_e);
      ASSERT( bbox.validate() );
   }
}

std::string GridSpec::str() const
{
	std::ostringstream oss;
	print(oss);
	return oss.str();
}

//------------------------------------------------------------------------------------------------------

} // namespace grid
} // namespace atlas

