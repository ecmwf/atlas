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
#include "eckit/parser/JSON.h"
#include "eckit/utils/MD5.h"

#include "atlas/GridSpec.h"

//------------------------------------------------------------------------------------------------------

namespace atlas {


//------------------------------------------------------------------------------------------------------

GridSpec::GridSpec(const std::string& grid_type)
{
   set("grid_type",Properties::property_t(grid_type));
   ASSERT(!empty());
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

Grid::uid_t GridSpec::uid() const
{
   Properties::property_t val = get("uid");
   if(val.isNil())
      throw eckit::SeriousBug("GridSpec with no short name specified", Here());

   return val.as<std::string>();
}

std::string GridSpec::hash() const
{
   return eckit::MD5( str() );
}

void GridSpec::print( std::ostream& s) const
{
   eckit::JSON js(s);
   js.precision(16);
   js << *this;
}

void GridSpec::set_latitudes(const std::vector<double>& latitudes)
{
   std::vector<eckit::Value> lats; lats.reserve(latitudes.size());
   for(size_t i = 0; i < latitudes.size(); ++i)
   {
      lats.push_back(eckit::Value(latitudes[i]));
   }
   set("latitudes",eckit::Value(lats) );
}

void GridSpec::set_npts_per_lat(const std::vector<int>& rgSpec_vec)
{
   std::vector<eckit::Value> rgSpec; rgSpec.reserve(rgSpec_vec.size());
   for(size_t i = 0; i < rgSpec_vec.size(); ++i)
   {
      rgSpec.push_back(eckit::Value(rgSpec_vec[i]));
   }
   set("npts_per_lat",eckit::Value(rgSpec) );
}

void GridSpec::set_bounding_box(const Grid::BoundBox& bbox )
{
   set("bbox_s", bbox.min().lat());
   set("bbox_w", bbox.min().lon());
   set("bbox_n", bbox.max().lat());
   set("bbox_e", bbox.max().lon());
}

void GridSpec::get_npts_per_lat(std::vector<int> &rgSpec) const
{
   if(has("npts_per_lat"))
   {
      eckit::ValueList vec = get("npts_per_lat");
      rgSpec.reserve(vec.size());
      for(size_t i=0; i < vec.size();++i)
      {
         rgSpec.push_back(vec[i]);
      }
   }
}

void GridSpec::get_bounding_box(Grid::BoundBox& bbox) const
{
   if(has("bbox_s"))
   {
      double area_s = get("bbox_s");
      double area_w = get("bbox_w");
      double area_n = get("bbox_n");
      double area_e = get("bbox_e");
      bbox.min().assign(area_w,area_s);
      bbox.max().assign(area_e,area_n);
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


} // namespace atlas

