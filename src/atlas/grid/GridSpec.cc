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


Grid::Ptr GridSpec::build( const GridSpec& spec)
{

   return Grid::Ptr();
}

//------------------------------------------------------------------------------------------------------

} // namespace grid
} // namespace atlas

