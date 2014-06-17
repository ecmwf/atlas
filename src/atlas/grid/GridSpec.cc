/*
 * (C) Copyright 1996-2014 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "atlas/grid/GridSpec.h"

//------------------------------------------------------------------------------------------------------

namespace atlas {
namespace grid {

//------------------------------------------------------------------------------------------------------


GridSpec::GridSpec(const std::string& the_grid_type, const std::string& the_short_name)
: the_grid_type_(the_grid_type),the_short_name_(the_short_name)
{
}

GridSpec::GridSpec(const std::string& the_grid_type)
: the_grid_type_(the_grid_type)
{
}

GridSpec::~GridSpec()
{
}

void GridSpec::add(const std::string& name, const eckit::Value& value)
{
   grid_spec_.insert(std::make_pair(name,value));
}

eckit::Value GridSpec::find(const std::string& key) const
{
   std::map<std::string,eckit::Value>::const_iterator i = grid_spec_.find(key);
   if (i != grid_spec_.end()) {
      return (*i).second;
   }
   return eckit::Value();
}

void GridSpec::print( std::ostream& s) const
{
   s << "GridSpec[ " << the_grid_type_;
   if (!the_short_name_.empty() ) s << ", " << the_short_name_;

   std::map<std::string,eckit::Value>::const_iterator i = grid_spec_.begin();
   for(i = grid_spec_.begin(); i != grid_spec_.end(); ++i) {
      s << ", " << (*i).first << ":" << (*i).second;
   }

   s << " ]";
}

//------------------------------------------------------------------------------------------------------

} // namespace grid
} // namespace atlas

