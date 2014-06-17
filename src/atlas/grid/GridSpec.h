#ifndef atlas_grid_GridSpec_H
#define atlas_grid_GridSpec_H
/*
 * (C) Copyright 1996-2014 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include <cstddef>
#include <map>

#include "eckit/value/Value.h"

//------------------------------------------------------------------------------------------------------

namespace atlas {
namespace grid {

//------------------------------------------------------------------------------------------------------

/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////8
/// A concrete class that holds specification that uniquely identifies a Grid
/// The description of the grid is added as name value pairs
/// This class will provides a short name for a GRID (i.e QG48_1)
/// This allows for easier matching with samples files.
///
/// Uses default copy constructor, assignment and equality operators
/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////8

class GridSpec  {
public:

   GridSpec(const std::string& the_grid_type);
   GridSpec(const std::string& the_grid_type, const std::string& the_short_name);
   ~GridSpec();

   /// returns the gridType. currently this matches grid _type found in GRIB
   std::string grid_type() const { return the_grid_type_; }

   /// returns short name description of this grid
   // 'LL'  -  regular lat/lon grid                             - other centres, wave data.
   // 'NP'  -  northern polar ster. projection (reg. lat/lon)
   // 'SP'  -  southern polar ster. projection (reg. lat/lon)
   // 'GG'  -  regular gaussian grid (surface)                  - Surface and some upper air fields.
   // 'SH'  -  spherical harmonics coefficients                 - Upper air fields.
   // 'QG'  -  quasi regular gaussian grid
   // 'RL'  -  rotated lat long
   // 'RedLL' = reduced lat long
   // Additionally we use _1/_2 to specify grib1/grib2 ??
   void set_short_name(const std::string& the_short_name) { the_short_name_ = the_short_name; }
   std::string short_name() const { return the_short_name_; }

   /// Used to build up description of a grid
   void add(const std::string&, const eckit::Value& );

   /// Find the key and return the value, if key NOT found returns a Nil value
   eckit::Value find(const std::string& key) const;

   friend std::ostream& operator<<( std::ostream& os, const GridSpec& v) { v.print(os); return os;}

private:

   void print( std::ostream& ) const;

private:

   std::string the_grid_type_;
   std::string the_short_name_;
   std::map<std::string,eckit::Value> grid_spec_;
};

//------------------------------------------------------------------------------------------------------

} // namespace grid
} // namespace atlas

#endif
