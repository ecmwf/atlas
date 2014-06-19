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
#include "eckit/value/Properties.h"

//------------------------------------------------------------------------------------------------------

namespace atlas {
namespace grid {

//------------------------------------------------------------------------------------------------------

/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////8
/// A concrete class that holds specification that uniquely identifies a Grid
/// The description of the grid is added as name value pairs
/// This class will provides a short name for a GRID (i.e QG48_1)
/// This allows for easier matching with samples files.
/// However this interface is independent of GRIB/NETCDF
///
/// Uses default copy constructor, assignment and equality operators
/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////8

class GridSpec : public eckit::Properties {
public:

   GridSpec(const std::string& the_grid_type);
   GridSpec(const std::string& the_grid_type, const std::string& the_short_name);
   ~GridSpec();

   /// returns the gridType. currently this matches grid _type found in GRIB
   std::string grid_type() const;

   /// returns short name description of this grid
   /// We use the following prefix:
   /// 'LL'  -  regular lat/lon grid                             - other centres, wave data.
   /// 'NP'  -  northern polar ster. projection (reg. lat/lon)
   /// 'SP'  -  southern polar ster. projection (reg. lat/lon)
   /// 'GG'  -  regular gaussian grid (surface)                  - Surface and some upper air fields.
   /// 'SH'  -  spherical harmonics coefficients                 - Upper air fields.
   /// 'QG'  -  quasi regular gaussian grid
   /// 'RL'  -  rotated lat long
   /// 'RedLL'- reduced lat long
   void set_short_name(const std::string& the_short_name);
   std::string short_name() const;

   friend std::ostream& operator<<( std::ostream& os, const GridSpec& v) { v.print(os); return os;}

private:

   void print( std::ostream& ) const;
};

//------------------------------------------------------------------------------------------------------

} // namespace grid
} // namespace atlas

#endif
