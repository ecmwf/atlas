/*
 * (C) Copyright 1996-2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

// ==================================================================================
// gribs use the following convention: (from Shahram)
//
// Horizontally:  Points scan in the +i (+x) direction
// Vertically:    Points scan in the -j (-y) direction
//
// The way I verified this was to look at our SAMPLE files (which IFS uses).
// I also verified that IFS does not modify the scanning modes
// so whatever the samples say, is the convention
// ==================================================================================

#include <stdexcept>

#include "eckit/log/Log.h"
#include "eckit/grib/GribAccessor.h"
#include "atlas/grid/GribRead.h"
#include "atlas/grid/GribGrid.h"

using namespace eckit;
using namespace std;

namespace atlas {
namespace grid {

//-----------------------------------------------------------------------------
// Area: Do we check the area.
// Area: Can we assume area is multiple of the grids ?

GribGrid::GribGrid(grib_handle* handle)
: editionNumber_(0),
  north_(0.0),south_(0.0),west_(0.0),east_(0.0),
  epsilon_(1e-6),
  numberOfDataPoints_(0)
{
   if (handle == NULL)
      throw std::runtime_error("NULL grib_handle");

   GRIB_CHECK(grib_get_long(handle,"editionNumber",&editionNumber_),0);
   epsilon_ = (editionNumber_ == 1) ? 1e-3 : 1e-6;


   hash_ = grib_hash(handle);

   GRIB_CHECK(grib_get_long(handle,"iScansNegatively",&iScansNegatively_),0);
   GRIB_CHECK(grib_get_long(handle,"jScansPositively",&jScansPositively_),0);

   GRIB_CHECK(grib_get_double(handle,"latitudeOfFirstGridPointInDegrees",&north_),0);
   GRIB_CHECK(grib_get_double(handle,"longitudeOfFirstGridPointInDegrees",&west_),0);
   GRIB_CHECK(grib_get_double(handle,"latitudeOfLastGridPointInDegrees",&south_),0);
   GRIB_CHECK(grib_get_double(handle,"longitudeOfLastGridPointInDegrees",&east_),0);


   grib_get_long(handle,"numberOfDataPoints",&numberOfDataPoints_);
}

GribGrid::~GribGrid(){}

//-----------------------------------------------------------------------------

} // namespace grid
} // namespace eckit
