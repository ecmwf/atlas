#ifndef atlas_grib_grid_H
#define atlas_grib_grid_H
/*
 * (C) Copyright 1996-2014 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "grib_api.h"

#include "atlas/grid/Grid.h"

//-----------------------------------------------------------------------------

namespace atlas {
namespace grid {

//-----------------------------------------------------------------------------

class GribGrid : public Grid {
public:
   virtual ~GribGrid();

   /// Overridden functions
   virtual std::string hash() const { return hash_;}

protected:
   GribGrid( grib_handle* h );

   // grib edition 1 - milli-degrees
   // grib edition 2 - micro-degrees or could be defined by the keys: "subdivisionsOfBasicAngle" and "basicAngleOfTheInitialProductionDomain"
   // Therefore the client needs access to this when dealing with double based comparisons (for tolerances)
   double epsilon() const { return epsilon_; }

protected:
   long   editionNumber_;           /// Grib 1 or Grib 2
   double north_;                   /// In degrees
   double south_;                   /// In degrees
   double west_;                    /// In degrees
   double east_;                    /// In degrees
   double epsilon_;                 /// Grib 1 or Grib 2
   long   numberOfDataPoints_;      ///
   long   iScansNegatively_;
   long   jScansPositively_;
   std::string hash_;
};

//-----------------------------------------------------------------------------

} // namespace grid
} // namespace eckit

#endif
