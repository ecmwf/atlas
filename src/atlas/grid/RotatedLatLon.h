/*
 * (C) Copyright 1996-2014 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef atlas_grid_RotatedLatLon_H
#define atlas_grid_RotatedLatLon_H

#include <cstddef>
#include <vector>

#include "atlas/grid/Grid.h"


//-----------------------------------------------------------------------------

namespace atlas {
namespace grid {

//-----------------------------------------------------------------------------

/// @note
/// gribs use the following convention: (from Shahram)
///
/// Horizontally:  Points scan in the +i (+x) direction
/// Vertically:    Points scan in the -j (-y) direction
///
/// The way I verified this was to look at our SAMPLE files (which IFS uses).
/// I also verified that IFS does not modify the scanning modes
/// so whatever the samples say, is the convention
///
/// @todo Do we check the area? Can we assume area is multiple of the grids ?

class RotatedLatLon : public Grid {

public: // methods

   RotatedLatLon( const eckit::Params& p );
   virtual ~RotatedLatLon();

   virtual std::string hash() const { return hash_;}
   virtual BoundBox boundingBox() const { return bbox_;}
   virtual size_t nPoints() const { return points_.size(); }

   virtual void coordinates( std::vector<double>& ) const;
   virtual void coordinates( std::vector<Point>& ) const;

   virtual std::string gridType() const { return std::string("rotated_ll"); }
   virtual GridSpec* spec() const;
   virtual bool same(const Grid&) const;

protected: // methods

   double rotated_latitude() const { return rotated_latitude_; }
   double rotated_longitude() const { return rotated_longitude_; }
   double rotated_angle() const { return rotated_angle_; }

   Point latLon(size_t lat, size_t lon) const;
   long rows() const { return nptsNS_;}
   long cols() const { return nptsWE_;}
   double incLat() const { return nsIncrement_; }
   double incLon() const { return weIncrement_; }

private: // members

   std::string hash_;
   BoundBox bbox_;
   double rotated_latitude_;
   double rotated_longitude_;
   double rotated_angle_;

   double nsIncrement_;             ///< In degrees
   double weIncrement_;             ///< In degrees
   long nptsNS_;
   long nptsWE_;

   std::vector< Point > points_;     ///< storage of coordinate points

};

//-----------------------------------------------------------------------------

} // namespace grid
} // namespace eckit

#endif
