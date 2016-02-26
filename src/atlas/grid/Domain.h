/*
 * (C) Copyright 1996-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

/// @author Tiago Quintino
/// @author Baudouin Raoult
/// @date Jun 2015

#ifndef atlas_Domain_H
#define atlas_Domain_H

#include "eckit/exception/Exceptions.h"
#include "eckit/geometry/Point2.h"

namespace eckit {
  class MD5;
  namespace geometry {
    class LLPoint2;
  }
}

namespace atlas {
namespace grid {
    class BoundBox;
} }

namespace atlas {
namespace grid {

//------------------------------------------------------------------------------------------------------

/// If is periodic then east() - west() == 360
/// It has North Pole is north() == 90
/// It has South Pole is south() == -90

class Domain {

public:  // methods

    /// Construct from a BoundBox
    Domain(const atlas::grid::BoundBox&);

    /// East and West are reduced to the interval [0,360[
    Domain(double north, double west, double south, double east);

    ~Domain(); // make it virtual once is a virtual base

    /// Adds to the MD5 the information
    void hash(eckit::MD5&) const;

    /// checks if the point is contained in the domain
    bool contains( const eckit::geometry::LLPoint2& ) const;

    /// checks if the point is contained in the domain
    bool contains( double lon, double lat ) const;

    void print(std::ostream&) const;

    static Domain makeGlobal();

    double north() const { return north_; }
    double west() const { return west_; }
    double south() const { return south_; }
    double east() const { return east_; }

    bool global() const;

private:  // methods

  /// normalises the constructor input
  void normalise();

  /// normalises the longitude of a query point
  double normalise(double lon) const;

  friend std::ostream& operator<<(std::ostream& s, const Domain& p) {
      p.print(s);
      return s;
  }

private: // members

    double north_;
    double west_;
    double south_;
    double east_;

};

//------------------------------------------------------------------------------------------------------

} // namespace grid
} // namespace atlas

#endif
