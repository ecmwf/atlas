/*
 * (C) Copyright 1996-2014 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef atlas_geometry_Intersect_h
#define atlas_geometry_Intersect_h

#include <iosfwd>

namespace atlas {
namespace geometry {

//----------------------------------------------------------------------------------------------------------------------

const double parametricEpsilon = 1e-12; ///< Epsilon used to compare weights and u,v's

/// Intersection data structure

struct Intersect {

  double u;
  double v;
  double t;

  Intersect();

  operator bool() const { return success_; }

  Intersect& success(bool s){ success_ = s; return *this; }

  void print(std::ostream& s) const;

  friend std::ostream& operator<<(std::ostream& s, const Intersect& p) {
    p.print(s);
    return s;
  }

private:

  bool success_;

};

//----------------------------------------------------------------------------------------------------------------------

}  // namespace geometry
}  // namespace atlas

#endif
