/*
 * (C) Copyright 1996-2015 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

/// @author Tiago Quintino
/// @date   May 2015

#ifndef atlas_BoundBox_h
#define atlas_BoundBox_h

#include "eckit/geometry/Point2.h"
#include "eckit/utils/MD5.h"

namespace atlas {

//------------------------------------------------------------------------------------------------------

class BoundBox : public eckit::geometry::LLBoundBox2 {

  typedef eckit::geometry::LLPoint2 Point;  ///< point type

public: // methods

  BoundBox();

  BoundBox( double north, double south, double east, double west );

  BoundBox( const Point& min, const Point& max );

  BoundBox( Point pt, double x_delta, double y_delta );

  void hash(eckit::MD5&) const;

};

//------------------------------------------------------------------------------------------------------

} // namespace atlas

#endif
