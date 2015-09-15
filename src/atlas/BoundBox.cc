/*
 * (C) Copyright 1996-2015 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "BoundBox.h"

using eckit::geometry::LLBoundBox2;
using eckit::MD5;

namespace atlas {

//------------------------------------------------------------------------------------------------------

BoundBox::BoundBox() : LLBoundBox2()
{}

BoundBox::BoundBox(double north, double south, double east, double west) : LLBoundBox2(north,south,east,west)
{}

BoundBox::BoundBox(const Point &min, const Point &max) : LLBoundBox2(min,max)
{}

BoundBox::BoundBox(Point pt, double x_delta, double y_delta) : LLBoundBox2(pt,x_delta,y_delta)
{}

void BoundBox::hash(eckit::MD5& md5) const {
  md5.add(min_.data(), min_.DIMS*sizeof(double));
  md5.add(max_.data(), max_.DIMS*sizeof(double));
}

//------------------------------------------------------------------------------------------------------

} // namespace atlas

