/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#pragma once

/// @file Point.h
///
/// This file contains classes and functions working on points.
/// The Point classes are inherited from eckit::geometry::Point2
/// or eckit::geometry::Point3.

#include <array>
#include <initializer_list>

#include "eckit/geometry/Point2.h"
#include "eckit/geometry/Point3.h"

#include "atlas/util/Earth.h"

namespace atlas {

using Point2 = eckit::geometry::Point2;
using Point3 = eckit::geometry::Point3;

inline bool operator==(const Point2& p1, const Point2& p2) {
    return eckit::geometry::points_equal(p1, p2);
}
inline bool operator!=(const Point2& p1, const Point2& p2) {
    return !eckit::geometry::points_equal(p1, p2);
}

inline bool operator==(const Point3& p1, const Point3& p2) {
    return eckit::geometry::points_equal(p1, p2);
}
inline bool operator!=(const Point3& p1, const Point3& p2) {
    return !eckit::geometry::points_equal(p1, p2);
}


/// @brief  Point in arbitrary XY-coordinate system
class PointXY : public eckit::geometry::Point2 {
    using array_t = std::array<double, 2>;

public:
    using Point2::Point2;

    PointXY(): Point2() {}

    // Allow initialization through PointXY xy = {0,0};
    PointXY(std::initializer_list<double> list): PointXY(list.begin()) {}
    PointXY(const array_t& arr): PointXY(arr.data()) {}

    double x() const { return x_[0]; }
    double y() const { return x_[1]; }
    double& x() { return x_[0]; }
    double& y() { return x_[1]; }

    using Point2::assign;

    void assign(double x, double y) {
        x_[0] = x;
        x_[1] = y;
    }
};

/// @brief Point in arbitrary XYZ-coordinate system
class PointXYZ : public eckit::geometry::Point3 {
    using array_t = std::array<double, 3>;

    PointXYZ(double, double) { /* No automatic converion allowed, otherwise
                                inherited from Point3 */
    }

public:
    using Point3::Point3;
    using Point3::x;

    PointXYZ(): Point3() {}

    // Allow initialization through PointXYZ xyz = {0,0,0};
    PointXYZ(std::initializer_list<double> list): PointXYZ(list.begin()) {}

    PointXYZ(const array_t& arr): PointXYZ(arr.data()) {}

    double x() const { return x_[0]; }
    double y() const { return x_[1]; }
    double z() const { return x_[2]; }
    double& x() { return x_[0]; }
    double& y() { return x_[1]; }
    double& z() { return x_[2]; }

    using Point3::assign;

    void assign(double x, double y, double z) {
        x_[0] = x;
        x_[1] = y;
        x_[2] = z;
    }

    PointXYZ& operator*=(double m) {
        x_[0] *= m;
        x_[1] *= m;
        x_[2] *= m;
        return *this;
    }

    PointXYZ& operator/=(double m) {
        x_[0] /= m;
        x_[1] /= m;
        x_[2] /= m;
        return *this;
    }
};

/// @brief Point in longitude-latitude coordinate system.
/// This class does *not* normalise the longitude by default,
/// but contains a normalise function.
class PointLonLat : public eckit::geometry::Point2 {
    using array_t = std::array<double, 2>;

public:
    using Point2::Point2;

    PointLonLat(): Point2() {}

    // Allow initialization through PointXY lonlat = {0,0};
    PointLonLat(std::initializer_list<double> list): PointLonLat(list.begin()) {}

    PointLonLat(const array_t& arr): PointLonLat(arr.data()) {}

    double lon() const { return x_[0]; }
    double lat() const { return x_[1]; }
    double& lon() { return x_[0]; }
    double& lat() { return x_[1]; }

    using Point2::assign;

    void assign(double lon, double lat) {
        x_[0] = lon;
        x_[1] = lat;
    }

    PointLonLat& operator*=(double a) {
        x_[0] *= a;
        x_[1] *= a;
        return *this;
    }

    void normalise();

    void normalise(double west);

    void normalise(double west, double east);
};

}  // namespace atlas

#ifndef DOXYGEN_SHOULD_SKIP_THIS
namespace eckit {
// When including "eckit/types/Types.h" we get an overload for
//    std::ostream& operator<<( std::ostream&, std::vector<T> )
// The default VectorPrintSelector is however the [VectorPrintContracted],
// which does not compile when [T={PointXY,PointLonLat,PointXYZ}]
// Following changes the default for these types to [VectorPrintSimple]
class VectorPrintSimple;
template <typename T>
struct VectorPrintSelector;
template <>
struct VectorPrintSelector<atlas::Point2> {
    typedef VectorPrintSimple selector;
};
template <>
struct VectorPrintSelector<atlas::PointXY> {
    typedef VectorPrintSimple selector;
};
template <>
struct VectorPrintSelector<atlas::PointLonLat> {
    typedef VectorPrintSimple selector;
};
template <>
struct VectorPrintSelector<atlas::Point3> {
    typedef VectorPrintSimple selector;
};
template <>
struct VectorPrintSelector<atlas::PointXYZ> {
    typedef VectorPrintSimple selector;
};
}  // namespace eckit
#endif
