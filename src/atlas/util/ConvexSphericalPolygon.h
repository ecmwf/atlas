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

#include <vector>

#include "eckit/deprecated.h"

#include "atlas/util/Point.h"
#include "atlas/util/Polygon.h"
#include "atlas/runtime/Log.h"

#define DEBUG_OUTPUT 1

namespace atlas {
namespace util {

//------------------------------------------------------------------------------------------------------

class ConvexSphericalPolygon {
public:
    static constexpr int MAX_GRIDCELL_EDGES = 4;
    static constexpr int MAX_SIZE           = 2 * MAX_GRIDCELL_EDGES + 1;

public:
    class GreatCircleSegment {
    public:
        GreatCircleSegment(const PointXYZ& p1, const PointXYZ& p2): p1_(p1), p2_(p2), cross_(PointXYZ::cross(p1, p2)) {}

        // For a given segment, the "left" hemisphere is defined on the left of the segment when walking from first() to second()
        inline bool inLeftHemisphere(const PointXYZ& P, const double offset = 0., std::ostream* out = NULL) const {
#if DEBUG_OUTPUT
           if (out) {
            *out << " inLeftHemi: " <<PointXYZ::dot(cross(), P) << " >= " << offset << " ? " 
                   << (PointXYZ::dot(cross(), P) >= offset) << std::endl;
           }
#endif
            return (PointXYZ::dot(cross(), P) >= offset); // has to have = included
        }

        PointXYZ intersect(const GreatCircleSegment&, std::ostream* f = NULL, double pointsSameEPS = std::numeric_limits<double>::epsilon()) const;

        const PointXYZ& first() const { return p1_; }

        const PointXYZ& second() const { return p2_; }

        const PointXYZ& cross() const { return cross_; }

    private:
        const PointXYZ p1_;
        const PointXYZ p2_;
        const PointXYZ cross_;
    };

public:
    ConvexSphericalPolygon() = default;

    template <class Points>
    using contains_PointLonLat = std::is_same<typename std::decay<typename Points::value_type>::type, PointLonLat>;

    template <class Points, typename std::enable_if<contains_PointLonLat<Points>::value, void>::type* = nullptr>
    ConvexSphericalPolygon(const Points& points): ConvexSphericalPolygon(points.data(), points.size()) {}

    ConvexSphericalPolygon(const PointLonLat points[], size_t size);

    ConvexSphericalPolygon(const PointXYZ& p1, const PointXYZ& p2, const PointXYZ& p3):
        ConvexSphericalPolygon(std::array<PointXYZ, 3>{p1, p2, p3}.data(), 3) {}

    ConvexSphericalPolygon(const PointXYZ points[], size_t size);

    operator bool() const { return valid_; }

    size_t size() const { return size_; }

    double area() const {
        if (not computed_area_) {
            compute_centroid();
        }
        return area_;
    }

    const PointXYZ& centroid() const {
        if (not computed_centroid_) {
            compute_centroid();
        }
        return centroid_;
    }

    double radius() const {
        if (not computed_radius_) {
            radius_          = compute_radius();
            computed_radius_ = true;
        }
        return radius_;
    }

    ConvexSphericalPolygon intersect(const ConvexSphericalPolygon& pol, std::ostream* f = nullptr, double pointsEqualEPS = std::numeric_limits<double>::epsilon()) const;

    /*
   * @brief check if two spherical polygons area equal
   * @param[in] P given point in (x,y,z) coordinates
   * @return true if equal vertices
   */
    bool equals(const ConvexSphericalPolygon& plg, const double deg_prec = 1e-10) const;

    void print(std::ostream&) const;

    std::string json(int precision = 0) const;

    friend std::ostream& operator<<(std::ostream& out, const ConvexSphericalPolygon& p) {
        p.print(out);
        return out;
    }

    const PointXYZ& operator[](idx_t n) const { return sph_coords_[n]; }

    int next(const int index) const { return (index == size_ - 1) ? 0 : index + 1; };

private:
    struct SubTriangle {
        PointXYZ centroid;
        double area;
    };
    struct SubTriangles : public std::array<SubTriangle, MAX_SIZE> {
    public:
        size_t size() const { return size_; }
        size_t& size() { return size_; }
        double area() const;
        std::array<SubTriangle, MAX_SIZE>::const_iterator end() const { return data() + size_; }

    private:
        size_t size_{0};
    };

    void compute_centroid() const;

    double compute_radius() const;

    SubTriangles triangulate() const;

    void clip(const GreatCircleSegment&, std::ostream* f = nullptr, double pointsSameEPS = std::numeric_limits<double>::epsilon());

    /*
   * @return true:polygon is convex
   */
    bool validate();

private:
    std::array<PointXYZ, MAX_SIZE> sph_coords_;
    mutable PointXYZ centroid_;
    mutable double area_{0};
    mutable double radius_{0};
    size_t size_{0};
    bool valid_{false};
    mutable bool computed_centroid_{false};
    mutable bool computed_radius_{false};
    mutable bool computed_area_{false};
};

//------------------------------------------------------------------------------------------------------

}  // namespace util
}  // namespace atlas
