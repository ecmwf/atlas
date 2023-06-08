/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <iomanip>
#include <iostream>

#include "eckit/geometry/Sphere.h"
#include "eckit/types/FloatCompare.h"

#include "atlas/runtime/Exception.h"
#include "atlas/runtime/Log.h"
#include "atlas/util/ConvexSphericalPolygon.h"
#include "atlas/util/CoordinateEnums.h"
#include "atlas/util/NormaliseLongitude.h"
#include "atlas/library/FloatingPointExceptions.h"

// #define DEBUG_OUTPUT 1

namespace atlas {
namespace util {

using GreatCircleSegment = ConvexSphericalPolygon::GreatCircleSegment;

namespace {

constexpr double EPS  = std::numeric_limits<double>::epsilon();
constexpr double EPS2 = EPS * EPS;

double distance2(const PointXYZ& p1, const PointXYZ& p2) {
    double dx = p2[0] - p1[0];
    double dy = p2[1] - p1[1];
    double dz = p2[2] - p1[2];
    return dx * dx + dy * dy + dz * dz;
}

double norm2(const PointXYZ& p) {
    return p[0] * p[0] + p[1] * p[1] + p[2] * p[2];
}

bool approx_eq(const double& v1, const double& v2) {
    return std::abs(v1 - v2) <= EPS;
}

bool approx_eq(const PointXYZ& v1, const PointXYZ& v2) {
    //return approx_eq( v1[0], v2[0], t ) && approx_eq( v1[1], v2[1], t ) && approx_eq( v1[2], v2[2], t );
    return distance2(v1, v2) <= EPS2;
}

void lonlat2xyz(const PointLonLat& lonlat, PointXYZ& xyz) {
    eckit::geometry::Sphere::convertSphericalToCartesian(1., lonlat, xyz);
}

void xyz2lonlat(const PointXYZ& xyz, PointLonLat& lonlat) {
    eckit::geometry::Sphere::convertCartesianToSpherical(1., xyz, lonlat);
}

PointLonLat xyz2lonlat(const PointXYZ& xyz) {
    PointLonLat lonlat;
    eckit::geometry::Sphere::convertCartesianToSpherical(1., xyz, lonlat);
    return lonlat;
}

#if 0
// NOTE: StackVector is not used
template <typename T, size_t MAX_SIZE>
struct StackVector {
private:
    using Wrapped = std::array<T,MAX_SIZE>;

public:
    using reference = typename Wrapped::reference;
    StackVector()   = default;
    StackVector(size_t size): size_(size) {}
#if ATLAS_VECTOR_BOUNDS_CHECKING
    reference operator[](size_t i) {
        if (i >= size_) {
            throw_OutOfRange(i, size_);
        }
        return wrapped_[i];
    }
#else
    reference operator[](size_t i) { return wrapped_[i]; }
#endif
    void push_back(const T& value) {
        wrapped_[size_++] = value;
        ATLAS_ASSERT(size_ < MAX_SIZE);
    }
    size_t size() const { return size_; }

private:
    size_t size_{0};
    Wrapped wrapped_;
};
#endif

}  // namespace

//------------------------------------------------------------------------------------------------------

ConvexSphericalPolygon::ConvexSphericalPolygon(const PointLonLat points[], size_t size): size_{size} {
    ATLAS_ASSERT(size_ > 2, "Polygon must have at least 3 points");
    ATLAS_ASSERT(size_ < MAX_SIZE, "Number of polygon points exceeds compile time MAX_SIZE");
    lonlat2xyz(points[0], sph_coords_[0]);
    size_t isp = 1;
    for (size_t i = 1; i < size_ - 1; ++i) {
        lonlat2xyz(points[i], sph_coords_[isp]);
//        Log::info() << " d : " << PointXYZ::distance(sph_coords_[isp], sph_coords_[isp - 1]) << std::endl;
        if (approx_eq(sph_coords_[isp], sph_coords_[isp - 1])) {
            continue;
        }
        ++isp;
    }
    lonlat2xyz(points[size_ - 1], sph_coords_[isp]);
    if (approx_eq(sph_coords_[isp], sph_coords_[0]) or approx_eq(sph_coords_[isp], sph_coords_[isp - 1])) {
    }
    else {
        ++isp;
    }
    size_  = isp;
    valid_ = size_ > 2;
    if (valid_) {
        ATLAS_ASSERT(validate());
    }
}

ConvexSphericalPolygon::ConvexSphericalPolygon(const PointXYZ points[], size_t size): size_{size} {
    ATLAS_ASSERT(size_ > 2, "Polygon must have at least 3 points");
    ATLAS_ASSERT(size_ < MAX_SIZE, "Number of polygon points exceeds compile time MAX_SIZE");
    sph_coords_[0] = points[0];
    size_t isp     = 1;
    for (size_t i = 1; i < size_ - 1; ++i) {
        sph_coords_[isp] = points[i];
        if (approx_eq(sph_coords_[isp], sph_coords_[isp - 1])) {
            continue;
        }
        ++isp;
    }
    sph_coords_[isp] = points[size_ - 1];
    if (approx_eq(sph_coords_[isp], sph_coords_[0]) or approx_eq(sph_coords_[isp], sph_coords_[isp - 1])) {
    }
    else {
        ++isp;
    }
    size_  = isp;
    valid_ = size_ > 2;
    if (valid_) {
        ATLAS_ASSERT(validate());
    }
}

void ConvexSphericalPolygon::compute_centroid() const {
    const auto triangles = triangulate();

    area_          = triangles.area();
    computed_area_ = true;

    // Compute centroid based on triangles rather than on vertices of polygon
    centroid_ = PointXYZ{0., 0., 0.};
    if (area_ > 0.) {
        for (auto& triangle : triangles) {
            for (size_t i = 0; i < 3; ++i) {
                centroid_[i] += triangle.centroid[i] * triangle.area;
            }
        }
        centroid_ /= PointXYZ::norm(centroid_);
    }
    computed_centroid_ = true;
}

bool ConvexSphericalPolygon::validate() {
    if (valid_) {
        for (int i = 0; i < size(); i++) {
            int ni                = next(i);
            int nni               = next(ni);
            const PointXYZ& P     = sph_coords_[i];
            const PointXYZ& nextP = sph_coords_[ni];
            ATLAS_ASSERT(std::abs(PointXYZ::dot(P, P) - 1.) < 10 * EPS);
            ATLAS_ASSERT(not approx_eq(P, PointXYZ::mul(nextP, -1.)));
            valid_ = valid_ && GreatCircleSegment{P, nextP}.inLeftHemisphere(sph_coords_[nni], -0.5*EPS);
        }
    }
    return valid_;
}

bool ConvexSphericalPolygon::equals(const ConvexSphericalPolygon& plg, const double deg_prec) const {
    if (size_ == 0 and plg.size_ == 0) {
        return true;
    }
    if (not valid_) {
        Log::info() << " ConvexSphericalPolygon::equals : this polygon is not valid\n";
        return false;
    }
    if (not plg.valid_) {
        Log::info() << " ConvexSphericalPolygon::equals : other polygon passed as argument is not valid\n";
        return false;
    }
    if (size_ != plg.size()) {
        Log::info() << " ConvexSphericalPolygon::equals : incompatible sizes: " << size_ << " != " << plg.size() << "\n";
        return false;
    }
    int offset       = 0;
    const double le  = 2. * std::sin(M_PI * deg_prec / 360.);
    const double le2 = le * le;
    for (; offset < size_; ++offset) {
        if (distance2(plg.sph_coords_[0], sph_coords_[offset]) < le2) {
            break;
        }
    }
    if (offset == size_) {
        Log::info() << "ConvexSphericalPolygon::equals : no point equal\n";
        return false;
    }

    for (int j = 0; j < size_; j++) {
        int idx    = (offset + j) % size_;
        auto dist2 = distance2(plg.sph_coords_[j], sph_coords_[idx]);
        if (dist2 > le2) {
            Log::info() << "  ConvexSphericalPolygon::equals : point distance " << std::sqrt(dist2) << " < " << le2 << "(= "<< deg_prec << " deg)"  << "\n";
            return false;
        }
    }
    return true;
}

// cf. Folke Eriksson, "On the Measure of Solid Angles", Mathematics Magazine, Vol. 63, No. 3, pp. 184-187 (1990)
ConvexSphericalPolygon::SubTriangles ConvexSphericalPolygon::triangulate() const {
    SubTriangles triangles;
    if (size_ < 3) {
        return triangles;
    }
    size_t itri{0};
    const PointXYZ& a = sph_coords_[0];
    for (size_t i = 1; i < size_ - 1; i++) {
        const PointXYZ& b   = sph_coords_[i];
        const PointXYZ& c   = sph_coords_[i + 1];
        triangles[itri].centroid = PointXYZ::normalize(a + b + c);
//        if (PointXYZ::distance(a, b) + PointXYZ::distance(b, c) < 1e-10) {
//            triangles[itri].area     = 0.5 * PointXYZ::norm(PointXYZ::cross(b - a, c - b));
//        }
//        else {
            auto abc            = PointXYZ::dot(a, b) + PointXYZ::dot(b, c) + PointXYZ::dot(c, a);
            auto a_bc           = PointXYZ::dot(a, PointXYZ::cross(b, c));
            triangles[itri].area     = 2. * std::atan(std::abs(a_bc) / (1. + abc));
//        }
        ++itri;
    }
    triangles.size() = itri;
    return triangles;

}


double ConvexSphericalPolygon::SubTriangles::area() const {
    double area = 0.;
    for (auto& triangle : *this) {
        area += triangle.area;
    }
    return area;
}

void ConvexSphericalPolygon::clip(const GreatCircleSegment& great_circle, std::ostream* out, double pointsSameEPS) {
    ATLAS_ASSERT(valid_);
    ATLAS_ASSERT(not approx_eq(great_circle.first(), great_circle.second()));
    std::vector<PointXYZ> clipped_sph_coords;
    clipped_sph_coords.reserve(ConvexSphericalPolygon::MAX_SIZE);
    auto invalidate_this_polygon = [&]() {
        size_  = 0;
        valid_ = false;
        area_  = 0.;
    };
#if DEBUG_OUTPUT
    int add_point_num = 0;
#endif
    bool first_in = great_circle.inLeftHemisphere(sph_coords_[0], -1.5 * EPS, out);
    for (int i = 0; i < size_; i++) {
        int in = (i+1) % size_;
        bool second_in = great_circle.inLeftHemisphere(sph_coords_[in], -1.5 * EPS, out);
#if DEBUG_OUTPUT
        if (out) {
            out->precision(18);
            *out << " ** first: " << xyz2lonlat(sph_coords_[i]) << ", in ? " << first_in << std::endl;
            *out << "    second: " << xyz2lonlat(sph_coords_[in]) << ", in ? " << second_in << std::endl;
        }
#endif
        if (first_in and second_in) {
            clipped_sph_coords.emplace_back(sph_coords_[in]);
#if DEBUG_OUTPUT
            if (out) {
                *out << " ((" << ++add_point_num << ")) both_in add second: " << xyz2lonlat(sph_coords_[in]) << std::endl;
            }
#endif
        }
        else if (not first_in and not second_in) {
            // continue to update first_in
        }
        else {
            const GreatCircleSegment segment(sph_coords_[i], sph_coords_[in]);
            PointXYZ ip = great_circle.intersect(segment, out, pointsSameEPS);
#if DEBUG_OUTPUT
            if (out) {
                *out << " ip : " << xyz2lonlat(ip) << std::endl;
            }
#endif
            if (ip[0] == 1 and ip[1] == 1 and ip[2] == 1) {
                // consider the segments parallel
#if DEBUG_OUTPUT
                if (out) {
                    *out << " ((" << ++add_point_num << ")) ip=(1,1,1) add second: "
                        << xyz2lonlat(sph_coords_[in]) << std::endl;
                }
#endif
                clipped_sph_coords.emplace_back(sph_coords_[in]);
                first_in = second_in;
                continue;
            }
            if (second_in) {
                int inn = (in+1) % size_;
                const GreatCircleSegment segment_n(sph_coords_[in], sph_coords_[inn]);
                if (segment.inLeftHemisphere(ip, -1.5 * EPS, out) and
                    segment_n.inLeftHemisphere(ip, -1.5 * EPS, out) and
                    (PointXYZ::distance(ip, sph_coords_[in]) > pointsSameEPS)) {
                    clipped_sph_coords.emplace_back(ip);
#if DEBUG_OUTPUT
                    if (out) {
                        *out << " ((" << ++add_point_num << ")) second_in add ip: " << xyz2lonlat(ip) << std::endl;
                    }
#endif
                }
                clipped_sph_coords.emplace_back(sph_coords_[in]);
#if DEBUG_OUTPUT
                if (out) {
                    *out << " ((" << ++add_point_num << ")) second_in add second: " << xyz2lonlat(sph_coords_[in]) << std::endl;
                }
#endif
            }
            else {
                if (PointXYZ::distance(ip, sph_coords_[i]) > pointsSameEPS) {
                    clipped_sph_coords.emplace_back(ip);
#if DEBUG_OUTPUT
                    if (out) {
                        *out << " ((" << ++add_point_num << ")) first_in add ip: " << xyz2lonlat(ip) << std::endl;
                    }
#endif
                }
            }
        }
        first_in = second_in;
    }

    // Update polygon
    {
        size_ = clipped_sph_coords.size();
        if (size_ < 3) {
            invalidate_this_polygon();
        }
        else {
            for (size_t i = 0; i < size_; i++) {
                sph_coords_[i] = clipped_sph_coords[i];
            }
        }
    }
}

// intersect a polygon with this polygon
// @param[in] pol clipping polygon
// @param[out] intersecting polygon
ConvexSphericalPolygon ConvexSphericalPolygon::intersect(const ConvexSphericalPolygon& plg, std::ostream* out, double pointsSameEPS) const {

    bool fpe_disabled = atlas::library::disable_floating_point_exception(FE_INVALID);
    auto restore_fpe = [fpe_disabled] {
        if (fpe_disabled) {
            atlas::library::enable_floating_point_exception(FE_INVALID);
        }
    };

    // the larger area polygon is the intersector
    ConvexSphericalPolygon intersection;
    ConvexSphericalPolygon intersector;
    std::string intor_id = "P1";
    std::string inted_id = "P2";

    bool this_intersector = (area() > plg.area());
    if (approx_eq(area(), plg.area())) {
        PointXYZ dc = centroid() - plg.centroid();
        if (dc[0] > 0. or (dc[0] == 0. and (dc[1] > 0. or (dc[1] == 0. and dc[2] > 0.)))) {
            this_intersector = true;
        }
    }
    if (this_intersector) {
        intersector = *this;
        intersection = plg;
    }
    else {
        intor_id = "P2";
        inted_id = "P1";
        intersector = plg;
        intersection = *this;
    }
#if DEBUG_OUTPUT
    if (out) {
        *out << inted_id << " : ";
        print(*out);
        *out << std::endl << intor_id << " : ";
        plg.print(*out);
        *out << std::endl;
    }
#endif
    if (intersection.valid_) {
        for (size_t i = 0; i < intersector.size_; i++) {
            const PointXYZ& s1 = intersector.sph_coords_[i];
            const PointXYZ& s2 = intersector.sph_coords_[(i != intersector.size_ - 1) ? i + 1 : 0];
#if DEBUG_OUTPUT
            if (out) {
                *out << std::endl << "Clip with [" << intor_id << "_" << i << ", " << intor_id << "_" 
                    << (i+1) % intersector.size_ << "]" << std::endl;
            }
#endif
            intersection.clip(GreatCircleSegment(s1, s2), out, pointsSameEPS);
            if (not intersection.valid_) {
                restore_fpe();
                return intersection;
            }
        }
    }
    intersection.computed_area_     = false;
    intersection.computed_radius_   = false;
    intersection.computed_centroid_ = false;
    restore_fpe();
    return intersection;
}

void ConvexSphericalPolygon::print(std::ostream& out) const {
    out << "{";
    for (size_t i = 0; i < size(); ++i) {
        if (i > 0) {
            out << ",";
        }
        PointLonLat ip_ll;
        xyz2lonlat(sph_coords_[i], ip_ll);
        out << ip_ll;
    }
    out << "}";
}

std::string ConvexSphericalPolygon::json(int precision) const {
    std::stringstream ss;
    if( precision ) {
        ss.precision(16);
    }
    ss << "[";
    for (size_t i = 0; i < size(); ++i) {
        if (i > 0) {
            ss << ",";
        }
        PointLonLat ip_ll;
        xyz2lonlat(sph_coords_[i], ip_ll);
        ss << "[" << ip_ll.lon() << "," << ip_ll.lat() << "]";
    }
    ss << "]";
    return ss.str();
}


double ConvexSphericalPolygon::compute_radius() const {
    double radius{0.};
    if (valid_) {
        if (not computed_centroid_) {
            compute_centroid();
        }
        for (size_t i = 0; i < size_; ++i) {
            radius = std::max(radius, PointXYZ::distance(sph_coords_[i], centroid_));
        }
    }
    return radius;
}


// 'this' great circle's intersection with the segment 'p': [p.first(), p.second())
PointXYZ ConvexSphericalPolygon::GreatCircleSegment::intersect(const GreatCircleSegment& p, std::ostream* out, double pointsSameEPS) const {
    const auto& s = *this;
    PointXYZ sp   = PointXYZ::cross(s.cross(), p.cross());

    double sp_norm = PointXYZ::norm(sp);
    bool gcircles_distinct = (sp_norm > EPS);
#if DEBUG_OUTPUT
    if (out) {
        *out << " Great circles distinct ? " << sp_norm << " > " << EPS << " ? " 
            << gcircles_distinct << std::endl;
    }
#endif
    if (gcircles_distinct) {
        sp /= sp_norm;
        auto sp_2 = sp * -1.;
        double d = distance2(p.first(), p.second());
        double d1 = std::max(distance2(sp, p.first()), distance2(sp, p.second()));
        double d2 = std::max(distance2(sp_2, p.first()), distance2(sp_2, p.second()));
        if (d1 < d2) {
            return sp;
        }
        else {
            return sp_2;
        }
    }
    else {
        return PointXYZ(1, 1, 1);
    }
}


//------------------------------------------------------------------------------------------------------

}  // namespace util
}  // namespace atlas
