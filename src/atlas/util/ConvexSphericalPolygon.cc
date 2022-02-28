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

#define DEBUG_OUTPUT 0
#define DEBUG_OUTPUT_DETAIL 0

namespace atlas {
namespace util {

using GreatCircleSegment = ConvexSphericalPolygon::GreatCircleSegment;

namespace {

constexpr double EPS  = std::numeric_limits<double>::epsilon();
constexpr double EPS2 = EPS * EPS;
constexpr double TOL  = 1.e4 * EPS;  // two points considered "same"
constexpr double TOL2 = TOL * TOL;

enum IntersectionType
{
    NO_INTERSECT = -100,
    OVERLAP
};

double distance2(const PointXYZ& p1, const PointXYZ& p2) {
    double dx = p2[0] - p1[0];
    double dy = p2[1] - p1[1];
    double dz = p2[2] - p1[2];
    return dx * dx + dy * dy + dz * dz;
}

double norm2(const PointXYZ& p) {
    return p[0] * p[0] + p[1] * p[1] + p[2] * p[2];
}

bool approx_eq(const double& v1, const double& v2, const double& tol) {
    return std::abs(v1 - v2) < tol;
}

bool approx_eq(const PointXYZ& v1, const PointXYZ& v2, const double& tol) {
    //return approx_eq( v1[0], v2[0], t ) && approx_eq( v1[1], v2[1], t ) && approx_eq( v1[2], v2[2], t );
    return distance2(v1, v2) < tol * tol;
}

bool approx_eq_null(const PointXYZ& v1, const double& tol) {
    //return approx_eq( v1[0], 0., t ) && approx_eq( v1[1], 0., t ) && approx_eq( v1[2], 0., t );
    double n = v1[0] * v1[0] + v1[1] * v1[1] + v1[2] * v1[2];
    return n < tol * tol;
}

void lonlat2xyz(const PointLonLat& lonlat, PointXYZ& xyz) {
    eckit::geometry::Sphere::convertSphericalToCartesian(1., lonlat, xyz);
}

void xyz2lonlat(const PointXYZ xyz, PointLonLat& lonlat) {
    eckit::geometry::Sphere::convertCartesianToSpherical(1., xyz, lonlat);
}

double norm_max(const PointXYZ& p, const PointXYZ& q) {
    double n01 = std::max(std::abs(p[0] - q[0]), std::abs(p[1] - q[1]));
    return std::max(n01, std::abs(p[2] - q[2]));
}

template <typename T>
struct StackVector {
private:
    using Wrapped = std::array<T, ConvexSphericalPolygon::MAX_SIZE>;

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
        ATLAS_ASSERT(size_ < ConvexSphericalPolygon::MAX_SIZE);
    }
    size_t size() const { return size_; }

private:
    size_t size_{0};
    Wrapped wrapped_;
};

struct PolygonEdgeIntersection {
    static constexpr int BEGIN  = 1;
    static constexpr int END    = 2;
    static constexpr int INSIDE = 3;

    PolygonEdgeIntersection(const ConvexSphericalPolygon& polygon, int edge_index, const PointXYZ& point) {
        auto matches = [](const PointXYZ& p1, const PointXYZ& p2) {
            return (distance2(p1, p2) < 1e-16);
            // We would like this to be TOL2 instead, but that gives bad results
        };

        ATLAS_ASSERT(edge_index >= 0);
        ATLAS_ASSERT(edge_index < polygon.size());

        const int node_index = edge_index;

        if (matches(point, polygon[node_index])) {
            location = BEGIN;
        }
        else if (matches(point, polygon[polygon.next(node_index)])) {
            location = END;
        }
        else {
            location = INSIDE;
        }
    }
    bool isPointAtBegin() const { return location == BEGIN; }
    bool isPointAtEnd() const { return location == END; }
    bool isPointInside() const { return location == INSIDE; }
    int location;
};

}  // namespace

//------------------------------------------------------------------------------------------------------

ConvexSphericalPolygon::ConvexSphericalPolygon(const PointLonLat points[], size_t size): size_{size} {
    ATLAS_ASSERT(size_ > 2, "Polygon must have at least 3 points");
    ATLAS_ASSERT(size_ < MAX_SIZE, "Number of polygon points exceeds compile time MAX_SIZE");
    lonlat2xyz(points[0], sph_coords_[0]);
    size_t isp = 1;
    for (size_t i = 1; i < size_ - 1; ++i) {
        lonlat2xyz(points[i], sph_coords_[isp]);
        if (approx_eq(sph_coords_[isp], sph_coords_[isp - 1], TOL)) {
            continue;
        }
        ++isp;
    }
    lonlat2xyz(points[size_ - 1], sph_coords_[isp]);
    if (approx_eq(sph_coords_[isp], sph_coords_[0], TOL) or approx_eq(sph_coords_[isp], sph_coords_[isp - 1], TOL)) {
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
        if (approx_eq(sph_coords_[isp], sph_coords_[isp - 1], TOL)) {
            continue;
        }
        ++isp;
    }
    sph_coords_[isp] = points[size_ - 1];
    if (approx_eq(sph_coords_[isp], sph_coords_[0], TOL) or approx_eq(sph_coords_[isp], sph_coords_[isp - 1], TOL)) {
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


void ConvexSphericalPolygon::simplify() {
    ATLAS_ASSERT(size_ < MAX_SIZE);
    if (size_ < 3) {
        size_  = 0;
        valid_ = false;
        return;
    }
    idx_t isp = 0;
    idx_t i   = 0;
    idx_t j;
    idx_t k;
    bool search_3pts = true;
    auto& points     = sph_coords_;
    for (; i < size_ && search_3pts; ++i) {
        const PointXYZ& P0 = points[i];
        for (j = i + 1; j < size_ && search_3pts; ++j) {
            const PointXYZ& P1 = points[j];
            if (approx_eq(P0, P1, 1.e-10)) {
                continue;
            }
            for (k = j + 1; k < size_ && search_3pts; ++k) {
                const PointXYZ& P2 = points[k];
                if (approx_eq(P1, P2, TOL) or approx_eq(P0, P2, TOL)) {
                    continue;
                }
                if (GreatCircleSegment{P0, P1}.inLeftHemisphere(P2, -EPS)) {
                    sph_coords_[isp++] = P0;
                    sph_coords_[isp++] = P1;
                    sph_coords_[isp++] = P2;
                    search_3pts        = false;
                }
            }
        }
    }
    if (search_3pts) {
        valid_ = false;
        size_  = 0;
        return;
    }
    for (; k < size_ - 1; ++k) {
        if (approx_eq(points[k], sph_coords_[isp - 1], TOL) or
            (not GreatCircleSegment{sph_coords_[isp - 2], sph_coords_[isp - 1]}.inLeftHemisphere(points[k], -EPS))) {
            continue;
        }
        sph_coords_[isp] = points[k];
        isp++;
    }
    const PointXYZ& Pl2 = sph_coords_[isp - 2];
    const PointXYZ& Pl1 = sph_coords_[isp - 1];
    const PointXYZ& P0  = sph_coords_[0];
    const PointXYZ& P   = points[size_ - 1];
    if ((not approx_eq(P, P0, EPS)) and (not approx_eq(P, Pl1, EPS)) and
        GreatCircleSegment{Pl2, Pl1}.inLeftHemisphere(P, -EPS) and
        GreatCircleSegment{Pl1, P}.inLeftHemisphere(P0, -EPS)) {
        sph_coords_[isp] = P;
        ++isp;
    }
    size_  = isp;
    valid_ = size_ > 2;

    computed_area_     = false;
    computed_radius_   = false;
    computed_centroid_ = false;
}

void ConvexSphericalPolygon::compute_centroid() const {
    const auto triangles = triangulate(radius());

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
            ATLAS_ASSERT(std::abs(PointXYZ::dot(P, P) - 1.) < 10. * EPS);
            ATLAS_ASSERT(not approx_eq(P, PointXYZ::mul(nextP, -1.), TOL));
            valid_ = valid_ && GreatCircleSegment{P, nextP}.inLeftHemisphere(sph_coords_[nni], -EPS);
        }
    }
    return valid_;
}

bool ConvexSphericalPolygon::equals(const ConvexSphericalPolygon& plg, const double deg_prec) const {
    if ((not plg.valid_) || (not valid_) || size_ != plg.size()) {
        Log::info() << " ConvexSphericalPolygon::equals == not compatible\n";
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
        Log::info() << "ConvexSphericalPolygon::equals == no point equal\n";
        return false;
    }

    for (int j = 0; j < size_; j++) {
        int idx    = (offset + j) % size_;
        auto dist2 = distance2(plg.sph_coords_[j], sph_coords_[idx]);
        if (dist2 > le2) {
            Log::info() << "  ConvexSphericalPolygon::equals == point distance " << std::sqrt(dist2) << "\n";
            return false;
        }
    }
    return true;
}

// note: unit sphere!
// I. Todhunter (1886), Paragr. 99
ConvexSphericalPolygon::SubTriangles ConvexSphericalPolygon::triangulate(const double cell_radius) const {
    SubTriangles triangles;
    if (size_ < 3) {
        return triangles;
    }

    size_t itri{0};
    if (cell_radius < 1.e-6) {  // plane area
        for (int i = 1; i < size_ - 1; i++) {
            const PointXYZ pl        = sph_coords_[i] - sph_coords_[0];
            const PointXYZ pr        = sph_coords_[i + 1] - sph_coords_[0];
            triangles[itri].centroid = PointXYZ::normalize(sph_coords_[0] + sph_coords_[i] + sph_coords_[i + 1]);
            triangles[itri].area     = 0.5 * PointXYZ::norm(PointXYZ::cross(pl, pr));
            ++itri;
        }
    }
    else {  // spherical area
        const PointXYZ& a = sph_coords_[0];
        for (size_t i = 1; i < size_ - 1; i++) {
            const PointXYZ& b    = sph_coords_[i];
            const PointXYZ& c    = sph_coords_[i + 1];
            auto ab              = PointXYZ::cross(a, b);
            auto bc              = PointXYZ::cross(b, c);
            auto ca              = PointXYZ::cross(c, a);
            const double ab_norm = PointXYZ::norm(ab);
            const double bc_norm = PointXYZ::norm(bc);
            const double ca_norm = PointXYZ::norm(ca);
            if (ab_norm < EPS or bc_norm < EPS or ca_norm < EPS) {
                continue;
            }
            double abc = -PointXYZ::dot(ab, bc) / (ab_norm * bc_norm);
            double bca = -PointXYZ::dot(bc, ca) / (bc_norm * ca_norm);
            double cab = -PointXYZ::dot(ca, ab) / (ca_norm * ab_norm);
            if (abc <= -1.) {
                abc = M_PI;
            }
            else if (abc < 1.) {
                abc = std::acos(abc);
            }
            else {
                abc = 0.;
            }
            if (bca <= -1.) {
                bca = M_PI;
            }
            else if (bca < 1.) {
                bca = std::acos(bca);
            }
            else {
                bca = 0.;
            }
            if (cab <= -1.) {
                cab = M_PI;
            }
            else if (cab < 1.) {
                cab = std::acos(cab);
            }
            else {
                cab = 0.;
            }
            triangles[itri].centroid = PointXYZ::normalize(a + b + c);
            triangles[itri].area     = abc + bca + cab - M_PI;
            ++itri;
        }
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

// @return lowest point id of this polygon's segment intersecting [s1,s2))
int ConvexSphericalPolygon::intersect(const int start, const GreatCircleSegment& s, PointXYZ& I) const {
    for (int i = start; i < size_; i++) {
        const int id0 = i;
        const int id1 = next(i);
        const GreatCircleSegment p(sph_coords_[id0], sph_coords_[id1]);
        I = s.intersect(p);
        if (I[0] == 0 && I[1] == 0 && I[2] == 0) {
            // intersection not on [p1,p2)
            continue;
        }
        if (I[0] == 1 && I[1] == 1) {
            return OVERLAP;
        }
        return id0;
    }
    return NO_INTERSECT;
}


void ConvexSphericalPolygon::clip(const GreatCircleSegment& great_circle) {
    ATLAS_ASSERT(valid_);
    ATLAS_ASSERT(distance2(great_circle.first(), great_circle.second()) > TOL2);

    auto invalidate_this_polygon = [&]() {
        size_  = 0;
        valid_ = false;
        area_  = 0.;
    };

    // Count and mark all vertices to be possibly considered in clipped polygon
    StackVector<int> vertex_in(size_);
    int num_vertices_in = 0;
    for (int i = 0; i < size_; i++) {
        vertex_in[i] = great_circle.inLeftHemisphere(sph_coords_[i], -EPS);
        num_vertices_in += vertex_in[i] ? 1 : 0;
    }

    PointXYZ i1;
    const int f1                               = intersect(0, great_circle, i1);
    const bool segment_only_touches_last_point = (f1 == size_ - 1);
    if (f1 == OVERLAP || f1 == NO_INTERSECT || segment_only_touches_last_point) {
        if (num_vertices_in < 3) {
            invalidate_this_polygon();
        }
        return;
    }
    PolygonEdgeIntersection intersection_1(*this, f1, i1);

    PointXYZ i2;  // second intersection point
    auto start2  = [&]() { return f1 + 1 + (intersection_1.isPointAtEnd() ? 1 : 0); };
    const int f2 = intersect(start2(), great_circle, i2);
    if (f2 == OVERLAP || f2 == NO_INTERSECT) {
        if (num_vertices_in < 3) {
            invalidate_this_polygon();
        }
        return;
    }
    PolygonEdgeIntersection intersection_2(*this, f2, i2);

    // Create new vector of clipped coordinates
    StackVector<PointXYZ> clipped_sph_coords;
    {
        auto keep_vertex  = [&](int index) { clipped_sph_coords.push_back(sph_coords_[index]); };
        auto insert_point = [&](const PointXYZ& p) { clipped_sph_coords.push_back(p); };

        for (int i = 0; i <= f1; i++) {
            if (vertex_in[i]) {
                keep_vertex(i);
            }
        }
        if ((not vertex_in[f1] and intersection_1.isPointAtBegin()) or
            (not vertex_in[next(f1)] and intersection_1.isPointAtEnd()) or intersection_1.isPointInside()) {
            insert_point(i1);
        }
        for (int i = f1 + 1; i <= f2; i++) {
            if (vertex_in[i]) {
                keep_vertex(i);
            }
        }
        if ((not vertex_in[f2] and intersection_2.isPointAtBegin()) or
            (not vertex_in[next(f2)] and intersection_2.isPointAtEnd()) or intersection_2.isPointInside()) {
            insert_point(i2);
        }
        for (int i = f2 + 1; i < size_; i++) {
            if (vertex_in[i]) {
                keep_vertex(i);
            }
        }
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
ConvexSphericalPolygon ConvexSphericalPolygon::intersect(const ConvexSphericalPolygon& plg) const {
    ConvexSphericalPolygon intersection = *this;
    if (intersection.valid_) {
        for (size_t i = 0; i < plg.size_; i++) {
            const PointXYZ& s1 = plg.sph_coords_[i];
            const PointXYZ& s2 = plg.sph_coords_[(i != plg.size_ - 1) ? i + 1 : 0];
            intersection.clip(GreatCircleSegment(s1, s2));
            if (not intersection.valid_) {
                return intersection;
            }
        }
    }
    intersection.simplify();
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

double ConvexSphericalPolygon::compute_radius() const {
    double radius{0.};
    if (valid_) {
        PointXYZ centroid;
        centroid   = sph_coords_[0];
        size_t isp = 1;
        for (size_t i = 1; i < size_ - 1; ++i) {
            if (approx_eq(sph_coords_[isp], sph_coords_[isp - 1], TOL)) {
                continue;
            }
            centroid = centroid + sph_coords_[isp];
            ++isp;
        }
        centroid = PointXYZ::div(centroid, PointXYZ::norm(centroid));

        for (size_t i = 0; i < size_; ++i) {
            radius = std::max(radius, PointXYZ::distance(sph_coords_[i], centroid));
        }
    }
    return radius;
}

bool ConvexSphericalPolygon::GreatCircleSegment::contains(const PointXYZ& p) const {
    /*
   * @brief Point-on-segment test on great circle segments
   * @param[in] P given point in (x,y,z) coordinates
   * @return
   */
    constexpr double eps_large = 1.e3 * EPS;

    // Case where p is one of the endpoints
    double pp1n2 = distance2(p, p1_);
    double pp2n2 = distance2(p, p2_);
    if (pp1n2 < EPS2 or pp2n2 < EPS2) {
        return true;
    }

    PointXYZ p12 = cross();
    double p12n2 = norm2(p12);
    double p12n  = std::sqrt(p12n2);
    p12 /= p12n;

    if (std::abs(PointXYZ::dot(p, p12)) > eps_large) {
        return false;
    }
    double pp  = PointXYZ::distance(p1_, p2_);
    double pp1 = PointXYZ::distance(p, p1_);
    double pp2 = PointXYZ::distance(p, p2_);
    return (std::min(pp - pp1, pp - pp2) > -eps_large);
}

PointXYZ ConvexSphericalPolygon::GreatCircleSegment::intersect(const GreatCircleSegment& p) const {
    const auto& s = *this;
    PointXYZ sp   = PointXYZ::cross(s.cross(), p.cross());

    double sp_norm = PointXYZ::norm(sp);
    if (sp_norm > EPS) {
        sp /= sp_norm;
        if (p.contains(sp)) {
            return sp;
        }
        sp *= -1.;
        if (p.contains(sp)) {
            return sp;
        }
        return PointXYZ(0, 0, 0);
    }
    else {
        return PointXYZ(1, 1, 1);
    }
}


//------------------------------------------------------------------------------------------------------

}  // namespace util
}  // namespace atlas
