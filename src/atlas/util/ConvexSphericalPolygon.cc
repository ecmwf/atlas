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
#include "atlas/util/ConvexSphericalPolygon.h"
#include "atlas/util/CoordinateEnums.h"
#include "atlas/util/NormaliseLongitude.h"

#define DEBUG_OUTPUT 1
#define DEBUG_OUTPUT_DETAIL 1

namespace atlas {
namespace util {

bool approx_eq(const double& v1, const double& v2, const double& tol) {
    return eckit::types::is_approximately_equal(v1, v2, tol);
}

bool approx_eq(const PointXYZ& v1, const PointXYZ& v2, const double& tol) {
    //return approx_eq( v1[0], v2[0], t ) && approx_eq( v1[1], v2[1], t ) && approx_eq( v1[2], v2[2], t );
    return PointXYZ::norm(v1 - v2) < tol;
}

bool approx_eq_null(const PointXYZ& v1, const double& tol) {
    //return approx_eq( v1[0], 0., t ) && approx_eq( v1[1], 0., t ) && approx_eq( v1[2], 0., t );
    return PointXYZ::norm(v1) < tol;
}

double cart_diff(const PointXYZ& v1, const PointXYZ& v2) {
    return PointXYZ::norm(v1 - v2);
}

PointLonLat sph_to_lonlat(const PointXYZ& p) {
    PointLonLat pp;
    eckit::geometry::Sphere::convertCartesianToSpherical(1., p, pp);
    return pp;
}

//------------------------------------------------------------------------------------------------------

ConvexSphericalPolygon::ConvexSphericalPolygon(): valid_(false), size_(0), area_(0), cell_radius_(0) {}

ConvexSphericalPolygon::ConvexSphericalPolygon(const std::vector<PointLonLat>& points): size_(points.size()), area_(0) {
    ATLAS_ASSERT(size_ < MAX_SIZE);
    eckit::geometry::Sphere::convertSphericalToCartesian(1., points[0], sph_coords_[0]);
    centroid_  = sph_coords_[0];
    size_t isp = 1;
    for (size_t i = 1; i < points.size() - 1; ++i) {
        eckit::geometry::Sphere::convertSphericalToCartesian(1., points[i], sph_coords_[isp]);
        if (approx_eq(sph_coords_[isp], sph_coords_[isp - 1], 1e-10)) {
            continue;
        }
        centroid_ = centroid_ + sph_coords_[isp];
        ++isp;
    }
    eckit::geometry::Sphere::convertSphericalToCartesian(1., points[points.size() - 1], sph_coords_[isp]);
    if (approx_eq(sph_coords_[isp], sph_coords_[0], 1e-10) or
        approx_eq(sph_coords_[isp], sph_coords_[isp - 1], 1e-10)) {
    }
    else {
        centroid_ = centroid_ + sph_coords_[isp];
        ++isp;
    }
    size_  = isp;
    valid_ = size_ > 2;
    if (valid_) {
        ATLAS_ASSERT(validate());
        ATLAS_ASSERT(not approx_eq_null(centroid_, 1e-10));
        centroid_    = PointXYZ::div(centroid_, PointXYZ::norm(centroid_));
        cell_radius_ = 0.;
        for (size_t i = 0; i < size_; ++i) {
            cell_radius_ = std::max(cell_radius_, cart_diff(sph_coords_[i], centroid_));
        }
        compute_area();
    }
}

ConvexSphericalPolygon::ConvexSphericalPolygon(const std::vector<PointXYZ>& points, const bool debug):
    size_(points.size()), area_(0), cell_radius_(0) {
#if DEBUG_OUTPUT
    if (debug) {
        std::cout << " \n\nConvexSphericalPolygon got points: ";
        for (int i = 0; i < points.size(); ++i) {
            std::cout << std::setprecision(10) << sph_to_lonlat(points[i]) << ",";
        }
        (std::cout << "\n").flush();
    }
#endif
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
    for (; i < points.size() && search_3pts; ++i) {
        const PointXYZ& P0 = points[i];
        for (j = i + 1; j < points.size() && search_3pts; ++j) {
            const PointXYZ& P1 = points[j];
            if (approx_eq(P0, P1, 1e-10)) {
                continue;
            }
            for (k = j + 1; k < points.size() && search_3pts; ++k) {
                const PointXYZ& P2 = points[k];
                if (approx_eq(P1, P2, 1e-10) or approx_eq(P0, P2, 1e-10)) {
                    continue;
                }
                if (leftOf(P2, P0, P1, 1e-14, debug)) {
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
#if DEBUG_OUTPUT
    if (debug) {
        std::cout << " ConvexSphericalPolygon 3 FIRST: " << std::setprecision(10) << sph_to_lonlat(sph_coords_[0])
                  << " " << sph_to_lonlat(sph_coords_[1]) << " " << sph_to_lonlat(sph_coords_[2]) << " "
                  << "\n";
    }
#endif

    centroid_ = sph_coords_[i] + sph_coords_[j] + sph_coords_[k];
    for (; k < points.size() - 1; ++k) {
        if (approx_eq(points[k], sph_coords_[isp - 1], 1e-10) or
            (not leftOf(points[k], sph_coords_[isp - 2], sph_coords_[isp - 1], 1e-14, debug))) {
            continue;
        }
        sph_coords_[isp] = points[k];
        centroid_        = centroid_ + sph_coords_[isp];
        isp++;
    }
    const PointXYZ& Pl2 = sph_coords_[isp - 2];
    const PointXYZ& Pl1 = sph_coords_[isp - 1];
    const PointXYZ& P0  = sph_coords_[0];
    const PointXYZ& P   = points[size_ - 1];
    if ((not approx_eq(P, P0, 1e-14)) and (not approx_eq(P, Pl1, 1e-14)) and leftOf(P, Pl2, Pl1, 1e-14, debug) and
        leftOf(P0, Pl1, P, 1e-14, debug)) {
        sph_coords_[isp] = P;
        centroid_        = centroid_ + P;
        ++isp;
    }
    size_  = isp;
    valid_ = size_ > 2;
    if (valid_) {
        ATLAS_ASSERT(not approx_eq_null(centroid_, 1e-10));
        centroid_    = PointXYZ::div(centroid_, PointXYZ::norm(centroid_));
        cell_radius_ = 0.;
        for (size_t i = 0; i < size_; ++i) {
            cell_radius_ = std::max(cell_radius_, cart_diff(sph_coords_[i], centroid_));
        }
        compute_area();
    }
#if DEBUG_OUTPUT
    if (debug) {
        std::cout << " Final polygon. Centroid : " << centroid_ << ", radius : " << cell_radius_ << ", area : " << area_
                  << "\n";
    }
#endif
}

bool ConvexSphericalPolygon::validate() {
    if (valid_) {
        for (int i = 0; i < size(); i++) {
            int ni                = (i != size() - 1 ? i + 1 : 0);
            int nni               = (ni != size() - 1 ? ni + 1 : 0);
            const PointXYZ& P     = sph_coords_[i];
            const PointXYZ& nextP = sph_coords_[ni];
            ATLAS_ASSERT(std::abs(PointXYZ::dot(P, P) - 1.) < 1e-14);
            ATLAS_ASSERT(not approx_eq(P, PointXYZ::mul(nextP, -1.), 1e-10));
            valid_ = valid_ && leftOf(sph_coords_[nni], P, nextP, 1e-14, 0);
        }
    }
    return valid_;
}

bool ConvexSphericalPolygon::equals(const ConvexSphericalPolygon& plg, const double deg_prec) const {
    const double le = 2. * std::sin(M_PI * deg_prec / 360.);
    if ((not plg.valid_) || (not valid_) || size_ != plg.size()) {
        (Log::info() << " ConvexSphericalPolygon::equals == not compatible\n").flush();
        return false;
    }
    int offset = 0;
    for (; offset < size_; ++offset) {
        if (PointXYZ::norm(plg.sph_coords_[0] - sph_coords_[offset]) < le) {
            break;
        }
    }
    if (offset == size_) {
        (Log::info() << "ConvexSphericalPolygon::equals == no point equal"
                     << "\n")
            .flush();
        return false;
    }

    for (int j = 0; j < size_; j++) {
        int idx   = (offset + j) % size_;
        auto dist = PointXYZ::norm(plg.sph_coords_[j] - sph_coords_[idx]);
        if (dist > le) {
            (Log::info() << "  ConvexSphericalPolygon::equals == point distance " << dist << "\n").flush();
            return false;
        }
    }
    return true;
}

// note: unit sphere!
// I. Todhunter (1886), Paragr. 99
void ConvexSphericalPolygon::compute_area(const int debug) {
    area_ = 0.;
    if (size_ < 3) {
        return;
    }
    std::vector<double> tri_area;
    std::vector<PointXYZ> tri_centre;
    if (cell_radius_ < 1e-6) {  // plane area
        for (int i = 1; i < size_ - 1; i++) {
            const PointXYZ& pl = sph_coords_[i] - sph_coords_[0];
            const PointXYZ& pr = sph_coords_[i + 1] - sph_coords_[0];
            auto ctr           = sph_coords_[0] + sph_coords_[i] + sph_coords_[i + 1];
            const double tarea = 0.5 * PointXYZ::norm(PointXYZ::cross(pl, pr));
            tri_centre.emplace_back(PointXYZ::div(ctr, PointXYZ::norm(ctr)));
            tri_area.emplace_back(tarea);
            area_ += tarea;
        }
    }
    else {  // spherical area
        const PointXYZ& a = sph_coords_[0];
        for (int i = 1; i < size_ - 1; i++) {
            const PointXYZ& b    = sph_coords_[i];
            const PointXYZ& c    = sph_coords_[i + 1];
            PointXYZ ab          = PointXYZ(PointXYZ::cross(a, b));
            PointXYZ bc          = PointXYZ(PointXYZ::cross(b, c));
            PointXYZ ca          = PointXYZ(PointXYZ::cross(c, a));
            const double ab_norm = PointXYZ::norm(ab);
            const double bc_norm = PointXYZ::norm(bc);
            const double ca_norm = PointXYZ::norm(ca);
            if (ab_norm < 1e-15 or bc_norm < 1e-15 or ca_norm < 1e-15) {
                continue;
            }
            double abc = -PointXYZ::dot(ab, bc) / (ab_norm * bc_norm);
            double bca = -PointXYZ::dot(bc, ca) / (bc_norm * ca_norm);
            double cab = -PointXYZ::dot(ca, ab) / (ca_norm * ab_norm);
#if DEBUG_OUTPUT_DETAIL
            if (debug) {
                (Log::info() << " === compute_area:: a, b, c: " << a << " " << b << " " << c << "\n").flush();
                (Log::info() << " === compute_area:: ab, bc, ca: " << ab << " " << bc << " " << ca << "\n").flush();
                (Log::info() << " === compute_area:: abc, bca, cab: " << std::abs(abc) - 1. << " " << std::abs(bca) - 1.
                             << " " << std::abs(cab) - 1. << "\n")
                    .flush();
                (Log::info() << " === compute_area::norms: " << ab_norm << " " << bc_norm << " " << ca_norm << "\n")
                    .flush();
            }
#endif
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
#if DEBUG_OUTPUT
            if (debug) {
                Log::info() << "\t\tabc, bca, cab, +area : " << abc << ", " << bca << ", " << cab << ", "
                            << abc + bca + cab - M_PI << "\n";
            }
#endif
            auto ctr           = a + b + c;
            const double tarea = abc + bca + cab - M_PI;
            tri_centre.emplace_back(PointXYZ::div(ctr, PointXYZ::norm(ctr)));
            tri_area.emplace_back(tarea);
            area_ += tarea;
        }
    }
    // fix centroid to be based on area rather than on vertices of polygon
    centroid_ = PointXYZ{0., 0., 0.};
    if (area_ > 0.) {
        for (int i = 1; i < size_ - 1; i++) {
            centroid_ = centroid_ + PointXYZ::mul(tri_centre[i - 1], tri_area[i - 1]);
        }
        centroid_ = PointXYZ::div(centroid_, PointXYZ::norm(centroid_));
    }
}

int ConvexSphericalPolygon::leftOf(const PointXYZ& P, const PointXYZ& p1, const PointXYZ& p2, const double tol,
                                   const int debug) {
    const PointXYZ& cp = PointXYZ(PointXYZ::cross(p1, p2));
    double cpP         = PointXYZ::dot(cp, P);
#if DEBUG_OUTPUT
    if (debug) {
        (Log::info() << "\tp, p1, p2: " << sph_to_lonlat(P) << ", " << sph_to_lonlat(p1) << ", " << sph_to_lonlat(p2)
                     << "\n")
            .flush();
        (Log::info() << "\tp1 x p2, p * (p1 x p2): " << cp << ", " << cpP << "\n").flush();
    }
#endif
    return (cpP > -tol);
}

double ConvexSphericalPolygon::norm_max(const PointXYZ& p, const PointXYZ& q) {
    double n01 = std::max(std::abs(p[0] - q[0]), std::abs(p[1] - q[1]));
    return std::max(n01, std::abs(p[2] - q[2]));
}

bool ConvexSphericalPolygon::between(const PointXYZ& p, const PointXYZ& p1, const PointXYZ& p2, const int debug) {
    PointXYZ p12 = PointXYZ::cross(p1, p2);
    double p12n  = PointXYZ::norm(p12) - std::numeric_limits<double>::epsilon();
    double pp1n  = PointXYZ::norm(p - p1) - std::numeric_limits<double>::epsilon();
    double pp2n  = PointXYZ::norm(p - p2) - std::numeric_limits<double>::epsilon();
    if (p12n < 0. && pp1n < 0.) {
        return true;
    }
    p12             = PointXYZ::div(p12, p12n);
    const double dp = 1e+3 * std::numeric_limits<double>::epsilon() - std::abs(PointXYZ::dot(p, p12));
    if (dp < 0.) {
#if DEBUG_OUTPUT_DETAIL
        if (debug) {
            (Log::info() << "NOT between: point not in the plane, " << dp << "\n").flush();
        }
#endif
        return false;
    }
    double pp = PointXYZ::norm(p1 - p2);
    pp        = std::min(pp - pp1n, pp - pp2n) + 5e-15;
#if DEBUG_OUTPUT_DETAIL
    if (debug) {
        (Log::info() << "  between pp, pp2n = " << pp << ", " << pp2n << "\n").flush();
    }
#endif
    return (pp > 0. && pp2n > 0.);
}

PointXYZ ConvexSphericalPolygon::common(const PointXYZ& s1, const PointXYZ& s2, const PointXYZ& p1, const PointXYZ& p2,
                                        const int debug) {
    PointXYZ s  = static_cast<PointXYZ>(PointXYZ::cross(s1, s2));
    PointXYZ p  = static_cast<PointXYZ>(PointXYZ::cross(p1, p2));
    PointXYZ sp = static_cast<PointXYZ>(PointXYZ::cross(s, p));

#if 0
	Log::info() << "s1-s2 = " << s1 - s2 << ", |s1-s2|: " << PointXYZ::norm( s1 - s2 ) << "\n";
	Log::info() << "p1-p2 = " << p1 - p2 << ", |p1-p2|: " << PointXYZ::norm( p1 - p2 ) << "\n";
	Log::info() << " s = " << s << ", |s| = " << PointXYZ::norm( s ) << "\n";
	Log::info() << " s = " << sph_to_lonlat(s) << "\n";
	Log::info() << " p = " << p << ", |p| = " << PointXYZ::norm( p ) << "\n";
	Log::info() << " p = " << sph_to_lonlat(p) << "\n";
	Log::info() << " sp = " << sp << ", |sp| = " << PointXYZ::norm( sp ) << "\n";
	Log::info() << " sp = " << sph_to_lonlat(sp) << "\n\n";
#endif

    double sp_norm = PointXYZ::norm(sp) - std::numeric_limits<double>::epsilon();
#if DEBUG_OUTPUT_DETAIL
    if (debug) {
        (Log::info() << " Parallel: " << sp_norm << " < 0 ?\n").flush();
    }
#endif
    if (sp_norm > 0.) {
        sp = PointXYZ::div(sp, PointXYZ::norm(sp));
#if DEBUG_OUTPUT
        if (debug) {
            (Log::info() << "Test sp = " << sph_to_lonlat(sp) << "\n").flush();
            (Log::info() << "	p1-sp-p2: " << between(sp, p1, p2) << "\n").flush();
        }
#endif
        if (between(sp, p1, p2, debug)) {
#if DEBUG_OUTPUT_DETAIL
            if (debug) {
                (Log::info() << ",	GOT  p: " << sph_to_lonlat(sp) << "\n").flush();
            }
#endif
            return sp;
        }
        sp = PointXYZ::mul(sp, -1);
#if DEBUG_OUTPUT_DETAIL
        if (debug) {
            (Log::info() << "Test sp = " << sph_to_lonlat(sp) << "\n").flush();
        }
#endif
        if (between(sp, p1, p2, debug)) {
#if DEBUG_OUTPUT_DETAIL
            if (debug) {
                (Log::info() << ",	GOT -p: " << sph_to_lonlat(sp) << "\n").flush();
            }
#endif
            return sp;
        }
#if DEBUG_OUTPUT_DETAIL
        if (debug) {
            (Log::info() << " Intersection not on [p1, p2).\n").flush();
        }
#endif
        return PointXYZ({0, 0, 0});
    }
    else {
#if DEBUG_OUTPUT_DETAIL
        if (debug) {
            (Log::info() << " Overlap. \n").flush();
        }
#endif
        return PointXYZ({1, 1, 1});
    }
}

// @return -1: overlap with one of polygon edges,
//  		0: no_intersect,
//          1 + (id of this polygon's segment intersecting [s1,s2)): otherwise
int ConvexSphericalPolygon::intersect(const PointXYZ& s1, const PointXYZ& s2, PointXYZ& I, int start,
                                      const bool debug) const {
    for (int i = start; i < size_; i++) {
        const int id0      = i;
        const int id1      = (id0 == size_ - 1) ? 0 : id0 + 1;
        const PointXYZ& p1 = sph_coords_[id0];
        const PointXYZ& p2 = sph_coords_[id1];
#if DEBUG_OUTPUT
        if (debug) {
            (Log::info() << "	** edge: " << sph_to_lonlat(p1) << " " << sph_to_lonlat(p2) << "\n").flush();
        }
#endif
        I = common(s1, s2, p1, p2, debug);
        if (I[0] == 0 && I[1] == 0 && I[2] == 0) {
            // intersection not on [p1,p2)
            continue;
        }
        if (I[0] == 1 && I[1] == 1) {
            // overlap
            return -1;
        }
        return 1 + id0;
    }
    return 0;
}

void ConvexSphericalPolygon::clip(const PointLonLat& s1, const PointLonLat& s2, const int debug) {
    PointXYZ p1;
    PointXYZ p2;
    eckit::geometry::Sphere::convertSphericalToCartesian(1., s1, p1);
    eckit::geometry::Sphere::convertSphericalToCartesian(1., s2, p2);
    clip(p1, p2, debug);
}

void ConvexSphericalPolygon::clip(const PointXYZ& s1, const PointXYZ& s2, const int debug) {
    if (PointXYZ::norm(s1 - s2) < 1e-14) {
        return;
    }
    PointXYZ i1    = PointXYZ({0, 0, 0});
    PointXYZ i2    = PointXYZ({0, 0, 0});
    int f1         = -1;
    int f2         = -1;
    f1             = -1 + intersect(s1, s2, i1, 0, debug);
    int f1n        = (f1 != size_ - 1) ? f1 + 1 : 0;
    bool i1_is_f1n = (PointXYZ::norm(i1 - sph_coords_[f1n]) < 1e-8);
    if (f1 >= 0) {
        f2 = -1 + intersect(s1, s2, i2, f1 + 1 + (i1_is_f1n ? 1 : 0), debug);
    }
    std::vector<bool> old_pt_in(size_);
    for (int i = 0; i < size_; i++) {
        old_pt_in[i] = false;
#if DEBUG_OUTPUT
        if (debug) {
            (Log::info() << " point " << sph_to_lonlat(sph_coords_[i]) << "\n").flush();
        }
#endif
        if (leftOf(sph_coords_[i], s1, s2, 1e-15, debug)) {
            old_pt_in[i] = true;
        }
    }
    bool add_i1 = false;
    bool add_i2 = false;
    if (f1 >= 0) {
        add_i1 = true;
        if (PointXYZ::norm(i1 - sph_coords_[f1]) < 1e-8) {
            add_i1 = add_i1 && (old_pt_in[f1] ? false : true);
        }
        if (i1_is_f1n) {
            add_i1 = add_i1 && (old_pt_in[f1n] ? false : true);
        }
        add_i1 = add_i1 && leftOf(i1, sph_coords_[f1], sph_coords_[f1n], 1e-10, debug);
    }
    else {
        f1 = -1;
    }
    if (f2 >= 0) {
        add_i2  = true;
        int f2n = (f2 != size_ - 1) ? f2 + 1 : 0;
        if (PointXYZ::norm(i2 - sph_coords_[f2]) < 1e-8) {
            add_i2 = add_i2 && (old_pt_in[f2] ? false : true);
        }
        if (PointXYZ::norm(i2 - sph_coords_[f2n]) < 1e-8) {
            add_i2 = add_i2 && (old_pt_in[f2n] ? false : true);
        }
        add_i2 = add_i2 && leftOf(i2, sph_coords_[f2], sph_coords_[f2n], 1e-10, debug);
    }
    else {
        f2 = -1;
    }
#if DEBUG_OUTPUT
    if (debug) {
        (Log::info() << "f1, f2, i1, i2: " << f1 << " " << f2 << " " << sph_to_lonlat(i1) << ", " << sph_to_lonlat(i2)
                     << "\n")
            .flush();
        (Log::info() << "add_i1, add_i2: " << add_i1 << " " << add_i2 << "\n").flush();
    }
#endif
    int f2_end = f2;
    if (f1 >= 0 && f2 == -1) {
        f2 = f1 + 1;
    }
    std::vector<PointXYZ> cp_sph_coords;
    cp_sph_coords.reserve(size_ + 2);
    for (int i = 0; i <= f1; i++) {
        if (old_pt_in[i]) {
            cp_sph_coords.emplace_back(sph_coords_[i]);
        }
    }
    if (add_i1) {
        cp_sph_coords.emplace_back(i1);
    }
    for (int i = f1 + 1; i <= f2; i++) {
        if (old_pt_in[i]) {
            cp_sph_coords.emplace_back(sph_coords_[i]);
        }
    }
    if (add_i2) {
        cp_sph_coords.emplace_back(i2);
    }
    for (int i = f2 + 1; i < size_; i++) {
        if (old_pt_in[i]) {
            cp_sph_coords.emplace_back(sph_coords_[i]);
        }
    }

    size_ = cp_sph_coords.size();
    for (int i = 0; i < size_; i++) {
        sph_coords_[i] = cp_sph_coords[i];
    }

    if (size_ < 3) {
        valid_ = false;
        area_  = 0.;
    }

#if DEBUG_OUTPUT
    if (debug) {
        (Log::info() << " new size " << size_ << "\n").flush();
        (Log::info() << " New plg: ").flush();
        for (int i = 0; i < size_; i++) {
            Log::info() << "," << sph_to_lonlat(sph_coords_[i]);
        }
        (Log::info() << "\n\n").flush();
    }
#endif
}

// intersect a polygon with this polygon
// @param[in] pol clipping polygon
// @param[out] intersecting polygon
ConvexSphericalPolygon ConvexSphericalPolygon::intersect(const ConvexSphericalPolygon& plg, const int debug) const {
    ConvexSphericalPolygon obj = *this;
    for (int i = 0; i < plg.size_ && bool(obj); i++) {
        const PointXYZ& s1 = plg.sph_coords_[i];
        const PointXYZ& s2 = plg.sph_coords_[(i != plg.size_ - 1) ? i + 1 : 0];
#if DEBUG_OUTPUT
        if (debug) {
            (Log::info() << std::setprecision(8) << "\n 	plg " << i << ": " << obj << "\n").flush();
            (Log::info() << " 	now clip with " << sph_to_lonlat(s1) << " " << sph_to_lonlat(s2) << "\n").flush();
        }
#endif
        obj.clip(s1, s2, debug);
    }
    std::vector<PointXYZ> pts;
    for (int i = 0; i < obj.size_; i++) {
        pts.emplace_back(obj.sph_coords_[i]);
    }
    return ConvexSphericalPolygon(pts, debug);
}

void ConvexSphericalPolygon::print(std::ostream& out) const {
    out << "{";
    for (size_t i = 0; i < size(); ++i) {
        if (i > 0) {
            out << ",";
        }
        PointLonLat ip_ll;
        eckit::geometry::Sphere::convertCartesianToSpherical(1., sph_coords_[i], ip_ll);
        out << ip_ll;
    }
    out << "}";
}


//------------------------------------------------------------------------------------------------------

}  // namespace util
}  // namespace atlas
