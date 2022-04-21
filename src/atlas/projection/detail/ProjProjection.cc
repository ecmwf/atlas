/*
 * (C) Copyright 1996- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */


#include "ProjProjection.h"

#include <proj.h>

#include "eckit/types/FloatCompare.h"
#include "eckit/utils/Hash.h"

#include "atlas/projection/detail/ProjectionFactory.h"
#include "atlas/runtime/Exception.h"
#include "atlas/util/CoordinateEnums.h"
#include "atlas/util/Point.h"

namespace atlas {
namespace projection {
namespace detail {

namespace {
struct pj_t : std::unique_ptr<PJ, decltype(&proj_destroy)> {
    using t = std::unique_ptr<PJ, decltype(&proj_destroy)>;
    explicit pj_t(PJ* ptr): t(ptr, &proj_destroy) {}
    operator PJ*() const { return t::get(); }
};

struct ctx_t : std::unique_ptr<PJ_CONTEXT, decltype(&proj_context_destroy)> {
    using t = std::unique_ptr<PJ_CONTEXT, decltype(&proj_context_destroy)>;
    explicit ctx_t(PJ_CONTEXT* ptr): t(ptr, &proj_context_destroy) {}
    operator PJ_CONTEXT*() const { return t::get(); }
};

bool proj_ellipsoid_params(PJ_CONTEXT* ctxt, const std::string& proj_str, double& a, double& b) {
    bool success = false;
    pj_t identity(proj_create_crs_to_crs(ctxt, proj_str.c_str(), proj_str.c_str(), nullptr));
    pj_t proj(proj_get_target_crs(ctxt, identity));
    pj_t proj_ellps(proj_get_ellipsoid(ctxt, proj));
    if (proj_ellps) {
        ATLAS_ASSERT(proj_ellipsoid_get_parameters(ctxt, proj_ellps, &a, &b, nullptr, nullptr));
        ATLAS_ASSERT(0 < b && b <= a);
        success = true;
    }
    return success;
};

std::string source_str(PJ_CONTEXT* ctxt, const std::string& proj_str) {
    double a, b;
    if (proj_ellipsoid_params(ctxt, proj_str, a, b)) {
        auto axes = b < a ? " +a=" + std::to_string(a) + " +b=" + std::to_string(b) : " +R=" + std::to_string(a);
        return "+proj=lonlat" + axes;
    }
    else {
        return "EPSG:4326";  // WGS84 (lat,lon)
    }
}


std::string geocentric_str(PJ_CONTEXT* ctxt, const std::string& proj_str) {
    double a, b;
    if (proj_ellipsoid_params(ctxt, proj_str, a, b)) {
        return b < a ? "+proj=geocent +a=" + std::to_string(a) + " +b=" + std::to_string(b)
                     : "+proj=geocent +R=" + std::to_string(a);
    }
    else {
        return "EPSG:4978";  // WGS84 (x,y,z)
    }
}

}  // namespace


struct Proj {
    pj_t sourceToTarget_{nullptr};
    pj_t sourceToGeocentric_{nullptr};
    ctx_t context_{PJ_DEFAULT_CTX};
};

ProjProjection::ProjProjection(const eckit::Parametrisation& param): normalise_(param), proj_(new Proj()) {
    ATLAS_ASSERT(param.get("proj", proj_string_) && !proj_string_.empty());
    source_encoded_     = param.get("proj_source", source_ = source_str(proj_->context_, proj_string_));
    geocentric_encoded_ = param.get("proj_geocentric", geocentric_ = geocentric_str(proj_->context_, proj_string_));

    // set x/y transformations to/from lon/lat and to/from geocentric coordinates
    {
        pj_t p1(proj_create_crs_to_crs(proj_->context_, source_.c_str(), proj_string_.c_str(), nullptr));
        ATLAS_ASSERT(p1);
        proj_->sourceToTarget_.reset(proj_normalize_for_visualization(proj_->context_, p1));
        ATLAS_ASSERT(proj_->sourceToTarget_);
    }

    {
        pj_t p2(proj_create_crs_to_crs(proj_->context_, source_.c_str(), geocentric_.c_str(), nullptr));
        ATLAS_ASSERT(p2);
        proj_->sourceToGeocentric_.reset(proj_normalize_for_visualization(proj_->context_, p2));
        ATLAS_ASSERT(proj_->sourceToGeocentric_);
    }

    // set semi-major/minor axis
    {
        double a, b;
        if (proj_ellipsoid_params(proj_->context_, proj_string_, a, b)) {
            eckit::types::is_approximately_equal(a, b) ? extraSpec_.set("radius", a)
                                                       : extraSpec_.set("semi_major_axis", a).set("semi_minor_axis", b);
        }
    }

    // set units
    {
        pj_t target(proj_get_target_crs(proj_->context_, proj_->sourceToTarget_));
        ATLAS_ASSERT(target);

        pj_t coord(proj_crs_get_coordinate_system(proj_->context_, target));
        ATLAS_ASSERT(coord);
        ATLAS_ASSERT(proj_cs_get_axis_count(proj_->context_, coord) > 0);

        const char* units_c_str;
        if (proj_cs_get_axis_info(nullptr, coord, 0, nullptr, nullptr, nullptr, nullptr, &units_c_str, nullptr,
                                  nullptr)) {
            std::string units(units_c_str);
            if (!units.empty()) {
                if (units == "metre") {
                    units = "meters";
                }
                extraSpec_.set("units", units);
            }
        }
    }
}  // namespace detail


void ProjProjection::xy2lonlat(double crd[]) const {
    PJ_COORD P = proj_coord(crd[XX], crd[YY], 0, 0);
    P          = proj_trans(proj_->sourceToTarget_, PJ_INV, P);

    //    std::memcpy(crd, &P, 2 * sizeof(double));
    crd[LON] = P.enu.e;
    crd[LAT] = P.enu.n;

    normalise_(crd);
}


void ProjProjection::lonlat2xy(double crd[]) const {
    PJ_COORD P = proj_coord(crd[LON], crd[LAT], 0, 0);
    P          = proj_trans(proj_->sourceToTarget_, PJ_FWD, P);

    //    std::memcpy(crd, &P, 2 * sizeof(double));
    crd[XX] = P.xy.x;
    crd[YY] = P.xy.y;
}


ProjectionImpl::Jacobian ProjProjection::jacobian(const PointLonLat&) const {
    throw_NotImplemented("ProjProjection::jacobian", Here());
}

PointXYZ ProjProjection::xyz(const PointLonLat& lonlat) const {
    PJ_COORD P = proj_coord(lonlat.lon(), lonlat.lat(), 0, 0);
    P          = proj_trans(proj_->sourceToGeocentric_, PJ_FWD, P);
    return {P.xyz.x, P.xyz.y, P.xyz.z};
}


RectangularLonLatDomain ProjProjection::lonlatBoundingBox(const Domain& domain) const {
    return ProjectionImpl::lonlatBoundingBox(domain);
}


ProjProjection::Spec ProjProjection::spec() const {
    Spec spec;
    spec.set("type", type());
    spec.set("proj", proj_string_);
    if (source_encoded_) {
        spec.set("proj_source", source_);
    }
    if (geocentric_encoded_) {
        spec.set("proj_geocentric", geocentric_);
    }
    normalise_.spec(spec);
    return spec | extraSpec_;
}


std::string ProjProjection::units() const {
    return extraSpec_.getString("units", "");
}


void ProjProjection::hash(eckit::Hash& h) const {
    h.add(type());
    h.add(proj_string_);
    if (source_encoded_) {
        h.add(source_);
    }
    if (geocentric_encoded_) {
        h.add(geocentric_);
    }
    normalise_.hash(h);
}


namespace {
static ProjectionBuilder<ProjProjection> register_1(ProjProjection::static_type());
}


}  // namespace detail
}  // namespace projection
}  // namespace atlas
