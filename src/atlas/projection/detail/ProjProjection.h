/*
 * (C) Copyright 1996- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */


#pragma once

#include <memory>

#include "atlas/projection/detail/ProjectionImpl.h"
#include "atlas/util/Config.h"

namespace atlas {
namespace projection {
namespace detail {

// Forward declaration of Proj, defined in ProjProjection.cc
struct Proj;

class ProjProjection final : public ProjectionImpl {
public:
    ProjProjection(const eckit::Parametrisation&);

    static std::string static_type() { return "proj"; }

    // -- Overridden methods

    std::string type() const override { return static_type(); }

    void xy2lonlat(double[]) const override;
    void lonlat2xy(double[]) const override;

    Jacobian jacobian(const PointLonLat&) const override;

    PointXYZ xyz(const PointLonLat&) const override;

    bool strictlyRegional() const override { return false; }
    RectangularLonLatDomain lonlatBoundingBox(const Domain&) const override;

    Spec spec() const override;
    std::string units() const override;
    void hash(eckit::Hash&) const override;

private:
    Normalise normalise_;
    std::string proj_string_;
    std::string source_;
    std::string geocentric_;
    bool source_encoded_;
    bool geocentric_encoded_;

    std::unique_ptr<Proj> proj_;

    Spec extraSpec_;
};


}  // namespace detail
}  // namespace projection
}  // namespace atlas
