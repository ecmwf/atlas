/**
* (C) Crown copyright 2021, Met Office
*
* This software is licensed under the terms of the Apache Licence Version 2.0
* which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
*/

#pragma once

#include <array>

#include "atlas/domain.h"
#include "atlas/projection/detail/ProjectionImpl.h"
#include "atlas/util/NormaliseLongitude.h"


namespace atlas {
namespace projection {
namespace detail {


template <typename Rotation>
class VariableResolutionProjectionT final : public ProjectionImpl {
public:
    ///<using Spec = ProjectionImpl::Spec;

    ///< constructor uses parametrisation and point to stretch
    VariableResolutionProjectionT(const eckit::Parametrisation& p);
    VariableResolutionProjectionT();

    ///< projection name
    static std::string static_type() { return Rotation::typePrefix() + "variable_resolution"; }
    std::string type() const override { return static_type(); }

    // projection and inverse projection

    void xy2lonlat(double crd[]) const override;
    void lonlat2xy(double crd[]) const override;

    ///< specification for stretching
    Spec spec() const override;

    std::string units() const override { return "degrees"; }

    void hash(eckit::Hash&) const override;

    ///< NOT for the moment
    Jacobian jacobian(const PointLonLat&) const override;

    bool strictlyRegional() const override { return true; }  ///< Stretch projection cannot be used for global grids
    RectangularLonLatDomain lonlatBoundingBox(const Domain& domain) const override {
        return ProjectionImpl::lonlatBoundingBox(domain);
    }

    void checkvalue(const double&, const double&) const;
    double general_stretch(const double, const bool, const int) const;
    double general_stretch_inv(const double, const bool, const int) const;

protected:
    double delta_outer;   ///< resolution of the external regular grid (rim) it should be larger than the last stretch
    double delta_inner;   ///< /< resolution of the regional model (regular grid)
    double var_ratio_;    ///< power used for the stretching
    double x_reg_start_;  ///< xstart of the internal regional grid
    double y_reg_start_;  ///< ystart of the internal regional grid
    double x_reg_end_;    ///< xend of the internal regional grid
    double y_reg_end_;    ///< yend of the internal regional grid
    double startx_;       ///< original domain startx
    double endx_;         ///< original domain endx
    double starty_;       ///< original domain starty
    double endy_;         ///< original domain endy
    double rim_widthx_;   ///< xsize of the rim
    double rim_widthy_;   ///< ysize of the rim


    //< variables derived from the configuration used for the projection
    double deltax_all, deltay_all;
    double add_xf_, add_yf_;
    int nx_stretched, ny_stretched, nx_rim, ny_rim;
    double check_x, check_y, check_st;
    double lam_hires_size;  ///< size regular grid x
    double phi_hires_size;  ///< size regular grid y
    double lambda_start;    ///< start grid x
    double phi_start;       ///< start grid y
    std::array<double, 2> new_ratio_;

    void setup(const eckit::Parametrisation& p);

private:
    Rotation rotation_;
};

using VariableResolutionProjection        = VariableResolutionProjectionT<NotRotated>;
using RotatedVariableResolutionProjection = VariableResolutionProjectionT<Rotated>;

}  // namespace detail
}  // namespace projection
}  // namespace atlas
