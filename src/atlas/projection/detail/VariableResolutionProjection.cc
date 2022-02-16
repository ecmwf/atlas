/**
* (C) Crown copyright 2021, Met Office
*
* This software is licensed under the terms of the Apache Licence Version 2.0
* which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
*/

#include <cmath>
#include <functional>
#include <sstream>

#include "eckit/config/Parametrisation.h"
#include "eckit/utils/Hash.h"

#include "atlas/functionspace/StructuredColumns.h"
#include "atlas/grid.h"
#include "atlas/grid/Grid.h"
#include "atlas/projection/detail/ProjectionFactory.h"
#include "atlas/projection/detail/VariableResolutionProjection.h"
#include "atlas/runtime/Exception.h"
#include "atlas/util/Config.h"
#include "atlas/util/Constants.h"
#include "eckit/testing/Test.h"

/**
* Projection for LAM stretching
*
* The theory used in this code can be found under:
* Lateral boundary conditions for limited area models
* Terry Davies 2014
*  https://doi.org/10.1002/qj.2127
*
* Some equations are under:
* The benefits of the Met Office variable resolution NWP model for forecasting convection
* Tang et al. 2013
* https://doi.org/10.1002/met.1300
*
* Generally: From a original regular grid the point in the middle of the grid are unchanged,
* based on the boundary of the new internal regular grid.
* The grid length is then stretched at a constant stretching (or inflation) factor r,
* that is Δsi+1 = rΔsi until all the points are mapped in the new grid.
* The domain is then extended to a further uniform coarse resolution region (rim) in the outer domain,
*
* The origin of the xy-system is at (lon,lat) = (0,0)
*
*/

#ifdef __NVCOMPILER
#define PREVENT_OPT volatile
#else
#define PREVENT_OPT
#endif

namespace atlas {
namespace projection {
namespace detail {

static double new_ratio(int n_stretched, double var_ratio) {
    /**
     *  compute ratio,
     *  change stretching factor so that high and low grids
     *  retain original sizes
     */
    ///< correction used to change from double to integer
    constexpr float epstest = std::numeric_limits<float>::epsilon();

    ///< number of variable (stretched) grid points in one side
    PREVENT_OPT int var_ints = (n_stretched + epstest) / 2.;
    double var_ints_f        = n_stretched / 2.;
    double logr              = std::log(var_ratio);
    double log_ratio         = (var_ints_f - 0.5) * logr;

    return std::exp(log_ratio / var_ints);
};


///< specification parameters
template <typename Rotation>
typename VariableResolutionProjectionT<Rotation>::Spec VariableResolutionProjectionT<Rotation>::spec() const {
    Spec proj_st;
    proj_st.set("type", static_type());
    proj_st.set("outer.dx", delta_outer );     ///< resolution of the external regular grid (rim)
    proj_st.set("inner.dx", delta_inner );      ///< resolution of the regional model (regular grid)
    proj_st.set("progression", var_ratio_ );      ///< power used for the stretching
    proj_st.set("inner.xmin", x_reg_start_);  ///< xstart of the internal regional grid
    proj_st.set("inner.ymin", y_reg_start_);  ///< ystart of the internal regional grid
    proj_st.set("inner.xend", x_reg_end_);      ///< xend of the regular part of stretched internal grid
    proj_st.set("inner.yend", y_reg_end_);      ///< yend of the regular part of stretched internal grid
    proj_st.set("outer.xmin", startx_);            ///< original domain startx
    proj_st.set("outer.xend", endx_);                ///< original domain endx
    proj_st.set("outer.ymin", starty_);            ///< original domain starty
    proj_st.set("outer.yend", endy_);                ///< original domain endy
    proj_st.set("rim_widthx", rim_widthx_);         ///< xsize of the rim
    proj_st.set("rim_widthy", rim_widthy_);         ///< ysize of the rim
    rotation_.spec(proj_st);
    return proj_st;
}


///< constructors
template <typename Rotation>
VariableResolutionProjectionT<Rotation>::VariableResolutionProjectionT(const eckit::Parametrisation& proj_st):
    ProjectionImpl(), rotation_(proj_st) {
    proj_st.get("outer.dx", delta_outer = 0.);     ///< resolution of the external regular grid (rim)
    proj_st.get("inner.dx", delta_inner = 0.);      ///< resolution of the regional model (regular grid)
    proj_st.get("progression", var_ratio_ = 0.);      ///< power used for the stretching
    proj_st.get("inner.xmin", x_reg_start_ = 0.);  ///< xstart of the internal regional grid
    proj_st.get("inner.ymin", y_reg_start_ = 0.);  ///< ystart of the internal regional grid
    proj_st.get("inner.xend", x_reg_end_ = 0.);      ///< xend of the regular part of stretched internal grid
    proj_st.get("inner.yend", y_reg_end_ = 0.);      ///< yend of the regular part of stretched internal grid
    proj_st.get("outer.xmin", startx_ = 0.);            ///< original domain startx
    proj_st.get("outer.xend", endx_ = 0.);                ///< original domain endx
    proj_st.get("outer.ymin", starty_ = 0.);            ///< original domain starty
    proj_st.get("outer.yend", endy_ = 0.);                ///< original domain endy
    proj_st.get("rim_widthx", rim_widthx_);         ///< xsize of the rim
    proj_st.get("rim_widthy", rim_widthy_);         ///< ysize of the rim

    constexpr float epsilon = std::numeric_limits<float>::epsilon();  ///< value used to check if the values are equal
    constexpr float epstest =
        std::numeric_limits<float>::epsilon();  ///< correction used to change from double to intege

    ///< original domain size includes the points for the rim
    deltax_all = (endx_ - startx_);
    deltay_all = (endy_ - starty_);


    nx_stretched = 0;
    ny_stretched = 0;
    nx_rim       = 0;
    ny_rim       = 0;

    if (var_ratio_ == 1) {
        lam_hires_size = deltax_all;
        phi_hires_size = deltay_all;
        lambda_start   = x_reg_start_;
        phi_start      = y_reg_start_;
    }
    else {
        lam_hires_size = x_reg_end_ - x_reg_start_;
        phi_hires_size = y_reg_end_ - y_reg_start_;

        /**
         *  check for start of the new regular LAM grid x and y
         *  in the middle of the previous regular grid
         */

        ///< distance end of the grid and internal regular grid
        add_xf_ = (deltax_all + epstest - lam_hires_size) / 2.;
        add_yf_ = (deltay_all + epstest - phi_hires_size) / 2.;
        /**
         *  Compute the number of points for different part of the grid
         *  internal regular grid high resolution
         *  stretched grid
         *  external regular grid low resolution
         *  +1 otherwise I have the intervals not the points
         */
        nx_rim       = rim_widthx_ / delta_outer;
        ny_rim       = rim_widthy_ / delta_outer;
        nx_stretched = ((deltax_all + epstest - lam_hires_size) / delta_inner) - nx_rim;
        ny_stretched = ((deltay_all + epstest - phi_hires_size) / delta_inner) - ny_rim;
        lambda_start = x_reg_start_;
        phi_start    = y_reg_start_;
        /**
         *  check if stretched grid is in the middle
         *  of the previous regular grid
         */
        check_x  = startx_ + add_xf_ - lambda_start;
        check_y  = starty_ + add_yf_ - phi_start;
        check_st = nx_stretched - ny_stretched;
        checkvalue(epsilon, check_x);
        checkvalue(epsilon, check_y);
        checkvalue(epsilon, check_st);
    }

    int nx_ = (deltax_all + epstest) / delta_inner + 1;
    int ny_ = (deltay_all + epstest) / delta_inner + 1;

    int nx_inner = (lam_hires_size + epstest) / delta_inner + 1;
    int ny_inner = (phi_hires_size + epstest) / delta_inner + 1;

    ATLAS_ASSERT((nx_ - 1) - nx_rim - (nx_inner - 1) == nx_stretched);
    ATLAS_ASSERT((ny_ - 1) - ny_rim - (ny_inner - 1) == ny_stretched);

    new_ratio_[0] = var_ratio_;
    new_ratio_[1] = var_ratio_;
    if (var_ratio_ != 1) {
        new_ratio_[0] = new_ratio(nx_stretched, var_ratio_);
        new_ratio_[1] = new_ratio(ny_stretched, var_ratio_);
    }
}


template <typename Rotation>
void VariableResolutionProjectionT<Rotation>::checkvalue(const double& epsilon, const double& value_check) const {
    //err_message = "USER defined limits not in the middle of the area " + str;
    if (value_check > epsilon || value_check < (-1. * epsilon)) {
        std::string err_message;
        std::string str = std::to_string(value_check);
        throw eckit::BadValue("USER defined limits not in the middle of the area " + str, Here());
    }
}

/**
 * General stretch from a point in regular grid to the
 * a correspective point in a variable grid
 */

template <typename Rotation>
double VariableResolutionProjectionT<Rotation>::general_stretch(const double crd, const bool L_long,
                                                                const int n_stretched) const {
    constexpr float epstest =
        std::numeric_limits<float>::epsilon();  ///< correction used to change from double to integer
    constexpr double epsrem =
        0.1 * std::numeric_limits<double>::epsilon() /
        std::numeric_limits<float>::epsilon();  ///< correction used to part the find a part of an integer

    double inner_size;   ///< number of new internal regular grid in double
    double inner_start;  ///< start of the regular grid
    double inner_end;    ///< end of the regular grid
    double point = crd;  ///< starting point

    auto normalised = [L_long](double p) {
        if (L_long) {
            p = (p < 180) ? p + 360.0 : p;
        }
        return p;
    };

    /*
     * regular grids
     */
    if (var_ratio_ == 1) {
        return normalised(point);
    }

    /*
     *  SECTION 1
     *  INTERNAL REGULAR GRID the point is mapped to the same point
     */

    if (L_long) {
        inner_start = x_reg_start_;
        inner_size  = lam_hires_size;
    }
    else {
        inner_start = y_reg_start_;
        inner_size  = phi_hires_size;
    }
    inner_end = inner_start + inner_size;

    if ((point >= inner_start) && (point <= inner_end)) {
        return normalised(point);
    }

    /* SECTION 2
     * Start variable res 'STRETCHED'
     * The point is mapped in a stretched grid related to
     * distance from the regular grid and stretch internal grid: delta_dist
     */

    double distance_to_inner = 0.;  ///< distance from point to reg. grid
    double delta_add = 0.;          ///< additional part in stretch different from internal high resolution
    int n_high = 0;                ///< number of points, from point to reg. grid
    int n_high_st = 0;             ///< number of stretched points, from point to reg grid
    int n_high_rim = 0;            ///< number of rim points, from point to reg grid
    double p_rem = 0.;              ///< remaining part in stretch if distance_to_inner not multiple of delta_high_
    double p_rem_low = 0.;          ///< remaining part in rim if distance_to_inner not multiple of delta_high_
    double new_ratio = new_ratio_[L_long ? 1 : 0];

    if (point < inner_start) {
        distance_to_inner = inner_start - point;
    }
    else if (point > inner_end) {
        distance_to_inner = point - inner_end;
    }

    /*
     * number of high resolution points intervals, that are not
     * in the internal regular grid on one side,
     * this work just for the part of the stretch not the rim
     */

    ///< always the lowest integer
    n_high = (distance_to_inner + epstest) / delta_inner;

    ///< only for the stretched part take out the rim part
    if (n_high > n_stretched / 2.) {
        n_high_st  = (n_stretched / 2.);
        n_high_rim = n_high - n_high_st;
        p_rem      = 0;
        p_rem_low  = std::fmod((distance_to_inner + epsrem), delta_inner);
    }
    else {
        n_high_st  = n_high;
        n_high_rim = 0;
        ///< part remaining, use modulo
        p_rem     = std::fmod((distance_to_inner + epsrem), delta_inner);
        p_rem_low = 0.;
    }

    ///< computation of the new stretched delta integer part
    double delta = delta_inner;
    ///< initialization if is not using the cycle (first interval)
    double delta_last = delta;
    double deltacheck = 0;
    /*
     * using difference in delta for stretch
     * The point stretched is not in the regular grid and took out the points for the rim
     */
    for (int i = 0; i < n_high_st; i += 1) {
        delta_last = delta * new_ratio;
        delta_add  = delta_last - delta_inner;
        delta      = delta_last;
        deltacheck += delta_add;
    }
    ///< recomputation of point for every interval
    if (point > inner_start) {
        /*
         * after the end of the internal high resolution grid
         * no delta_add in the last points as they are rim
         */
        point += deltacheck;
    }
    else {
        /*
         * before the begin of the internal high resolution grid
         * no delta_add in the first points as they are rim
         */
        point -= deltacheck;
    }

    ///< SECTION 3 last part of stretch adding the remaing non integer part with the same ratio as in the stretching
    double delta_r    = p_rem * std::pow(new_ratio, (n_high_st + 1));
    double delta_addr = delta_r - p_rem;

    if (point > inner_start) {
        point += delta_addr;
    }
    else {
        point -= delta_addr;
    }
    ///< SECTION 4 rim area
    if (n_high > n_stretched / 2.) {
        double delta_l_h_ = 0;
        for (int i = 0; i < n_high_rim; i += 1) {
            delta_l_h_ += (delta_outer - delta_inner);
        }
        if (point > inner_start) {
            point += (delta_l_h_ + p_rem_low * (delta_outer - delta_inner));
        }
        else {
            point -= (delta_l_h_ + p_rem_low * (delta_outer - delta_inner));
        }
    }

    return normalised(point);
}


/**
 * Inverse of general stretch from a point in variable grid to the
 * a correspective point in a regular grid
 */

template <typename Rotation>
double VariableResolutionProjectionT<Rotation>::general_stretch_inv(const double crd, const bool L_long,
                                                                    const int n_stretched) const {
    constexpr float epstest =
        std::numeric_limits<float>::epsilon();  ///< correction used to change from double to integer

    double inner_size;       ///< number of new internal regular grid in double
    double inner_start;      ///< start of the regular grid
    double inner_end;        ///< end of the regular grid
    double point_st  = crd;  ///< input point in the variational grid
    double point_reg = 0.;   ///< point in the regular grid


    auto normalised = [L_long](double p) {
        if (L_long) {
            p = (p < 180) ? p + 360.0 : p;
        }
        return p;
    };

    point_st = normalised(point_st);

    /*
     * regular grids
     */
    if (var_ratio_ == 1) {
        point_reg = point_st;
        return normalised(point_reg);
    }

    /*
     *  SECTION 1
     *  INTERNAL REGULAR GRID the point is mapped to the same point
     */


    if (L_long) {
        inner_start = x_reg_start_;
        inner_size  = lam_hires_size;
    }
    else {
        inner_start = y_reg_start_;
        inner_size  = phi_hires_size;
    }

    inner_end = inner_start + inner_size;

    //< points in the internal regular grid
    if ((point_st >= inner_start - epstest) && (point_st <= inner_end + epstest)) {
        point_reg = point_st;
        return normalised(point_reg);
    }


    /* SECTION 2
     * Start variable res 'STRETCHED'
     * The point is mapped in a regular grid
     * from stretch internal grid: delta_dist
     * simply using delta_high
     */

    double distance_to_inner;  ///< distance from point to reg. grid
    double delta_add;          ///< additional part in stretch different from internal high resolution
    double new_ratio = new_ratio_[L_long ? 1 : 0];

    if (point_st < inner_start) {
        distance_to_inner = inner_start - point_st;
    }
    else if (point_st > inner_end) {
        distance_to_inner = point_st - inner_end;
    }

    /*
     * number of high resolution points intervals, that are not
     * in the internal regular grid on one side,
     * this work just for the part of the stretch not the rim
     */

    ///< computation of the new stretched delta integer part
    double delta = delta_inner;
    ///< initialization if is not using the cycle (first interval)
    double delta_last = delta;
    double deltacheck = 0;
    /*
     * using difference in delta for stretch
     * The point stretched is not in the regular grid and took out the points for the rim
     */

    double point_var = 0.;
    for (int i = 1; i < n_stretched / 2.; i += 1) {
        delta_last = delta * new_ratio;
        delta_add  = delta_last - delta_inner;
        delta      = delta_last;
        deltacheck += delta_add;
        if (point_st > inner_start) {
            point_reg = inner_end + (delta_inner * i);
            point_var = point_reg + deltacheck;
            if ((point_st <= point_var + epstest) && (point_st >= point_var - epstest)) {
                return normalised(point_reg);
            }
        }
        else {
            point_reg = inner_start - (delta_inner * i);
            point_var = point_reg - deltacheck;
            if ((point_st <= point_var + epstest) && (point_st >= point_var - epstest)) {
                return normalised(point_reg);
            }
        }
    }

    ///< SECTION 3_inv rim area
    int n_rim;  ///< number of rim points
    if (point_st > point_var) {
        n_rim     = (point_st - point_var) / delta_outer;
        point_reg = inner_end + (delta_inner * (n_stretched / 2 + n_rim));
        return normalised(point_reg);
    }
    else {
        if (point_st < point_var) {
            n_rim     = (point_var - point_st) / delta_outer;
            point_reg = inner_start - (delta_inner * (n_stretched / 2 + n_rim));
            return normalised(point_reg);
        }
    }
    return normalised(point_reg);
}


///< xy unstretched, only unrotation
template <typename Rotation>
void VariableResolutionProjectionT<Rotation>::lonlat2xy(double crd[]) const {
    ///< unrotate
    rotation_.unrotate(crd);
    //< normalization
    crd[0] = (crd[0] < 0) ? crd[0] + 360.0 : crd[0];
    ///< PUT the unstretch
    crd[0] = general_stretch_inv(crd[0], true, nx_stretched);
    crd[1] = general_stretch_inv(crd[1], false, ny_stretched);
}


///< From unstretched to stretched
template <typename Rotation>
void VariableResolutionProjectionT<Rotation>::xy2lonlat(double crd[]) const {
    /** PUT here the stretch for a point that come input
    * give number of power as input,work it out using as from the start of regular grid.
    * atlas::PointXY unstretchedXY = crd;
    * atlas::PointXY stretchedXY = stretch_LAM_gen(unstretchedXY);
    * stretch_LAM_gen(crd[]);
    */

    crd[0] = general_stretch(crd[0], true, nx_stretched);
    crd[1] = general_stretch(crd[1], false, ny_stretched);

    rotation_.rotate(crd);
}

template <typename Rotation>
ProjectionImpl::Jacobian VariableResolutionProjectionT<Rotation>::jacobian(const PointLonLat&) const {
    throw_NotImplemented("VariableResolution::jacobian", Here());
}

template <typename Rotation>
void VariableResolutionProjectionT<Rotation>::hash(eckit::Hash& hsh) const {
    hsh.add(static_type());
    rotation_.hash(hsh);
    //hsh.add( radius_ );
}

template class VariableResolutionProjectionT<NotRotated>;
template class VariableResolutionProjectionT<Rotated>;

namespace {
static ProjectionBuilder<VariableResolutionProjection> register_1(VariableResolutionProjection::static_type());
static ProjectionBuilder<RotatedVariableResolutionProjection> register_2(
    RotatedVariableResolutionProjection::static_type());
}  // namespace

}  // namespace detail
}  // namespace projection
}  // namespace atlas
