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
#include "atlas/projection/detail/StretchRegular.h"
#include "atlas/runtime/Exception.h"
#include "atlas/util/Config.h"
#include "atlas/util/Constants.h"
#include "eckit/testing/Test.h"

/**
* Projection for LAM stretching
*
* The theory used in this code can be found under:
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

namespace atlas {
namespace projection {
namespace detail {


///< specification parameters
template <typename Rotation>
typename StretchLAM<Rotation>::Spec StretchLAM<Rotation>::spec() const {
    Spec proj_st;
    proj_st.set("delta_low", delta_low_);      ///< resolution of the external regular grid (rim)
    proj_st.set("delta_high", delta_high_);    ///< resolution of the regional model (regular grid)
    proj_st.set("var_ratio", var_ratio_);      ///< power used for the stretching
    proj_st.set("x_reg_start", x_reg_start_);  ///< xstart of the internal regional grid
    proj_st.set("x_reg_end", x_reg_end_);      ///< xend of the internal regional grid
    proj_st.set("y_reg_start", y_reg_start_);  ///< ystart of the internal regional grid
    proj_st.set("y_reg_end", y_reg_end_);      ///< yend of the internal regional grid
    proj_st.set("startx", startx_);            ///< original domain startx
    proj_st.set("endx", endx_);                ///< original domain endx
    proj_st.set("starty", starty_);            ///< original domain starty
    proj_st.set("endy", endy_);                ///< original domain endy
    proj_st.set("rim_widthx", rim_widthx_);    ///< xsize of the rim
    proj_st.set("rim_widthy", rim_widthy_);    ///< ysize of the rim

    return proj_st;
}


///< constructors
template <typename Rotation>
StretchLAM<Rotation>::StretchLAM(const eckit::Parametrisation& proj_st): ProjectionImpl(), rotation_(proj_st) {
    proj_st.get("delta_low", delta_low_ = 0.);      ///< resolution of the external regular grid (rim)
    proj_st.get("delta_hi", delta_high_ = 0.);      ///< resolution of the regional model (regular grid)
    proj_st.get("var_ratio", var_ratio_ = 0.);      ///< power used for the stretching
    proj_st.get("x_reg_start", x_reg_start_ = 0.);  ///< xstart of the internal regional grid
    proj_st.get("y_reg_start", y_reg_start_ = 0.);  ///< ystart of the internal regional grid
    proj_st.get("x_reg_end", x_reg_end_ = 0.);      ///< xend of the regular part of stretched internal grid
    proj_st.get("y_reg_end", y_reg_end_ = 0.);      ///< yend of the regular part of stretched internal grid
    proj_st.get("startx", startx_ = 0.);            ///< original domain startx
    proj_st.get("endx", endx_ = 0.);                ///< original domain endx
    proj_st.get("starty", starty_ = 0.);            ///< original domain starty
    proj_st.get("endy", endy_ = 0.);                ///< original domain endy
    proj_st.get("rim_widthx", rim_widthx_);         ///< xsize of the rim
    proj_st.get("rim_widthy", rim_widthy_);         ///< ysize of the rim


    constexpr float epsilon = std::numeric_limits<float>::epsilon();  ///< value used to check if the values are equal
    constexpr float epstest =
        std::numeric_limits<float>::epsilon();  ///< correction used to change from double to intege

    ///< original domain size includes the points for the rim
    deltax_all = (endx_ - startx_);
    deltay_all = (endy_ - starty_);

    n_stretchedx_ = 0;
    n_stretchedy_ = 0;
    n_x_rim_      = 0;
    n_y_rim_      = 0;
    nx_           = 0;
    ny_           = 0;

    if (var_ratio_ == 1) {
        lam_hires_size = deltax_all;
        phi_hires_size = deltay_all;
        lambda_start   = x_reg_start_;
        phi_start      = y_reg_start_;
    }
    else {
        lam_hires_size = x_reg_end_ - x_reg_start_;
        phi_hires_size = y_reg_end_ - y_reg_start_;

        /** check for start of the new regular LAM grid x and y
                      * in the middle of the previous regular grid
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
        n_x_rim_      = rim_widthx_ / delta_low_;
        n_y_rim_      = rim_widthy_ / delta_low_;
        n_stretchedx_ = ((deltax_all + epstest - lam_hires_size) / delta_high_) - n_x_rim_;
        n_stretchedy_ = ((deltay_all + epstest - phi_hires_size) / delta_high_) - n_y_rim_;
        lambda_start  = x_reg_start_;
        phi_start     = y_reg_start_;
        /**
                       *  check if stretched grid is in the middle
                       *  of the previous regular grid
                      */
        check_x  = startx_ + add_xf_ - lambda_start;
        check_y  = starty_ + add_yf_ - phi_start;
        check_st = n_stretchedx_ - n_stretchedy_;
        checkvalue(epsilon, check_x);
        checkvalue(epsilon, check_y);
        checkvalue(epsilon, check_st);
    }

    nx_ = ((deltax_all + epstest) / delta_high_) + 1;
    ny_ = ((deltay_all + epstest) / delta_high_) + 1;
}


template <typename Rotation>
void StretchLAM<Rotation>::checkvalue(const double& epsilon, const double& value_check) const {
    //err_message = "USER defined limits not in the middle of the area " + str;
    if (value_check > epsilon || value_check < (-1. * epsilon)) {
        std::string err_message;
        std::string str = std::to_string(value_check);
        throw eckit::BadValue("USER defined limits not in the middle of the area " + str, Here());
    }
}

/**
 * General stretch from a point in regular grid to the
 * a correspective point in a stretched grid
 */

template <typename Rotation>
double StretchLAM<Rotation>::general_stretch(const double crd, const bool L_long, int n_int, const int n_stretched_,
                                             const int n_rim_) const {
    double high_size;     ///< number of new internal regular grid in double
    double lamphi_start;  ///< start of the regular grid
    double lamphi_end;    ///< end of the regular grid
    double lamphi = crd;
    double point  = lamphi;  ///< starting point
    constexpr float epstest =
        std::numeric_limits<float>::epsilon();  ///< correction used to change from double to integer
    constexpr double epsrem =
        0.1 * std::numeric_limits<double>::epsilon() /
        std::numeric_limits<float>::epsilon();  ///< correction used to part the find a part of an integer


    if ((n_int > 0) && (var_ratio_ > 1)) {
        double delta_dist;  ///< distance from point to reg. grid
        double delta_add;   ///< additional part in stretch different from internal high resolution
        int n_high;         ///< number of points, from point to reg. grid
        int n_high_st;      ///< number of stretched points, from point to reg grid
        int n_high_rim;     ///< number of rim points, from point to reg grid
        double p_rem;       ///< remaining part in stretch if delta_dist not multiple of delta_high_
        double p_rem_low;   ///< remaining part in rim if delta_dist not multiple of delta_high_

        n_int -= 1;  ///< input number of points whole grid, output intervals
        ///< n_rim_ are the number of points in the rim area
        high_size = (n_int - n_stretched_ - n_rim_) * delta_high_;
        ///< number of regular internal grid points, integer
        int ints_high = (high_size + epstest) / delta_high_;

        ///< number of variable (stretched) grid points in one side, integer
        int var_ints = (n_int + epstest - n_rim_ - ints_high) / 2.;
        /** compute ratio,
           *  change stretching factor so that high and low grids
           *   retain original sizes
           */
        ///< number of variable (stretched) grid points in one side,double
        double var_ints_f = (n_int - n_rim_ - ints_high) / 2.;
        double logr       = std::log(var_ratio_);
        double log_ratio  = (var_ints_f - 0.5) * logr;
        double new_ratio  = std::exp(log_ratio / var_ints);


        /**
           *  SECTION 1
           *  INTERNAL REGULAR GRID the point is mapped to the same point
           */

        if (L_long) {
            lamphi_start = x_reg_start_;
        }
        else {
            lamphi_start = y_reg_start_;
        }
        lamphi_end = lamphi_start + high_size;

        if ((point >= lamphi_start) || (point <= lamphi_end)) {
            lamphi = point;
        }

        /** SECTION 2
          * Start variable res 'STRETCHED'
          * The point is mapped in a stretched grid related to
          * distance from the regular grid and stretch internal grid: delta_dist
          */

        if (point < lamphi_start) {
            delta_dist = lamphi_start - point;
        }
        else if (point > lamphi_end) {
            delta_dist = point - lamphi_end;
        }

        if ((point < lamphi_start) || (point > lamphi_end)) {
            /**
              * number of high resolution points intervals, that are not
              * in the internal regular grid on one side,
              * this work just for the part of the stretch not the rim
              */

            ///< always the lowest integer
            n_high = (delta_dist + epstest) / delta_high_;

            ///< only for the stretched part take out the rim part
            if (n_high > n_stretched_ / 2.) {
                n_high_st  = (n_stretched_ / 2.);
                n_high_rim = n_high - n_high_st;
                p_rem      = 0;
                p_rem_low  = std::fmod((delta_dist + epsrem), delta_high_);
            }
            else {
                n_high_st  = n_high;
                n_high_rim = 0;
                ///< part remaining, use modulo
                p_rem     = std::fmod((delta_dist + epsrem), delta_high_);
                p_rem_low = 0.;
            }

            ///< computation of the new stretched delta integer part
            double delta = delta_high_;
            ///< initialization if is not using the cycle (first interval)
            double delta_last = delta;
            double deltacheck = 0;
            /**
             * using difference in delta for stretch
             * The point stretched is not in the regular grid and took out the points for the rim
             */
            for (int i = 0; i < n_high_st; i += 1) {
                delta_last = delta * new_ratio;
                delta_add  = delta_last - delta_high_;
                delta      = delta_last;
                deltacheck += delta_add;
            }
            ///< recomputation of point for every interval
            if (point > lamphi_start) {
                /**
                  * after the end of the internal high resolution grid
                  * no delta_add in the last points as they are rim
                  */
                point = point + deltacheck;
            }
            else {
                /**
                 * before the begin of the internal high resolution grid
                 * no delta_add in the first points as they are rim
                 */
                point = point - deltacheck;
            }

            ///< SECTION 3 last part of stretch adding the remaing non integer part with the same ratio as in the stretching
            double delta_r    = p_rem * std::pow(new_ratio, (n_high_st + 1));
            double delta_addr = delta_r - p_rem;

            if (point > lamphi_start) {
                point = point + delta_addr;
            }
            else {
                point = point - delta_addr;
            }
            ///< SECTION 4 rim area
            if (n_high > n_stretched_ / 2.) {
                double delta_l_h_ = 0;
                for (int i = 0; i < n_high_rim; i += 1) {
                    delta_l_h_ = delta_l_h_ + (delta_low_ - delta_high_);
                }
                if (point > lamphi_start) {
                    point = point + delta_l_h_ + p_rem_low * (delta_low_ - delta_high_);
                }
                else {
                    point = point - delta_l_h_ - p_rem_low * (delta_low_ - delta_high_);
                }
            }

            lamphi = point;
            lamphi = (lamphi >= 360) ? lamphi - 360.0 : lamphi;
        }


        /**
          * SECTION 5: Reset lat/lons
          * If longitude domain crosses the meridian (0 degrees) then
          * add 360 degrees to points East of the meridian
          */
        if (L_long) {
            lamphi = (lamphi < 180) ? lamphi + 360.0 : lamphi;
        }
    }
    else {
        /*
            * SECTION 6: regular grids
            * var_ratio = 1.0
            */
        lamphi = point;
        if (L_long) {
            lamphi = (lamphi < 180) ? lamphi + 360.0 : lamphi;
        }
    }
    return lamphi;
}

///< xy unstretched, only unrotation
template <typename Rotation>
void StretchLAM<Rotation>::lonlat2xy(double crd[]) const {
    ///<unrotate
    rotation_.rotate(crd);

    ///< PUT the unstretch, I don't have it now
    ATLAS_NOTIMPLEMENTED;
}

///< From unstretched to stretched
template <typename Rotation>
void StretchLAM<Rotation>::xy2lonlat(double crd[]) const {
    /** PUT here the stretch for a point that come input
    * give number of power as input,work it out using as from the start of regular grid.
    * atlas::PointXY unstretchedXY = crd;
    * atlas::PointXY stretchedXY = stretch_LAM_gen(unstretchedXY);
    * stretch_LAM_gen(crd[]);
    */


    crd[0] = general_stretch(crd[0], true, nx_, n_stretchedx_, n_x_rim_);
    crd[1] = general_stretch(crd[1], false, ny_, n_stretchedy_, n_y_rim_);

    ///< rotate
    rotation_.unrotate(crd);
}

template <typename Rotation>
ProjectionImpl::Jacobian StretchLAM<Rotation>::jacobian(const PointLonLat&) const {
    throw_NotImplemented("StretchRegularT::jacobian", Here());
}

template <typename Rotation>
void StretchLAM<Rotation>::hash(eckit::Hash& hsh) const {
    hsh.add(static_type());
    rotation_.hash(hsh);
    //hsh.add( radius_ );
}

template class StretchLAM<NotRotated>;
template class StretchLAM<Rotated>;

namespace {
static ProjectionBuilder<StretchRegular> register_1(StretchRegular::static_type());
static ProjectionBuilder<RotatedStretchRegular> register_2(RotatedStretchRegular::static_type());
}  // namespace

}  // namespace detail
}  // namespace projection
}  // namespace atlas
