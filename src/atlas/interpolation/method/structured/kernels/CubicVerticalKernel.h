/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction. and Interpolation
 */

#pragma once
#include <cmath>

#include "atlas/grid/Stencil.h"
#include "atlas/grid/StencilComputer.h"
#include "atlas/grid/Vertical.h"
#include "atlas/library/config.h"
#include "atlas/util/Point.h"

namespace atlas {
namespace interpolation {
namespace method {

class CubicVerticalKernel {
    grid::ComputeVerticalStencil compute_vertical_stencil_;
    Vertical vertical_;
    static constexpr idx_t stencil_width() { return 4; }
    static constexpr idx_t stencil_size() { return stencil_width() * stencil_width(); }
    idx_t first_level_;
    idx_t last_level_;
    bool limiter_;

    static constexpr bool match_IFS() { return true; }

public:
    CubicVerticalKernel() = default;

    CubicVerticalKernel(const Vertical& vertical, const eckit::Configuration& config = util::NoConfig()):
        compute_vertical_stencil_(vertical, stencil_width()),
        vertical_(vertical),
        first_level_(vertical_.k_begin()),
        last_level_(vertical_.k_end() - 1) {
        limiter_ = config.getBool("limiter", false);
    }
    struct Weights {
        std::array<double, 4> weights_k;
    };
    using Stencil = grid::VerticalStencil<4>;

    template <typename stencil_t>
    void compute_stencil(const double z, stencil_t& stencil) const {
        compute_vertical_stencil_(z, stencil);
    }

    template <typename stencil_t, typename weights_t>
    void compute_weights(const double z, const stencil_t& stencil, weights_t& weights) const {
        auto& w = weights.weights_k;

        std::array<double, 4> zvec;
        for (idx_t k = 0; k < 4; ++k) {
            zvec[k] = vertical_(stencil.k(k));
        }

        //        auto quadratic_interpolation = [z]( const double zvec[], double w[] ) {
        //            double d01 = zvec[0] - zvec[1];
        //            double d02 = zvec[0] - zvec[2];
        //            double d12 = zvec[1] - zvec[2];
        //            double dc0 = d01 * d02;
        //            double dc1 = -d01 * d12;
        //            double d0  = z - zvec[0];
        //            double d1  = z - zvec[1];
        //            double d2  = z - zvec[2];
        //            w[0]       = ( d1 * d2 ) / dc0;
        //            w[1]       = ( d0 * d2 ) / dc1;
        //            w[2]       = 1. - w[0] - w[1];
        //        };

        auto cubic_interpolation = [z](const double zvec[], double w[]) {
            double d01 = zvec[0] - zvec[1];
            double d02 = zvec[0] - zvec[2];
            double d03 = zvec[0] - zvec[3];
            double d12 = zvec[1] - zvec[2];
            double d13 = zvec[1] - zvec[3];
            double d23 = zvec[2] - zvec[3];
            double dc0 = d01 * d02 * d03;
            double dc1 = -d01 * d12 * d13;
            double dc2 = d02 * d12 * d23;

            double d0 = z - zvec[0];
            double d1 = z - zvec[1];
            double d2 = z - zvec[2];
            double d3 = z - zvec[3];

            w[0] = (d1 * d2 * d3) / dc0;
#if defined(_CRAYC) && ATLAS_BUILD_TYPE_RELEASE
            // prevents FE_INVALID somehow (tested with Cray 8.7)
            ATLAS_ASSERT(!std::isnan(w[0]));
#endif
            w[1] = (d0 * d2 * d3) / dc1;
            w[2] = (d0 * d1 * d3) / dc2;
            w[3] = 1. - w[0] - w[1] - w[2];
        };

        if (stencil.k_interval() == 1) {
            // lev(k+0)   lev(k+1)   lev(k+2)   lev(k+3)
            //    |          |     x    |          |
            cubic_interpolation(zvec.data(), w.data());
            return;
        }

        if (stencil.k_interval() == 0) {
            if (match_IFS()) {
                // linear interpolation
                // lev0   lev1   lev2   lev3
                //  |  +   |      |      |
                //               w=0    w=0
                const double alpha = (zvec[1] - z) / (zvec[1] - zvec[0]);
                w[0]               = alpha;
                w[1]               = 1. - alpha;
                w[2]               = 0.;
                w[3]               = 0.;
            }
            else {
                cubic_interpolation(zvec.data(), w.data());
            }
            return;
        }

        if (stencil.k_interval() == 2) {
            if (match_IFS()) {
                // linear interpolation
                // lev(n-4)  lev(n-3)  lev(n-2)  lev(n-1)
                //   |         |         |    +    |
                //  w=0       w=0
                const double alpha = (zvec[3] - z) / (zvec[3] - zvec[2]);
                w[0]               = 0.;
                w[1]               = 0.;
                w[2]               = alpha;
                w[3]               = 1. - alpha;
            }
            else {
                cubic_interpolation(zvec.data(), w.data());
            }
            return;
        }

        if (stencil.k_interval() == -1) {
            // constant extrapolation
            //        lev0   lev1   lev2   lev3
            //      +  |------X------X------X
            //        w=1    w=0    w=0    w=0
            w[0] = 1.;
            w[1] = 0.;
            w[2] = 0.;
            w[3] = 0.;
            return;
        }

        if (stencil.k_interval() == 3) {
            // constant extrapolation
            //   lev(n-4)  lev(n-3)  lev(n-2)  lev(n-1)
            //      X---------X---------X---------|   +
            //     w=0      w=0       w=0       w=1
            w[0] = 0.;
            w[1] = 0.;
            w[2] = 0.;
            w[3] = 1.;
            return;
        }
        ATLAS_NOTIMPLEMENTED;  // should never be here
    }

    template <typename stencil_t, typename weights_t, typename array_t>
    void interpolate(const stencil_t& stencil, const weights_t& weights, const array_t& input, double& output) const {
        output        = 0.;
        const auto& w = weights.weights_k;
        for (idx_t k = 0; k < stencil_width(); ++k) {
            output += w[k] * input[stencil.k(k)];
        }


        if (limiter_) {
            idx_t k = stencil.k_interval();
            idx_t k1, k2;
            if (k < 1) {
                k1 = 0;
                k2 = 1;
            }
            else if (k > 1) {
                k1 = 2;
                k2 = 3;
            }
            else {
                k1 = 1;
                k2 = 2;
            }
            double f1     = input[stencil.k(k1)];
            double f2     = input[stencil.k(k2)];
            double maxval = std::max(f1, f2);
            double minval = std::min(f1, f2);
            output        = std::min(maxval, std::max(minval, output));
        }
    }

    template <typename array_t>
    double operator()(const double z, const array_t& input) const {
        Stencil stencil;
        compute_vertical_stencil_(z, stencil);
        Weights weights;
        compute_weights(z, stencil, weights);
        double output;
        interpolate(stencil, weights, input, output);
        return output;
    }
};

}  // namespace method
}  // namespace interpolation
}  // namespace atlas
