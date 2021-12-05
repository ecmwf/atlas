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
#include "atlas/util/Point.h"

namespace atlas {
namespace interpolation {
namespace method {

class LinearVerticalKernel {
    grid::ComputeVerticalStencil compute_vertical_stencil_;
    Vertical vertical_;
    static constexpr idx_t stencil_width() { return 2; }
    static constexpr idx_t stencil_size() { return stencil_width() * stencil_width(); }
    idx_t first_level_;
    idx_t last_level_;

public:
    LinearVerticalKernel() = default;

    LinearVerticalKernel(const Vertical& vertical, const eckit::Configuration&):
        compute_vertical_stencil_(vertical, stencil_width()),
        vertical_(vertical),
        first_level_(vertical_.k_begin()),
        last_level_(vertical_.k_end() - 1) {}
    struct Weights {
        std::array<double, 2> weights_k;
    };
    using Stencil = grid::VerticalStencil<2>;

    template <typename stencil_t>
    void compute_stencil(const double z, stencil_t& stencil) const {
        compute_vertical_stencil_(z, stencil);
    }

    template <typename stencil_t, typename weights_t>
    void compute_weights(const double z, const stencil_t& stencil, weights_t& weights) const {
        auto& w = weights.weights_k;

        std::array<double, stencil_width()> zvec;
        for (idx_t k = 0; k < stencil_width(); ++k) {
            zvec[k] = vertical_(stencil.k(k));
        }

        if (stencil.k_interval() == -1) {
            // constant extrapolation
            //        lev0   lev1
            //      +  |------X
            //        w=1    w=0
            w[0] = 1.;
            w[1] = 0.;
            return;
        }
        else if (stencil.k_interval() == 1) {
            // constant extrapolation
            //   lev(n-2)  lev(n-1)
            //      X---------|   +
            //     w=0       w=1
            w[0] = 0.;
            w[1] = 1.;
            return;
        }

        // Linear interpolation
        // lev(k+0)   lev(k+1)
        //    |     x    |
        const double alpha = (zvec[1] - z) / (zvec[1] - zvec[0]);
        w[0]               = alpha;
        w[1]               = 1. - alpha;
    }

    template <typename stencil_t, typename weights_t, typename array_t>
    void interpolate(const stencil_t& stencil, const weights_t& weights, const array_t& input, double& output) const {
        output        = 0.;
        const auto& w = weights.weights_k;
        for (idx_t k = 0; k < stencil_width(); ++k) {
            output += w[k] * input[stencil.k(k)];
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
