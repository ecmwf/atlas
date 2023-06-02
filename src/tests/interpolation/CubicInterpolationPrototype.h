/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "eckit/linalg/SparseMatrix.h"
#include "eckit/types/Types.h"

#include "atlas/array/ArrayView.h"
#include "atlas/array/MakeView.h"
#include "atlas/functionspace/StructuredColumns.h"
#include "atlas/grid/Grid.h"
#include "atlas/grid/Stencil.h"
#include "atlas/grid/StencilComputer.h"
#include "atlas/grid/Vertical.h"
#include "atlas/runtime/Exception.h"

using namespace eckit;
using namespace atlas::functionspace;
using namespace atlas::util;

namespace atlas {
namespace test {

//-----------------------------------------------------------------------------

class CubicVerticalInterpolation {
    grid::ComputeVerticalStencil compute_vertical_stencil_;
    Vertical vertical_;
    static constexpr idx_t stencil_width() { return 4; }
    static constexpr idx_t stencil_size() { return stencil_width() * stencil_width(); }
    idx_t first_level_;
    idx_t last_level_;
    bool limiter_;

public:
    CubicVerticalInterpolation(const Vertical& vertical, const eckit::Configuration& config = util::NoConfig()):
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
        else if (stencil.k_interval() == 3) {
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
        //        else if ( stencil.k_interval() == 0 ) {
        //            // quadratic interpolation
        //            // lev0   lev1   lev2   lev3
        //            //  |  +   |      |      |
        //            //                      w=0
        //            quadratic_interpolation( zvec.data(), w.data() );
        //            w[3] = 0.;
        //            return;
        //        }
        //        else if ( stencil.k_interval() == 2 ) {
        //            // quadratic interpolation
        //            // lev(n-4)  lev(n-3)  lev(n-2)  lev(n-1)
        //            //   |         |         |    +    |
        //            //  w=0
        //            quadratic_interpolation( zvec.data() + 1, w.data() + 1 );
        //            w[0] = 0.;
        //            return;
        //        }

        // cubic interpolation
        // lev(k+0)   lev(k+1)   lev(k+2)   lev(k+3)
        //    |          |     x    |          |
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
#if defined(_CRAYC) || defined(__NVCOMPILER) && ATLAS_BUILD_TYPE_RELEASE
        // prevents FE_INVALID somehow (tested with Cray 8.7)
        ATLAS_ASSERT(!std::isnan(w[0]));
#endif
        w[1] = (d0 * d2 * d3) / dc1;
        w[2] = (d0 * d1 * d3) / dc2;
        w[3] = 1. - w[0] - w[1] - w[2];
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
            if (k < 0) {
                k1 = k2 = 0;
            }
            else if (k > 2) {
                k1 = k2 = 3;
            }
            else {
                k1 = k;
                k2 = k + 1;
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
        grid::VerticalStencil<stencil_width()> stencil;
        compute_vertical_stencil_(z, stencil);
        Weights weights;
        compute_weights(z, stencil, weights);
        double output;
        interpolate(stencil, weights, input, output);
        return output;
    }
};

//-----------------------------------------------------------------------------

class CubicHorizontalInterpolation {
    functionspace::StructuredColumns fs_;
    grid::ComputeHorizontalStencil compute_horizontal_stencil_;
    static constexpr idx_t stencil_width() { return 4; }
    static constexpr idx_t stencil_size() { return stencil_width() * stencil_width(); }
    bool limiter_{false};

public:
    using Stencil = grid::HorizontalStencil<4>;

public:
    CubicHorizontalInterpolation(const functionspace::StructuredColumns& fs):
        fs_(fs), compute_horizontal_stencil_(fs.grid(), stencil_width()) {}
    template <typename weights_t>
    void compute_weights(const double x, const double y, weights_t& weights) const {
        grid::HorizontalStencil<stencil_width()> stencil;
        compute_horizontal_stencil_(x, y, stencil);
        compute_weights(x, y, stencil, weights);
    }

    struct Weights {
        std::array<std::array<double, 4>, 4> weights_i;
        std::array<double, 4> weights_j;
    };

    template <typename stencil_t>
    void compute_stencil(const double x, const double y, stencil_t& stencil) const {
        compute_horizontal_stencil_(x, y, stencil);
    }

    template <typename stencil_t, typename weights_t>
    void compute_weights(const double x, const double y, const stencil_t& stencil, weights_t& weights) const {
        PointXY P1, P2;
        std::array<double, 4> yvec;
        for (idx_t j = 0; j < stencil_width(); ++j) {
            auto& weights_i = weights.weights_i[j];
            fs_.compute_xy(stencil.i(1, j), stencil.j(j), P1);
            fs_.compute_xy(stencil.i(2, j), stencil.j(j), P2);
            double alpha               = (P2.x() - x) / (P2.x() - P1.x());
            double alpha_sqr           = alpha * alpha;
            double two_minus_alpha     = 2. - alpha;
            double one_minus_alpha_sqr = 1. - alpha_sqr;
            weights_i[0]               = -alpha * one_minus_alpha_sqr / 6.;
            weights_i[1]               = 0.5 * alpha * (1. + alpha) * two_minus_alpha;
            weights_i[2]               = 0.5 * one_minus_alpha_sqr * two_minus_alpha;
            weights_i[3]               = 1. - weights_i[0] - weights_i[1] - weights_i[2];
            yvec[j]                    = P1.y();
        }
        double dl12 = yvec[0] - yvec[1];
        double dl13 = yvec[0] - yvec[2];
        double dl14 = yvec[0] - yvec[3];
        double dl23 = yvec[1] - yvec[2];
        double dl24 = yvec[1] - yvec[3];
        double dl34 = yvec[2] - yvec[3];
        double dcl1 = dl12 * dl13 * dl14;
        double dcl2 = -dl12 * dl23 * dl24;
        double dcl3 = dl13 * dl23 * dl34;

        double dl1 = y - yvec[0];
        double dl2 = y - yvec[1];
        double dl3 = y - yvec[2];
        double dl4 = y - yvec[3];

        auto& weights_j = weights.weights_j;
        weights_j[0]    = (dl2 * dl3 * dl4) / dcl1;
#if defined(_CRAYC) || defined(__NVCOMPILER) && ATLAS_BUILD_TYPE_RELEASE
        // prevents FE_INVALID somehow (tested with Cray 8.7)
        ATLAS_ASSERT(!std::isnan(weights_j[0]));
#endif
        weights_j[1] = (dl1 * dl3 * dl4) / dcl2;
        weights_j[2] = (dl1 * dl2 * dl4) / dcl3;
        weights_j[3] = 1. - weights_j[0] - weights_j[1] - weights_j[2];
    }

    template <typename stencil_t, typename weights_t, typename array_t>
    void interpolate(const stencil_t& stencil, const weights_t& weights, const array_t& input, double& output) {
        std::array<std::array<idx_t, stencil_width()>, stencil_width()> index;
        const auto& weights_j = weights.weights_j;
        output                = 0.;
        for (idx_t j = 0; j < stencil_width(); ++j) {
            const auto& weights_i = weights.weights_i[j];
            for (idx_t i = 0; i < stencil_width(); ++i) {
                idx_t n = fs_.index(stencil.i(i, j), stencil.j(j));
                output += weights_i[i] * weights_j[j] * input[n];
                index[j][i] = n;
            }
        }

        if (limiter_) {
            limit(output, index, input);
        }
    }

    template <typename array_t>
    void limit(double& output, const std::array<std::array<idx_t, 4>, 4>& index, const array_t& input) {
        // Limit output to max/min of values in stencil marked by '*'
        //         x        x        x         x
        //              x     *-----*     x
        //                   /   P  |
        //          x       *------ *        x
        //        x        x        x         x
        double maxval = std::numeric_limits<double>::lowest();
        double minval = std::numeric_limits<double>::max();
        for (idx_t j = 1; j < 3; ++j) {
            for (idx_t i = 1; i < 3; ++i) {
                idx_t n    = index[j][i];
                double val = input[n];
                maxval     = std::max(maxval, val);
                minval     = std::min(minval, val);
            }
        }
        output = std::min(maxval, std::max(minval, output));
    }


    template <typename array_t>
    double operator()(const double x, const double y, const array_t& input) {
        grid::HorizontalStencil<stencil_width()> stencil;
        compute_horizontal_stencil_(x, y, stencil);
        Weights weights;
        compute_weights(x, y, stencil, weights);
        double output;
        interpolate(stencil, weights, input, output);
        return output;
    }

    struct WorkSpace {
        grid::HorizontalStencil<4> stencil;
        Weights weights;
    };

    using Triplet  = eckit::linalg::Triplet;
    using Triplets = std::vector<Triplet>;

    // Thread private workspace
    Triplets compute_triplets(const idx_t row, const double x, const double y, WorkSpace& ws) {
        Triplets triplets;
        triplets.reserve(stencil_size());
        insert_triplets(row, x, y, triplets, ws);
        return triplets;
    }

    Triplets reserve_triplets(size_t N) {
        Triplets triplets;
        triplets.reserve(N * stencil_size());
        return triplets;
    }

    void insert_triplets(const idx_t row, const double x, const double y, Triplets& triplets, WorkSpace& ws) {
        compute_horizontal_stencil_(x, y, ws.stencil);
        compute_weights(x, y, ws.stencil, ws.weights);
        const auto& wj = ws.weights.weights_j;
        for (idx_t j = 0; j < stencil_width(); ++j) {
            const auto& wi = ws.weights.weights_i[j];
            for (idx_t i = 0; i < stencil_width(); ++i) {
                idx_t col = fs_.index(ws.stencil.i(i, j), ws.stencil.j(j));
                double w  = wi[i] * wj[j];
                triplets.emplace_back(row, col, w);
            }
        }
    }
};

//-----------------------------------------------------------------------------

class Cubic3DInterpolation {
    functionspace::StructuredColumns fs_;
    CubicHorizontalInterpolation horizontal_interpolation_;
    CubicVerticalInterpolation vertical_interpolation_;
    static constexpr idx_t stencil_width() { return 4; }
    static constexpr idx_t stencil_size() { return stencil_width() * stencil_width(); }

public:
    struct Weights {
        std::array<std::array<double, 4>, 4> weights_i;
        std::array<double, 4> weights_j;
        std::array<double, 4> weights_k;
    };

    using Stencil = grid::Stencil3D<4>;

    Cubic3DInterpolation(const functionspace::StructuredColumns& fs):
        fs_(fs), horizontal_interpolation_(fs), vertical_interpolation_(fs.vertical()) {}

    template <typename stencil_t>
    void compute_stencil(const double x, const double y, const double z, stencil_t& stencil) const {
        horizontal_interpolation_.compute_stencil(x, y, stencil);
        vertical_interpolation_.compute_stencil(z, stencil);
    }

    template <typename weights_t>
    void compute_weights(const double x, const double y, const double z, weights_t& weights) const {
        Stencil stencil;
        compute_stencil(x, y, z, stencil);
        compute_weights(x, y, z, stencil, weights);
    }

    template <typename stencil_t, typename weights_t>
    void compute_weights(const double x, const double y, const double z, const stencil_t& stencil,
                         weights_t& weights) const {
        horizontal_interpolation_.compute_weights(x, y, stencil, weights);
        vertical_interpolation_.compute_weights(z, stencil, weights);
    }

    template <typename stencil_t, typename weights_t, typename array_t>
    void interpolate(const stencil_t& stencil, const weights_t& weights, const array_t& input, double& output) {
        output         = 0.;
        const auto& wj = weights.weights_j;
        const auto& wk = weights.weights_k;

        for (idx_t j = 0; j < stencil_width(); ++j) {
            const auto& wi = weights.weights_i[j];
            for (idx_t i = 0; i < stencil_width(); ++i) {
                idx_t n = fs_.index(stencil.i(i, j), stencil.j(j));
                for (idx_t k = 0; k < stencil_width(); ++k) {
                    output += wi[i] * wj[j] * wk[k] * input(n, stencil.k(k));
                }
            }
        }
    }

    template <typename array_t>
    double operator()(const double x, const double y, const double z, const array_t& input) {
        Stencil stencil;
        Weights weights;
        compute_stencil(x, y, z, stencil);
        compute_weights(x, y, z, stencil, weights);
        double output;
        interpolate(stencil, weights, input, output);
        return output;
    }

    template <typename point_t, typename array_t>
    double operator()(const point_t& p, const array_t& input) {
        Stencil stencil;
        Weights weights;
        compute_stencil(p[0], p[1], p[2], stencil);
        compute_weights(p[0], p[1], p[2], stencil, weights);
        double output;
        interpolate(stencil, weights, input, output);
        return output;
    }
};

//-----------------------------------------------------------------------------

}  // namespace test
}  // namespace atlas
