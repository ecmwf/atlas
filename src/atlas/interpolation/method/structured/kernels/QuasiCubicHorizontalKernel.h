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
#include <limits>

#include "CubicHorizontalLimiter.h"

#include "eckit/linalg/Triplet.h"

#include "atlas/array/ArrayView.h"
#include "atlas/functionspace/StructuredColumns.h"
#include "atlas/grid/Stencil.h"
#include "atlas/grid/StencilComputer.h"
#include "atlas/library/config.h"
#include "atlas/runtime/Exception.h"
#include "atlas/util/CoordinateEnums.h"
#include "atlas/util/Point.h"

namespace atlas {
namespace interpolation {
namespace method {

class QuasiCubicHorizontalKernel {
    using Limiter  = CubicHorizontalLimiter;
    using Triplet  = eckit::linalg::Triplet;
    using Triplets = std::vector<Triplet>;

public:
    QuasiCubicHorizontalKernel() = default;

    QuasiCubicHorizontalKernel(const functionspace::StructuredColumns& fs,
                               const util::Config& config = util::NoConfig()) {
        src_ = fs;
        ATLAS_ASSERT(src_);
        ATLAS_ASSERT(src_.halo() >= 2);
        compute_horizontal_stencil_ = grid::ComputeHorizontalStencil(src_.grid(), stencil_width());
        limiter_                    = config.getBool("limiter", false);
    }

protected:
    functionspace::StructuredColumns src_;
    grid::ComputeHorizontalStencil compute_horizontal_stencil_;
    bool limiter_{false};

public:
    static std::string className() { return "CubicHorizontalKernel"; }
    static constexpr idx_t stencil_width() { return 4; }
    static constexpr idx_t stencil_size() { return 12; }
    static constexpr idx_t stencil_halo() {
        return static_cast<idx_t>(static_cast<double>(stencil_width()) / 2. + 0.5);
    }

public:
    using Stencil = grid::HorizontalStencil<4>;
    struct Weights {
        std::array<std::array<double, 4>, 4> weights_i;
        std::array<double, 4> weights_j;
    };

public:
    struct WorkSpace {
        Stencil stencil;
        Weights weights;
    };

    template <typename stencil_t>
    void compute_stencil(const double x, const double y, stencil_t& stencil) const {
        compute_horizontal_stencil_(x, y, stencil);
    }

    template <typename stencil_t>
    void make_valid_stencil(double& x, double y, stencil_t& stencil, bool retry = true) const {
        for (idx_t j = 0; j < stencil_width(); ++j) {
            idx_t imin = stencil.i(0, j);
            idx_t imax = stencil.i(stencil_width()-1, j);
            if (imin < src_.i_begin_halo(stencil.j(j))) {
                if (retry ) {
                    x += 360.;
                    compute_stencil(x, y, stencil);
                    return make_valid_stencil(x, y, stencil, false);
                }
                else {
                    Log::error() << "Stencil out of bounds" << std::endl;
                    ATLAS_THROW_EXCEPTION("stencil out of bounds");
                }
            }
            if (imax >= src_.i_end_halo(stencil.j(j))) {
                if (retry ) {
                    x -= 360.;
                    compute_stencil(x, y, stencil);
                    return make_valid_stencil(x, y, stencil, false);
                }
                else {
                    Log::error() << "Stencil out of bounds" << std::endl;
                    ATLAS_THROW_EXCEPTION("Stencil out of bounds");
                }
            }
        }

    }

    template <typename weights_t>
    void compute_weights(double x, double y, weights_t& weights) const {
        Stencil stencil;
        compute_stencil(x, y, stencil);
        compute_weights(x, y, stencil, weights);
    }

    template <typename stencil_t, typename weights_t>
    void compute_weights(const double x, const double y, const stencil_t& stencil, weights_t& weights) const {
        PointXY P1, P2;
        std::array<double, 4> yvec;

        // Compute x-direction weights LINEAR for outer rows  ( j = {0,3} )
        for (idx_t j = 0; j < 4; j += 3) {
            auto& weights_i = weights.weights_i[j];
            src_.compute_xy(stencil.i(1, j), stencil.j(j), P1);
            src_.compute_xy(stencil.i(2, j), stencil.j(j), P2);
            const double alpha = (P2.x() - x) / (P2.x() - P1.x());
            weights_i[1]       = alpha;
            weights_i[2]       = 1. - alpha;
            yvec[j]            = P1.y();
        }

        // Compute x-direction weights CUBIC for inner rows ( j = {1,2} )
        for (idx_t j = 1; j < 3; ++j) {
            auto& weights_i = weights.weights_i[j];
            src_.compute_xy(stencil.i(1, j), stencil.j(j), P1);
            src_.compute_xy(stencil.i(2, j), stencil.j(j), P2);
            const double alpha               = (P2.x() - x) / (P2.x() - P1.x());
            const double alpha_sqr           = alpha * alpha;
            const double two_minus_alpha     = 2. - alpha;
            const double one_minus_alpha_sqr = 1. - alpha_sqr;
            weights_i[0]                     = -alpha * one_minus_alpha_sqr / 6.;
            weights_i[1]                     = 0.5 * alpha * (1. + alpha) * two_minus_alpha;
            weights_i[2]                     = 0.5 * one_minus_alpha_sqr * two_minus_alpha;
            weights_i[3]                     = 1. - weights_i[0] - weights_i[1] - weights_i[2];
            yvec[j]                          = P1.y();
        }

        // Compute weights in y-direction
        const double dl12 = yvec[0] - yvec[1];
        const double dl13 = yvec[0] - yvec[2];
        const double dl14 = yvec[0] - yvec[3];
        const double dl23 = yvec[1] - yvec[2];
        const double dl24 = yvec[1] - yvec[3];
        const double dl34 = yvec[2] - yvec[3];
        const double dcl1 = dl12 * dl13 * dl14;
        const double dcl2 = -dl12 * dl23 * dl24;
        const double dcl3 = dl13 * dl23 * dl34;

        const double dl1 = y - yvec[0];
        const double dl2 = y - yvec[1];
        const double dl3 = y - yvec[2];
        const double dl4 = y - yvec[3];

        auto& weights_j = weights.weights_j;
        weights_j[0]    = (dl2 * dl3 * dl4) / dcl1;
#if defined(_CRAYC) && ATLAS_BUILD_TYPE_RELEASE
        // prevents FE_INVALID somehow (tested with Cray 8.7)
        ATLAS_ASSERT(!std::isnan(weights_j[0]));
#endif
        weights_j[1] = (dl1 * dl3 * dl4) / dcl2;
        weights_j[2] = (dl1 * dl2 * dl4) / dcl3;
        weights_j[3] = 1. - weights_j[0] - weights_j[1] - weights_j[2];
    }

    template <typename stencil_t, typename weights_t, typename array_t>
    typename array_t::value_type interpolate(const stencil_t& stencil, const weights_t& weights,
                                             const array_t& input) const {
        using Value = typename array_t::value_type;

        std::array<std::array<idx_t, stencil_width()>, stencil_width()> index;
        const auto& weights_j = weights.weights_j;
        Value output          = 0.;
        // LINEAR for outer rows  ( j = {0,3} )
        for (idx_t j = 0; j < 4; j += 3) {
            const auto& weights_i = weights.weights_i[j];
            for (idx_t i : {1, 2}) {
                idx_t n = src_.index(stencil.i(i, j), stencil.j(j));
                Value w = weights_i[i] * weights_j[j];
                output += w * input[n];
                index[j][i] = n;
            }
        }
        // CUBIC for inner rows ( j = {1,2} )
        for (idx_t j = 1; j < 3; ++j) {
            const auto& weights_i = weights.weights_i[j];
            for (idx_t i = 0; i < stencil_width(); ++i) {
                idx_t n = src_.index(stencil.i(i, j), stencil.j(j));
                Value w = weights_i[i] * weights_j[j];
                output += w * input[n];
                index[j][i] = n;
            }
        }

        if (limiter_) {
            Limiter::limit(output, index, input);
        }
        return output;
    }

    template <typename stencil_t, typename weights_t, typename Value, int Rank>
    typename std::enable_if<(Rank == 1), void>::type interpolate(const stencil_t& stencil, const weights_t& weights,
                                                                 const array::ArrayView<const Value, Rank>& input,
                                                                 array::ArrayView<Value, Rank>& output, idx_t r) const {
        std::array<std::array<idx_t, stencil_width()>, stencil_width()> index;
        const auto& weights_j = weights.weights_j;
        output(r)             = 0.;

        // LINEAR for outer rows  ( j = {0,3} )
        for (idx_t j = 0; j < 4; j += 3) {
            const auto& weights_i = weights.weights_i[j];
            for (idx_t i = 1; i < 3; ++i) {  // i = {1,2}
                idx_t n = src_.index(stencil.i(i, j), stencil.j(j));
                Value w = weights_i[i] * weights_j[j];
                output(r) += w * input[n];
                index[j][i] = n;
            }
        }
        // CUBIC for inner rows ( j = {1,2} )
        for (idx_t j = 1; j < 3; ++j) {
            const auto& weights_i = weights.weights_i[j];
            for (idx_t i = 0; i < stencil_width(); ++i) {
                idx_t n = src_.index(stencil.i(i, j), stencil.j(j));
                Value w = weights_i[i] * weights_j[j];
                output(r) += w * input[n];
                index[j][i] = n;
            }
        }

        if (limiter_) {
            Limiter::limit(index, input, output, r);
        }
    }

    template <typename stencil_t, typename weights_t, typename Value, int Rank>
    typename std::enable_if<(Rank == 2), void>::type interpolate(const stencil_t& stencil, const weights_t& weights,
                                                                 const array::ArrayView<const Value, Rank>& input,
                                                                 array::ArrayView<Value, Rank>& output, idx_t r) const {
        std::array<std::array<idx_t, stencil_width()>, stencil_width()> index;
        const auto& weights_j = weights.weights_j;
        const idx_t Nk        = output.shape(1);
        for (idx_t k = 0; k < Nk; ++k) {
            output(r, k) = 0.;
        }

        // LINEAR for outer rows  ( j = {0,3} )
        for (idx_t j = 0; j < 4; j += 3) {
            const auto& weights_i = weights.weights_i[j];
            for (idx_t i = 1; i < 3; ++i) {  // i = {1,2}
                idx_t n = src_.index(stencil.i(i, j), stencil.j(j));
                Value w = weights_i[i] * weights_j[j];
                for (idx_t k = 0; k < Nk; ++k) {
                    output(r, k) += w * input(n, k);
                }
                index[j][i] = n;
            }
        }
        // CUBIC for inner rows ( j = {1,2} )
        for (idx_t j = 1; j < 3; ++j) {
            const auto& weights_i = weights.weights_i[j];
            for (idx_t i = 0; i < stencil_width(); ++i) {
                idx_t n = src_.index(stencil.i(i, j), stencil.j(j));
                Value w = weights_i[i] * weights_j[j];
                for (idx_t k = 0; k < Nk; ++k) {
                    output(r, k) += w * input(n, k);
                }
                index[j][i] = n;
            }
        }

        if (limiter_) {
            Limiter::limit(index, input, output, r);
        }
    }


    template <typename stencil_t, typename weights_t, typename Value, int Rank>
    typename std::enable_if<(Rank == 3), void>::type interpolate(const stencil_t& stencil, const weights_t& weights,
                                                                 const array::ArrayView<const Value, Rank>& input,
                                                                 array::ArrayView<Value, Rank>& output, idx_t r) const {
        std::array<std::array<idx_t, stencil_width()>, stencil_width()> index;
        const auto& weights_j = weights.weights_j;
        const idx_t Nk        = output.shape(1);
        const idx_t Nl        = output.shape(2);
        for (idx_t k = 0; k < Nk; ++k) {
            for (idx_t l = 0; l < Nl; ++l) {
                output(r, k, l) = 0.;
            }
        }

        // LINEAR for outer rows  ( j = {0,3} )
        for (idx_t j = 0; j < 4; j += 3) {
            const auto& weights_i = weights.weights_i[j];
            for (idx_t i = 1; i < 3; ++i) {  // i = {1,2}
                idx_t n = src_.index(stencil.i(i, j), stencil.j(j));
                Value w = weights_i[i] * weights_j[j];
                for (idx_t k = 0; k < Nk; ++k) {
                    for (idx_t l = 0; l < Nl; ++l) {
                        output(r, k, l) += w * input(n, k, l);
                    }
                }
                index[j][i] = n;
            }
        }
        // CUBIC for inner rows ( j = {1,2} )
        for (idx_t j = 1; j < 3; ++j) {
            const auto& weights_i = weights.weights_i[j];
            for (idx_t i = 0; i < stencil_width(); ++i) {
                idx_t n = src_.index(stencil.i(i, j), stencil.j(j));
                Value w = weights_i[i] * weights_j[j];
                for (idx_t k = 0; k < Nk; ++k) {
                    for (idx_t l = 0; l < Nl; ++l) {
                        output(r, k, l) += w * input(n, k, l);
                    }
                }
                index[j][i] = n;
            }
        }

        if (limiter_) {
            Limiter::limit(index, input, output, r);
        }
    }



    template <typename array_t>
    typename array_t::value_type operator()(double x, double y, const array_t& input) const {
        Stencil stencil;
        compute_stencil(x, y, stencil);
        Weights weights;
        compute_weights(x, y, stencil, weights);
        make_valid_stencil(x, y, stencil);
        return interpolate(stencil, weights, input);
    }

    template <typename array_t>
    typename array_t::value_type interpolate(PointXY p, const array_t& input, WorkSpace& ws) const {
        compute_stencil(p.x(), p.y(), ws.stencil);
        compute_weights(p.x(), p.y(), ws.stencil, ws.weights);
        make_valid_stencil(p.x(), p.y(), ws.stencil);
        return interpolate(ws.stencil, ws.weights, input);
    }

    // Thread private workspace
    Triplets compute_triplets(const idx_t row, const double x, const double y, WorkSpace& ws) const {
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

    Triplets allocate_triplets(size_t N) { return Triplets(N * stencil_size()); }

    void insert_triplets(const idx_t row, const PointXY& p, Triplets& triplets, WorkSpace& ws) const {
        insert_triplets(row, p.x(), p.y(), triplets, ws);
    }

    void insert_triplets(const idx_t row, double x, double y, Triplets& triplets, WorkSpace& ws) const {
        compute_stencil(x, y, ws.stencil);
        compute_weights(x, y, ws.stencil, ws.weights);

        make_valid_stencil(x, y, ws.stencil);
        const auto& wj = ws.weights.weights_j;

        idx_t pos = row * stencil_size();

        // LINEAR for outer rows  ( j = {0,3} )
        for (idx_t j = 0; j < 4; j += 3) {
            const auto& wi = ws.weights.weights_i[j];
            for (idx_t i = 1; i < 3; ++i) {  // i = {1,2}
                idx_t col       = src_.index(ws.stencil.i(i, j), ws.stencil.j(j));
                double w        = wi[i] * wj[j];
                triplets[pos++] = Triplet(row, col, w);
            }
        }

        // CUBIC for inner rows ( j = {1,2} )
        for (idx_t j = 1; j < 3; ++j) {
            const auto& wi = ws.weights.weights_i[j];
            for (idx_t i = 0; i < stencil_width(); ++i) {
                idx_t col       = src_.index(ws.stencil.i(i, j), ws.stencil.j(j));
                double w        = wi[i] * wj[j];
                triplets[pos++] = Triplet(row, col, w);
            }
        }
    }
};

}  // namespace method
}  // namespace interpolation
}  // namespace atlas
