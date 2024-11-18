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

#include "eckit/linalg/Triplet.h"

#include "atlas/array/ArrayView.h"
#include "atlas/functionspace/StructuredColumns.h"
#include "atlas/grid/Stencil.h"
#include "atlas/grid/StencilComputer.h"
#include "atlas/runtime/Exception.h"
#include "atlas/util/CoordinateEnums.h"
#include "atlas/util/Point.h"

namespace atlas {
namespace interpolation {
namespace method {

class LinearHorizontalKernel {
    using Triplet  = eckit::linalg::Triplet;
    using Triplets = std::vector<Triplet>;

public:
    LinearHorizontalKernel() = default;

    LinearHorizontalKernel(const functionspace::StructuredColumns& fs, const util::Config& = util::NoConfig()) {
        src_ = fs;
        ATLAS_ASSERT(src_);
        compute_horizontal_stencil_ = grid::ComputeHorizontalStencil(src_.grid(), stencil_width());
    }

private:
    functionspace::StructuredColumns src_;
    grid::ComputeHorizontalStencil compute_horizontal_stencil_;

public:
    static std::string className() { return "LinearHorizontalKernel"; }
    static constexpr idx_t stencil_width() { return 2; }
    static constexpr idx_t stencil_size() { return stencil_width() * stencil_width(); }
    static constexpr idx_t stencil_halo() { return 0; }

public:
    using Stencil = grid::HorizontalStencil<2>;
    struct Weights {
        std::array<std::array<double, 2>, 2> weights_i;
        std::array<double, 2> weights_j;
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
    void make_valid_stencil(double& x, const double y, stencil_t& stencil, bool retry = true) const {
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
    void compute_weights(const double x, const double y, weights_t& weights) const {
        Stencil stencil;
        compute_stencil(x, y, stencil);
        compute_weights(x, y, stencil, weights);
    }


    template <typename stencil_t, typename weights_t>
    void compute_weights(const double x, const double y, const stencil_t& stencil, weights_t& weights) const {
        PointXY P1, P2;
        std::array<double, stencil_width()> yvec;
        // Compute weights for each constant-Y
        for (idx_t j = 0; j < stencil_width(); ++j) {
            auto& weights_i = weights.weights_i[j];
            src_.compute_xy(stencil.i(0, j), stencil.j(j), P1);
            src_.compute_xy(stencil.i(1, j), stencil.j(j), P2);
            const double alpha = (P2.x() - x) / (P2.x() - P1.x());
            weights_i[0]       = alpha;
            weights_i[1]       = 1. - alpha;
            yvec[j]            = P1.y();
        }
        // Compute weights in y-direction
        {
            auto& weights_j    = weights.weights_j;
            const double alpha = (yvec[1] - y) / (yvec[1] - yvec[0]);
            weights_j[0]       = alpha;
            weights_j[1]       = 1. - alpha;
        }
    }

    template <typename stencil_t, typename weights_t, typename array_t>
    typename array_t::value_type interpolate(const stencil_t& stencil, const weights_t& weights,
                                             const array_t& input) const {
        using Value = typename array_t::value_type;

        const auto& weights_j = weights.weights_j;
        Value output          = 0.;
        for (idx_t j = 0; j < stencil_width(); ++j) {
            const auto& weights_i = weights.weights_i[j];
            for (idx_t i = 0; i < stencil_width(); ++i) {
                idx_t n = src_.index(stencil.i(i, j), stencil.j(j));
                Value w = weights_i[i] * weights_j[j];
                output += w * input[n];
            }
        }

        return output;
    }

    template <typename stencil_t, typename weights_t, typename Value, int Rank>
    typename std::enable_if<(Rank == 1), void>::type interpolate(const stencil_t& stencil, const weights_t& weights,
                                                                 const array::ArrayView<const Value, Rank>& input,
                                                                 array::ArrayView<Value, Rank>& output, idx_t r) const {
        const auto& weights_j = weights.weights_j;
        output(r)             = 0.;
        for (idx_t j = 0; j < stencil_width(); ++j) {
            const auto& weights_i = weights.weights_i[j];
            for (idx_t i = 0; i < stencil_width(); ++i) {
                idx_t n = src_.index(stencil.i(i, j), stencil.j(j));
                Value w = static_cast<Value>(weights_i[i] * weights_j[j]);
                output(r) += w * input[n];
            }
        }
    }

    template <typename stencil_t, typename weights_t, typename Value, int Rank>
    typename std::enable_if<(Rank == 2), void>::type interpolate(const stencil_t& stencil, const weights_t& weights,
                                                                 const array::ArrayView<const Value, Rank>& input,
                                                                 array::ArrayView<Value, Rank>& output, idx_t r) const {
        const auto& weights_j = weights.weights_j;
        const idx_t Nk        = output.shape(1);
        for (idx_t k = 0; k < Nk; ++k) {
            output(r, k) = 0.;
        }
        for (idx_t j = 0; j < stencil_width(); ++j) {
            const auto& weights_i = weights.weights_i[j];
            for (idx_t i = 0; i < stencil_width(); ++i) {
                idx_t n = src_.index(stencil.i(i, j), stencil.j(j));
                Value w = static_cast<Value>(weights_i[i] * weights_j[j]);
                for (idx_t k = 0; k < Nk; ++k) {
                    output(r, k) += w * input(n, k);
                }
            }
        }
    }

    template <typename stencil_t, typename weights_t, typename Value, int Rank>
    typename std::enable_if<(Rank == 3), void>::type interpolate(const stencil_t& stencil, const weights_t& weights,
                                                                 const array::ArrayView<const Value, Rank>& input,
                                                                 array::ArrayView<Value, Rank>& output, idx_t r) const {
        const auto& weights_j = weights.weights_j;
        const idx_t Nk        = output.shape(1);
        const idx_t Nl        = output.shape(2);
        for (idx_t k = 0; k < Nk; ++k) {
            for (idx_t l = 0; l < Nl; ++l) {
                output(r, k, l) = 0.;
            }
        }
        for (idx_t j = 0; j < stencil_width(); ++j) {
            const auto& weights_i = weights.weights_i[j];
            for (idx_t i = 0; i < stencil_width(); ++i) {
                idx_t n = src_.index(stencil.i(i, j), stencil.j(j));
                Value w = static_cast<Value>(weights_i[i] * weights_j[j]);
                for (idx_t k = 0; k < Nk; ++k) {
                    for (idx_t l = 0; l < Nl; ++l) {
                        output(r, k, l) += w * input(n, k, l);
                    }
                }
            }
        }
    }

    template <typename array_t>
    typename array_t::value_type operator()(double x, double y, const array_t& input) const {
        Stencil stencil;
        Weights weights;
        compute_stencil(x, y, stencil);
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
        for (idx_t j = 0; j < stencil_width(); ++j) {
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
