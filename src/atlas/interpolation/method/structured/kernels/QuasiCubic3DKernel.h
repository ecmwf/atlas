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

#include "atlas/array/ArrayView.h"
#include "atlas/functionspace/StructuredColumns.h"
#include "atlas/grid/Stencil.h"
#include "atlas/grid/StencilComputer.h"
#include "atlas/runtime/Exception.h"
#include "atlas/util/CoordinateEnums.h"
#include "atlas/util/Point.h"

#include "Cubic3DLimiter.h"
#include "CubicVerticalKernel.h"
#include "LinearHorizontalKernel.h"
#include "QuasiCubicHorizontalKernel.h"

namespace atlas {
namespace interpolation {
namespace method {

struct QuasiCubicLinearPoints {
    constexpr QuasiCubicLinearPoints() {}
    static constexpr std::array<idx_t, 2> j{{1, 2}};
    static constexpr std::array<idx_t, 2> jj{{0, 3}};
    static constexpr std::array<idx_t, 2> jw{{4, 5}};
    static constexpr std::array<idx_t, 2> i{{1, 2}};
    static constexpr std::array<idx_t, 2> ii{{0, 3}};
};

class QuasiCubic3DKernel {
    using Limiter = Cubic3DLimiter;

public:
    QuasiCubic3DKernel(const functionspace::StructuredColumns& fs, const util::Config& config = util::NoConfig()) {
        src_ = fs;
        ATLAS_ASSERT(src_);
        ATLAS_ASSERT(src_.halo() >= 2);
        ATLAS_ASSERT(src_.vertical().size());
        quasi_cubic_horizontal_interpolation_ = QuasiCubicHorizontalKernel(src_, config);
        linear_horizontal_interpolation_      = LinearHorizontalKernel(src_, config);
        vertical_interpolation_               = CubicVerticalKernel(fs.vertical(), config);
        limiter_                              = config.getBool("limiter", false);
    }

private:
    functionspace::StructuredColumns src_;
    QuasiCubicHorizontalKernel quasi_cubic_horizontal_interpolation_;
    LinearHorizontalKernel linear_horizontal_interpolation_;
    CubicVerticalKernel vertical_interpolation_;
    bool limiter_{false};

public:
    static std::string className() { return "QuasiCubic3DKernel"; }
    static constexpr idx_t stencil_width() { return 4; }
    static constexpr idx_t stencil_size() { return 32; }
    static constexpr idx_t stencil_halo() {
        return static_cast<idx_t>(static_cast<double>(stencil_width()) / 2. + 0.5);
    }

public:
    using Stencil = grid::Stencil3D<4>;
    struct Weights {
        std::array<std::array<double, 4>, 4> weights_i;
        std::array<double, 6> weights_j;
        std::array<double, 4> weights_k;
        struct LinearWeights {
            std::array<std::array<double, 2>, 2> weights_i;
            std::array<double, 2> weights_j;
        } linear;
    };


public:
    struct WorkSpace {
        Stencil stencil;
        Weights weights;
    };

    template <typename stencil_t>
    void compute_stencil(double& x, const double y, const double z, stencil_t& stencil) const {
        quasi_cubic_horizontal_interpolation_.compute_stencil(x, y, stencil);
        vertical_interpolation_.compute_stencil(z, stencil);
    }

    template <typename weights_t>
    void compute_weights(double x, double y, double z, weights_t& weights) const {
        Stencil stencil;
        compute_stencil(x, y, z, stencil);
        compute_weights(x, y, z, stencil, weights);
    }


    template <typename stencil_t, typename weights_t>
    void compute_weights(const double x, const double y, const double z, const stencil_t& stencil,
                         weights_t& weights) const {
        quasi_cubic_horizontal_interpolation_.compute_weights(x, y, stencil, weights);


        // Insert more linear weights in available slots (weights_i[0], weights_i[3], weights_j[4], weights_j[5])
        {
            PointXY P1, P2;
            std::array<double, 2> yvec;
            constexpr QuasiCubicLinearPoints pts{};
            // Top and bottom row x-direction
            for (idx_t l = 0; l < 2; ++l) {
                idx_t j         = pts.j[l];   // index in stencil
                idx_t jj        = pts.jj[l];  // row index in weights_i
                auto& weights_i = weights.weights_i[jj];
                src_.compute_xy(stencil.i(pts.i[0], j), stencil.j(j), P1);
                src_.compute_xy(stencil.i(pts.i[1], j), stencil.j(j), P2);
                const double alpha   = (P2.x() - x) / (P2.x() - P1.x());
                weights_i[pts.ii[0]] = alpha;
                weights_i[pts.ii[1]] = 1. - alpha;
                yvec[l]              = P1.y();
            }
            // Compute weights in y-direction
            {
                auto& weights_j    = weights.weights_j;
                const double alpha = (yvec[1] - y) / (yvec[1] - yvec[0]);
                weights_j[4]       = alpha;
                weights_j[5]       = 1. - alpha;
            }
        }

        vertical_interpolation_.compute_weights(z, stencil, weights);
    }

    template <typename stencil_t, typename weights_t, typename array_t>
    typename std::enable_if<(array_t::RANK == 2), typename array_t::value_type>::type interpolate(
        const stencil_t& stencil, const weights_t& weights, const array_t& input) const {
        using Value = typename std::remove_const<typename array_t::value_type>::type;

        std::array<std::array<idx_t, stencil_width()>, stencil_width()> index;
        const auto& wk = weights.weights_k;

        Value output = 0.;

        // Horizontally quasi-cubic part for inner levels ( k = {1,2} )
        {
            // Inner levels, inner rows (cubic in i, cubic in j)   --> 16 points
            const auto& wj = weights.weights_j;
            for (idx_t j = 1; j < 3; ++j) {  // j = {1,2}
                const auto& wi = weights.weights_i[j];
                for (idx_t i = 0; i < 4; ++i) {  // i = {0,1,2,3}
                    idx_t n   = src_.index(stencil.i(i, j), stencil.j(j));
                    Value wij = wi[i] * wj[j];
                    for (idx_t k = 1; k < 3; ++k) {  // k = {1,2}
                        Value w = wij * wk[k];
                        output += w * input(n, stencil.k(k));
                    }
                    index[j][i] = n;
                }
            }
            // Inner levels, outer rows: (linear in i, cubic in j)  --> 8 points
            for (idx_t j = 0; j < 4; j += 3) {  // j = {0,3}
                const auto& wi = weights.weights_i[j];
                for (idx_t i = 1; i < 3; ++i) {  // i = {1,2}
                    idx_t n   = src_.index(stencil.i(i, j), stencil.j(j));
                    Value wij = wi[i] * wj[j];
                    for (idx_t k = 1; k < 3; ++k) {  // k = {1,2}
                        Value w = wij * wk[k];
                        output += w * input(n, stencil.k(k));
                    }
                    index[j][i] = n;
                }
            }
        }
        // Horizontally Linear part for outer levels ( k = {0,3} )
        {
            constexpr QuasiCubicLinearPoints pts{};
            // Outer levels: (linear in i, linear in j) -- > 8 points
            for (idx_t m = 0; m < 2; ++m) {
                idx_t j         = pts.j[m];   // index in stencil ( j = {1,2} )
                idx_t jj        = pts.jj[m];  // row index in weights_i ( jj = {0,3} )
                idx_t jw        = pts.jw[m];  // jw = {4,5};
                const auto& wi  = weights.weights_i[jj];
                const double wj = weights.weights_j[jw];
                for (idx_t l = 0; l < 2; ++l) {
                    idx_t i   = pts.i[l];   // i = {1,2}
                    idx_t ii  = pts.ii[l];  // ii = {0,3}
                    idx_t n   = src_.index(stencil.i(i, j), stencil.j(j));
                    Value wij = wi[ii] * wj;
                    for (idx_t k = 0; k < 4; k += 3) {  // k = {0,3}
                        Value w = wij * wk[k];
                        output += w * input(n, stencil.k(k));
                    }
                }
            }
        }


        if (limiter_) {
            Limiter::limit_scalar(output, index, stencil, input);
        }
        return output;
    }

    template <typename Value>
    struct OutputView1D {
        template <typename Int>
        Value& operator()(Int v) {
            return data_[v];
        }
        template <typename Int>
        Value& operator[](Int v) {
            return data_[v];
        }
        static constexpr int RANK{1};
        OutputView1D(Value* data): data_(data) {}
        using value_type = Value;

        Value* data_;
    };

    template <typename Value>
    OutputView1D<Value> make_outputview(Value* data) const {
        return OutputView1D<Value>(data);
    }

    template <typename stencil_t, typename weights_t, typename InputArray, typename OutputArray>
    typename std::enable_if<(InputArray::RANK == 3), void>::type interpolate_vars(const stencil_t& stencil,
                                                                                  const weights_t& weights,
                                                                                  const InputArray& input,
                                                                                  OutputArray& output,
                                                                                  const idx_t nvar) const {
        using Value = typename OutputArray::value_type;

        std::array<std::array<idx_t, stencil_width()>, stencil_width()> index;
        const auto& wk = weights.weights_k;

        const Value* _input_;

        for (idx_t v = 0; v < nvar; ++v) {
            output[v] = 0.;
        }

        // Horizontally quasi-cubic part for inner levels ( k = {1,2} )
        {
            // Inner levels, inner rows (cubic in i, cubic in j)   --> 16 points
            const auto& wj = weights.weights_j;
            for (idx_t j = 1; j < 3; ++j) {  // j = {1,2}
                const auto& wi = weights.weights_i[j];
                for (idx_t i = 0; i < 4; ++i) {  // i = {0,1,2,3}
                    idx_t n   = src_.index(stencil.i(i, j), stencil.j(j));
                    Value wij = wi[i] * wj[j];
                    for (idx_t k = 1; k < 3; ++k) {  // k = {1,2}
                        Value w        = wij * wk[k];
                        const idx_t kk = stencil.k(k);
                        _input_        = &(input(n, kk, 0));  // Assumption that input.stride(2) == 1
                        for (idx_t v = 0; v < nvar; ++v) {
                            output[v] += w * _input_[v];
                        }
                    }
                    index[j][i] = n;
                }
            }
            // Inner levels, outer rows: (linear in i, cubic in j)  --> 8 points
            for (idx_t j = 0; j < 4; j += 3) {  // j = {0,3}
                const auto& wi = weights.weights_i[j];
                for (idx_t i = 1; i < 3; ++i) {  // i = {1,2}
                    idx_t n   = src_.index(stencil.i(i, j), stencil.j(j));
                    Value wij = wi[i] * wj[j];
                    for (idx_t k = 1; k < 3; ++k) {  // k = {1,2}
                        Value w        = wij * wk[k];
                        const idx_t kk = stencil.k(k);
                        _input_        = &(input(n, kk, 0));  // Assumption that input.stride(2) == 1
                        for (idx_t v = 0; v < nvar; ++v) {
                            output[v] += w * _input_[v];
                        }
                    }
                    index[j][i] = n;
                }
            }
        }
        // Horizontally Linear part for outer levels ( k = {0,3} )
        {
            constexpr QuasiCubicLinearPoints pts{};
            // Outer levels: (linear in i, linear in j) -- > 8 points
            for (idx_t m = 0; m < 2; ++m) {
                idx_t j         = pts.j[m];   // index in stencil ( j = {1,2} )
                idx_t jj        = pts.jj[m];  // row index in weights_i ( jj = {0,3} )
                idx_t jw        = pts.jw[m];  // jw = {4,5};
                const auto& wi  = weights.weights_i[jj];
                const double wj = weights.weights_j[jw];
                for (idx_t l = 0; l < 2; ++l) {
                    idx_t i   = pts.i[l];   // i = {1,2}
                    idx_t ii  = pts.ii[l];  // ii = {0,3}
                    idx_t n   = src_.index(stencil.i(i, j), stencil.j(j));
                    Value wij = wi[ii] * wj;
                    for (idx_t k = 0; k < 4; k += 3) {  // k = {0,3}
                        Value w        = wij * wk[k];
                        const idx_t kk = stencil.k(k);
                        _input_        = &(input(n, kk, 0));  // Assumption that input.stride(2) == 1
                        for (idx_t v = 0; v < nvar; ++v) {
                            output[v] += w * _input_[v];
                        }
                    }
                }
            }
        }


        if (limiter_) {
            Limiter::limit_vars(index, stencil, input, output, nvar);
        }
    }


    template <typename stencil_t, typename weights_t, typename InputArray, typename OutputArray>
    typename std::enable_if<(InputArray::RANK == 2 && OutputArray::RANK == 1), void>::type interpolate(
        const stencil_t& stencil, const weights_t& weights, const InputArray& input, OutputArray& output,
        idx_t r) const {
        output(r) = interpolate(stencil, weights, input);
    }

    template <typename stencil_t, typename weights_t, typename InputArray, typename OutputArray>
    typename std::enable_if<(InputArray::RANK == 2 && OutputArray::RANK == 2), void>::type interpolate(
        const stencil_t& stencil, const weights_t& weights, const InputArray& input, OutputArray& output, idx_t r,
        idx_t k) const {
        output(r, k) = interpolate(stencil, weights, input);
    }

    template <typename stencil_t, typename weights_t, typename InputArray, typename OutputArray>
    typename std::enable_if<(InputArray::RANK == 3 && OutputArray::RANK == 3), void>::type interpolate(
        const stencil_t& stencil, const weights_t& weights, const InputArray& input, OutputArray& output, idx_t r,
        idx_t k) const {
        auto output_vars = make_outputview(&output(r, k, 0));
        interpolate_vars(stencil, weights, input, output_vars, output.shape(2));
    }

    template <typename stencil_t, typename weights_t, typename InputArray, typename OutputArray>
    typename std::enable_if<(InputArray::RANK == 2 && OutputArray::RANK == 3), void>::type interpolate(
        const stencil_t&, const weights_t&, const InputArray&, OutputArray&, idx_t /*r*/, idx_t /*k*/) const {
        ATLAS_NOTIMPLEMENTED;
    }

    template <typename stencil_t, typename weights_t, typename InputArray, typename OutputArray>
    typename std::enable_if<(InputArray::RANK == 3 && OutputArray::RANK == 1), void>::type interpolate(
        const stencil_t&, const weights_t&, const InputArray&, OutputArray&, idx_t /*r*/) const {
        ATLAS_NOTIMPLEMENTED;
    }

    template <typename stencil_t, typename weights_t, typename InputArray, typename OutputArray>
    typename std::enable_if<(InputArray::RANK == 3 && OutputArray::RANK == 1), void>::type interpolate(
        const stencil_t&, const weights_t&, const InputArray&, OutputArray&, idx_t /*r*/, idx_t /*k*/) const {
        ATLAS_NOTIMPLEMENTED;
    }

    template <typename stencil_t, typename weights_t, typename InputArray, typename OutputArray>
    typename std::enable_if<(InputArray::RANK == 3 && OutputArray::RANK == 2), void>::type interpolate(
        const stencil_t&, const weights_t&, const InputArray&, OutputArray&, idx_t /*r*/, idx_t /*k*/) const {
        ATLAS_NOTIMPLEMENTED;
    }
};  // namespace method

}  // namespace method
}  // namespace interpolation
}  // namespace atlas
