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

#include <array>
#include <cmath>
#include <limits>

namespace atlas {
namespace interpolation {
namespace method {


class Cubic3DLimiter {
public:
    Cubic3DLimiter() = default;

    template <typename array_t, typename stencil_t>
    static typename std::enable_if<(array_t::RANK == 2), void>::type limit_scalar(
        typename std::remove_const<typename array_t::value_type>::type& output,
        const std::array<std::array<idx_t, 4>, 4>& index, const stencil_t& stencil, const array_t& input) {
        using Scalar = typename std::remove_const<typename array_t::value_type>::type;
        // Limit output to max/min of values in stencil marked by '*'
        //         x        x        x         x
        //              x     *-----*     x
        //                   /   P  |
        //          x       *------ *        x
        //        x        x        x         x
        idx_t k = stencil.k_interval();
        idx_t k1, k2;
        if (k < 1) {
            k1 = stencil.k(0);
            k2 = stencil.k(1);
        }
        else if (k > 1) {
            k1 = stencil.k(2);
            k2 = stencil.k(3);
        }
        else {
            k1 = stencil.k(1);
            k2 = stencil.k(2);
        }

        Scalar maxval = std::numeric_limits<Scalar>::lowest();
        Scalar minval = std::numeric_limits<Scalar>::max();
        for (idx_t j = 1; j < 3; ++j) {
            for (idx_t i = 1; i < 3; ++i) {
                idx_t n = index[j][i];

                Scalar f1 = input(n, k1);
                Scalar f2 = input(n, k2);

                maxval = std::max(maxval, f1);
                maxval = std::max(maxval, f2);
                minval = std::min(minval, f1);
                minval = std::min(minval, f2);
            }
        }
        if (output < minval) {
            output = minval;
        }
        else if (output > maxval) {
            output = maxval;
        }
    }

    template <typename InputArray, typename OutputArray, typename stencil_t>
    static typename std::enable_if<(InputArray::RANK == 3), void>::type limit_vars(
        const std::array<std::array<idx_t, 4>, 4>& index, const stencil_t& stencil, const InputArray& input,
        OutputArray& output, const idx_t nvar) {
        // Limit output to max/min of values in stencil marked by '*'
        //         x        x        x         x
        //              x     *-----*     x
        //                   /   P  |
        //          x       *------ *        x
        //        x        x        x         x

        using Value = typename OutputArray::value_type;

        const idx_t k = stencil.k_interval();
        idx_t k1, k2;
        if (k < 1) {
            k1 = stencil.k(0);
            k2 = stencil.k(1);
        }
        else if (k > 1) {
            k1 = stencil.k(2);
            k2 = stencil.k(3);
        }
        else {
            k1 = stencil.k(1);
            k2 = stencil.k(2);
        }

        for (idx_t v = 0; v < nvar; ++v) {
            Value limited = output[v];
            Value maxval  = std::numeric_limits<Value>::lowest();
            Value minval  = std::numeric_limits<Value>::max();
            for (idx_t j = 1; j < 3; ++j) {
                for (idx_t i = 1; i < 3; ++i) {
                    idx_t n = index[j][i];

                    Value f1 = input(n, k1, v);
                    Value f2 = input(n, k2, v);

                    maxval = std::max(maxval, f1);
                    maxval = std::max(maxval, f2);
                    minval = std::min(minval, f1);
                    minval = std::min(minval, f2);
                }
            }
            if (limited < minval) {
                limited = minval;
            }
            else if (limited > maxval) {
                limited = maxval;
            }
            output[v] = limited;
        }
    }
};

}  // namespace method
}  // namespace interpolation
}  // namespace atlas
