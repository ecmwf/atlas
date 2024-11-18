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

class CubicHorizontalLimiter {
public:
    CubicHorizontalLimiter() = default;

    template <typename array_t>
    static void limit(typename array_t::value_type& output, const std::array<std::array<idx_t, 4>, 4>& index,
                      const array_t& input) {
        using Scalar = typename array_t::value_type;
        // Limit output to max/min of values in stencil marked by '*'
        //         x        x        x         x
        //              x     *-----*     x
        //                   /   P  |
        //          x       *------ *        x
        //        x        x        x         x
        Scalar maxval = std::numeric_limits<Scalar>::lowest();
        Scalar minval = std::numeric_limits<Scalar>::max();
        for (idx_t j = 1; j < 3; ++j) {
            for (idx_t i = 1; i < 3; ++i) {
                idx_t n    = index[j][i];
                Scalar val = input[n];
                maxval     = std::max(maxval, val);
                minval     = std::min(minval, val);
            }
        }
        if (output < minval) {
            output = minval;
        }
        else if (output > maxval) {
            output = maxval;
        }
    }

    template <typename Value, int Rank>
    static typename std::enable_if<(Rank == 1), void>::type limit(const std::array<std::array<idx_t, 4>, 4>& index,
                                                                  const array::ArrayView<const Value, Rank>& input,
                                                                  array::ArrayView<Value, Rank>& output, idx_t r) {
        // Limit output to max/min of values in stencil marked by '*'
        //         x        x        x         x
        //              x     *-----*     x
        //                   /   P  |
        //          x       *------ *        x
        //        x        x        x         x
        Value maxval = std::numeric_limits<Value>::lowest();
        Value minval = std::numeric_limits<Value>::max();
        for (idx_t j = 1; j < 3; ++j) {
            for (idx_t i = 1; i < 3; ++i) {
                idx_t n   = index[j][i];
                Value val = input[n];
                maxval    = std::max(maxval, val);
                minval    = std::min(minval, val);
            }
        }
        if (output(r) < minval) {
            output(r) = minval;
        }
        else if (output(r) > maxval) {
            output(r) = maxval;
        }
    }


    template <typename Value, int Rank>
    static typename std::enable_if<(Rank == 2), void>::type limit(const std::array<std::array<idx_t, 4>, 4>& index,
                                                                  const array::ArrayView<const Value, Rank>& input,
                                                                  array::ArrayView<Value, Rank>& output, idx_t r) {
        // Limit output to max/min of values in stencil marked by '*'
        //         x        x        x         x
        //              x     *-----*     x
        //                   /   P  |
        //          x       *------ *        x
        //        x        x        x         x
        for (idx_t k = 0; k < output.shape(1); ++k) {
            Value maxval = std::numeric_limits<Value>::lowest();
            Value minval = std::numeric_limits<Value>::max();
            for (idx_t j = 1; j < 3; ++j) {
                for (idx_t i = 1; i < 3; ++i) {
                    idx_t n   = index[j][i];
                    Value val = input(n, k);
                    maxval    = std::max(maxval, val);
                    minval    = std::min(minval, val);
                }
            }
            if (output(r, k) < minval) {
                output(r, k) = minval;
            }
            else if (output(r, k) > maxval) {
                output(r, k) = maxval;
            }
        }
    }

    template <typename Value, int Rank>
    static typename std::enable_if<(Rank == 3), void>::type limit(const std::array<std::array<idx_t, 4>, 4>& index,
                                                                  const array::ArrayView<const Value, Rank>& input,
                                                                  array::ArrayView<Value, Rank>& output, idx_t r) {
        // Limit output to max/min of values in stencil marked by '*'
        //         x        x        x         x
        //              x     *-----*     x
        //                   /   P  |
        //          x       *------ *        x
        //        x        x        x         x
        for (idx_t k = 0; k < output.shape(1); ++k) {
            for (idx_t l = 0; l < output.shape(2); ++l) {
                Value maxval = std::numeric_limits<Value>::lowest();
                Value minval = std::numeric_limits<Value>::max();
                for (idx_t j = 1; j < 3; ++j) {
                    for (idx_t i = 1; i < 3; ++i) {
                        idx_t n   = index[j][i];
                        Value val = input(n, k, l);
                        maxval    = std::max(maxval, val);
                        minval    = std::min(minval, val);
                    }
                }
                if (output(r, k, l) < minval) {
                    output(r, k, l) = minval;
                }
                else if (output(r, k, l) > maxval) {
                    output(r, k, l) = maxval;
                }
            }
        }
    }

};

}  // namespace method
}  // namespace interpolation
}  // namespace atlas
