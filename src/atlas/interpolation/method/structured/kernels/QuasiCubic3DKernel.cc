/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction. and Interpolation
 */

#include "QuasiCubic3DKernel.h"

namespace atlas {
namespace interpolation {
namespace method {

// Note: Following symbols should no longer be necessary from C++17 onwards
constexpr std::array<idx_t, 2> QuasiCubicLinearPoints::j;
constexpr std::array<idx_t, 2> QuasiCubicLinearPoints::jj;
constexpr std::array<idx_t, 2> QuasiCubicLinearPoints::jw;
constexpr std::array<idx_t, 2> QuasiCubicLinearPoints::i;
constexpr std::array<idx_t, 2> QuasiCubicLinearPoints::ii;

}  // namespace method
}  // namespace interpolation
}  // namespace atlas
