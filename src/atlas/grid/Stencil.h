/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#pragma once

#include <array>
#include "atlas/library/config.h"

namespace atlas {
namespace grid {

//---------------------------------------------------------------------------------------------------------------------

template <idx_t StencilWidth>
class HorizontalStencil {
    friend class ComputeHorizontalStencil;
    std::array<idx_t, StencilWidth> i_begin_;
    idx_t j_begin_;

public:
    idx_t i(idx_t offset_i, idx_t offset_j) const { return i_begin_[offset_j] + offset_i; }
    idx_t j(idx_t offset) const { return j_begin_ + offset; }
    constexpr idx_t width() const { return StencilWidth; }
};

//-----------------------------------------------------------------------------

template <idx_t StencilWidth>
class VerticalStencil {
    friend class ComputeVerticalStencil;
    idx_t k_begin_;
    idx_t k_interval_;
    // within the stencil, a point falls in a certain interval.
    // e.g. for a cubic stencil:
    // +  |-----|-----|-----|      --> k_interval_ = -1
    //    |--+--|-----|-----|      --> k_interval_ = 0
    //    |-----|--+--|-----|      --> k_interval_ = 1 (the centred case)
    //    |-----|-----|--+--|      --> k_interval_ = 2
    //    |-----|-----|-----|  +   --> k_interval_ = 3
public:
    idx_t k(idx_t offset) const { return k_begin_ + offset; }
    constexpr idx_t width() const { return StencilWidth; }
    idx_t k_interval() const { return k_interval_; }
};

//-----------------------------------------------------------------------------

template <idx_t StencilWidth>
class Stencil3D {
    friend class ComputeHorizontalStencil;
    friend class ComputeVerticalStencil;
    std::array<idx_t, StencilWidth> i_begin_;
    idx_t j_begin_;
    idx_t k_begin_;
    idx_t k_interval_;

public:
    idx_t i(idx_t offset_i, idx_t offset_j) const { return i_begin_[offset_j] + offset_i; }
    idx_t j(idx_t offset) const { return j_begin_ + offset; }
    idx_t k(idx_t offset) const { return k_begin_ + offset; }
    constexpr idx_t width() const { return StencilWidth; }
    idx_t k_interval() const { return k_interval_; }
};

//---------------------------------------------------------------------------------------------------------------------

}  // namespace grid
}  // namespace atlas
