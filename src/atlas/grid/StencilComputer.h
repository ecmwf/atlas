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

#include <cmath>
#include <vector>

#include "atlas/grid/Vertical.h"
#include "atlas/library/config.h"
#include "atlas/runtime/Exception.h"

namespace atlas {
class StructuredGrid;
}  // namespace atlas

namespace atlas {
namespace grid {

class ComputeLower {
    std::vector<double> z_;
    std::vector<idx_t> nvaux_;
    idx_t nlev_;
    idx_t nlevaux_;
    double rlevaux_;

public:
    ComputeLower() = default;

    ComputeLower(const Vertical& z);

    idx_t operator()(double z) const {
        idx_t idx = static_cast<idx_t>(std::floor(z * rlevaux_));
#ifndef NDEBUG
        ATLAS_ASSERT(idx < static_cast<idx_t>(nvaux_.size()) && idx >= 0);
#endif
        idx = nvaux_[idx];
        if (idx < nlev_ - 1 && z > z_[idx + 1]) {
            ++idx;
        }
        return idx;
    }
};

//-----------------------------------------------------------------------------

class ComputeNorth {
    std::vector<double> y_;
    double dy_;
    idx_t halo_;
    idx_t ny_;
    static constexpr double tol() { return 0.5e-6; }

public:
    ComputeNorth() = default;

    ComputeNorth(const StructuredGrid& grid, idx_t halo);

    idx_t operator()(double y) const {
        idx_t j = static_cast<idx_t>(std::floor((y_[halo_ + 0] - y) / dy_));
        j       = std::max<idx_t>(halo_, std::min<idx_t>(j, halo_ + ny_ - 1));
        while (y_[halo_ + j] > y) {
            ++j;
        }
        do {
            --j;
        } while (y_[halo_ + j] < y);

        return j;
    }
};

//-----------------------------------------------------------------------------

class ComputeWest {
    std::vector<double> dx;
    std::vector<double> xref;
    idx_t halo_;  // halo in north-south direction
    idx_t ny_;
    static constexpr double tol() { return 0.5e-6; }

public:
    ComputeWest() = default;

    ComputeWest(const StructuredGrid& grid, idx_t halo = 0);

    idx_t operator()(const double& x, idx_t j) const {
        idx_t jj = halo_ + j;
        idx_t i  = static_cast<idx_t>(std::floor((x - xref[jj]) / dx[jj]));
        return i;
    }
};


//-----------------------------------------------------------------------------

/// @class ComputeHorizontalStencil
/// @brief Compute stencil in horizontal direction (i,j)
///
/// @details
/// Given a stencil width, the stencil for a given P{x,y} is:
/// @code
///        i[0]     i[1]     i[2]    i[3]
///         x        x        x         x       j + 0
///          x       x       x        x         j + 1
///                     P
///          x       x       x        x         j + 2
///         x        x        x         x       j + 3
/// @endcode
///
/// In case the x-component of P is aligned with any
/// stencil, gridpoint, the stencil will assume the grid point
/// is on the point P's left side:
/// @code
///        i[0]     i[1]     i[2]    i[3]
///         x        x        x         x       j + 0
///          x       x       x        x         j + 1
///                  P
///          x       x       x        x         j + 2
///         x        x        x         x       j + 3
/// @endcode
class ComputeHorizontalStencil {
    idx_t halo_;
    ComputeNorth compute_north_;
    ComputeWest compute_west_;
    idx_t stencil_width_;
    idx_t stencil_begin_;

public:
    ComputeHorizontalStencil() = default;

    ComputeHorizontalStencil(const StructuredGrid& grid, idx_t stencil_width);

    template <typename stencil_t>
    void operator()(const double& x, const double& y, stencil_t& stencil) const {
        stencil.j_begin_ = compute_north_(y) - stencil_begin_;
        for (idx_t jj = 0; jj < stencil_width_; ++jj) {
            stencil.i_begin_[jj] = compute_west_(x, stencil.j_begin_ + jj) - stencil_begin_;
        }
    }
};


//-----------------------------------------------------------------------------

//---------------------------------------------------------------------------------------------------------------------
/// @class ComputeVerticalStencil
/// @brief Compute lower vertical level index for given coordinate
/// zcoord:
/// @verbatim
///   0----1----2----3--...--(n-1)----(n)----(n+1)
///   --->|<---|<---|<--...-|<--------------------
/// @endverbatim
/// If coordinate falls on vertical level (+- epsilon), that level is returned
/// If coordinate falls in range [0,1) or [n,n+1],
/// the index is snapped to 1 and (n-1) respectively. This allows reliably that
/// the returned index can be used for stencil operations.
///
/// IFS full levels don't have a level at the boundaries (0.,1.)
/// It is the "half" levels that contain (0.,1.). For reasons of boundary conditions
/// however, the full levels also have 0. prepended and 1. appended.
///
/// Example IFS full levels for regular distribution dz ( level 0 and n+1 are added for boundary conditions )
///  0      :  0.0
///  jlev   :  jlev*dz - 0.5*dz
///  nlev   :  nlev*dz - 0.5*dz
///  nlev+1 :  1.0
class ComputeVerticalStencil {
    ComputeLower compute_lower_;
    idx_t stencil_width_;
    idx_t stencil_begin_;

    idx_t clip_begin_;
    idx_t clip_end_;
    double vertical_min_;
    double vertical_max_;

public:
    ComputeVerticalStencil() = default;

    ComputeVerticalStencil(const Vertical& vertical, idx_t stencil_width);

    template <typename stencil_t>
    void operator()(const double& z, stencil_t& stencil) const {
        idx_t k_begin = compute_lower_(z) - stencil_begin_;
        idx_t k_end   = k_begin + stencil_width_;
        idx_t move    = 0;

        if (k_begin < clip_begin_) {
            move = k_begin - clip_begin_;
            if (z < vertical_min_) {
                --k_begin;
                --move;
            }
        }
        else if (k_end > clip_end_) {
            move = k_end - clip_end_;
        }
        stencil.k_begin_    = k_begin - move;
        stencil.k_interval_ = stencil_begin_ + move;
    }
};

//---------------------------------------------------------------------------------------------------------------------

}  // namespace grid
}  // namespace atlas
