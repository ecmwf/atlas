/*
 * (C) Copyright 2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <array>
#include <functional>
#include <iostream>
#include <memory>
#include <numeric>
#include <vector>

#include "eckit/types/Types.h"

#include "atlas/grid/Spacing.h"
#include "atlas/grid/Tiles.h"
#include "atlas/grid/detail/grid/Grid.h"
#include "atlas/library/config.h"
#include "atlas/projection/detail/CubedSphereProjectionBase.h"
#include "atlas/runtime/Exception.h"
#include "atlas/runtime/Log.h"
#include "atlas/util/CoordinateEnums.h"
#include "atlas/util/Object.h"
#include "atlas/util/ObjectHandle.h"
#include "atlas/util/Point.h"

namespace atlas {
namespace grid {
namespace detail {
namespace grid {


/**
 * @brief CubedSphere Grid
 *
 * This class is a base class for all grids that can be described as
 * a cubed sphere.
 *
 * For more detail on this implementation see atlas/grid/CubedSphereGrid.h
 */

using atlas::projection::detail::CubedSphereProjectionBase;

class CubedSphere : public Grid {
private:
    // Get the position in the xy plane and return as PointXY object
    struct ComputePointXY {
        ComputePointXY(const CubedSphere& grid): grid_(grid) {}
        void operator()(idx_t i, idx_t j, idx_t t, PointXY& point) { grid_.xy(i, j, t, point.data()); }
        const CubedSphere& grid_;
    };

    // Get the lonlat and return as PointLonLat object
    struct ComputePointLonLat {
        ComputePointLonLat(const CubedSphere& grid): grid_(grid) {}
        void operator()(idx_t i, idx_t j, idx_t t, PointLonLat& point) { grid_.lonlat(i, j, t, point.data()); }
        const CubedSphere& grid_;
    };

    class PointTIJ : public std::array<idx_t, 3> {
    public:
        using std::array<idx_t, 3>::array;
        idx_t t() const { return data()[0]; }
        idx_t i() const { return data()[1]; }
        idx_t j() const { return data()[2]; }
        idx_t& t() { return data()[0]; }
        idx_t& i() { return data()[1]; }
        idx_t& j() { return data()[2]; }
    };

    struct ComputePointTIJ {
        ComputePointTIJ(const CubedSphere& grid): grid_(grid) {}
        void operator()(idx_t i, idx_t j, idx_t t, PointTIJ& point) {
            point.t() = t;
            point.i() = i;
            point.j() = j;
        }
        const CubedSphere& grid_;
    };

    // -----------------------------------------------------------------------------------------------

    template <typename Base, typename ComputePoint>
    class CubedSphereIterator : public Base {
    public:
        // Create an iterator and return the first or last point. If begin is true it starts at the
        // beginning of the iterator, otherwise at the end. Class is templated and point can be xy or
        // lonlat.
        CubedSphereIterator(const CubedSphere& grid, bool begin = true):
            grid_(grid),
            i_(begin ? 0 : grid_.N()),
            j_(begin ? 0 : grid_.N()),
            t_(begin ? 0 : 5),
            size_(grid_.size()),
            n_(begin ? 0 : size_),
            compute_point{grid_} {
            // Check that point lies in grid and if so return the xy/lonlat
            if (grid_.inGrid(i_, j_, t_)) {
                compute_point(i_, j_, t_, point_);
            }
        }

        // Return the point and move iterator to the next location
        virtual bool next(typename Base::value_type& point) {
            if (n_ != size_) {
                compute_point(i_, j_, t_, point);
                std::unique_ptr<int[]> ijt = grid_.nextElement(i_, j_, t_);
                i_                         = ijt[0];
                j_                         = ijt[1];
                t_                         = ijt[2];
                ++n_;
                return true;
            }
            return false;
        }

        // * operator
        virtual const typename Base::reference operator*() const { return point_; }

        // ++ operator, move to next element in grid iterator and return point
        virtual const Base& operator++() {
            std::unique_ptr<int[]> ijt = grid_.nextElement(i_, j_, t_);
            i_                         = ijt[0];
            j_                         = ijt[1];
            t_                         = ijt[2];
            ++n_;
            if (n_ != size_) {
                compute_point(i_, j_, t_, point_);
            }
            return *this;
        }

        // += operator, move some distance d through the iterator and return point
        virtual const Base& operator+=(typename Base::difference_type distance) {
            idx_t d = distance;
            // Following loop could be optimised to not iterate through every point,
            // but rather jump through a tile at a time if possible.
            // Then OpenMP algorithms can be made much quicker.
            for (int n = 0; n < d; n++) {
                std::unique_ptr<int[]> ijt = grid_.nextElement(i_, j_, t_);
                i_                         = ijt[0];
                j_                         = ijt[1];
                t_                         = ijt[2];
            }
            n_ += d;
            if (n_ != size_) {
                compute_point(i_, j_, t_, point_);
            }
            return *this;
        }

        // Given two positions in the grid iterator return the distance, which for the cubed-sphere
        // grid is just the number of grid points between the two points.
        virtual typename Base::difference_type distance(const Base& other) const {
            const auto& _other               = static_cast<const CubedSphereIterator&>(other);
            typename Base::difference_type d = 0;
            idx_t i                          = i_;
            idx_t j                          = j_;
            idx_t t                          = t_;
            bool found                       = false;
            for (int n = 0; n < grid_.size(); n++) {
                if (i == _other.i_ && j == _other.j_ && t == _other.t_) {
                    found = true;
                    break;
                }
                std::unique_ptr<int[]> ijt = grid_.nextElement(i, j, t);
                i                          = ijt[0];
                j                          = ijt[1];
                t                          = ijt[2];
                ++d;
            }
            ATLAS_ASSERT(!found, "CubedSphereIterator.distance: cycled entire grid without finding other");
            return d;
        }

        // == operator for checking two positions in the iterator are equal
        virtual bool operator==(const Base& other) const {
            return i_ == static_cast<const CubedSphereIterator&>(other).i_ &&
                   j_ == static_cast<const CubedSphereIterator&>(other).j_ &&
                   t_ == static_cast<const CubedSphereIterator&>(other).t_;
        }

        // != operator for checking that two positions in the iterator are not equal
        virtual bool operator!=(const Base& other) const {
            return i_ != static_cast<const CubedSphereIterator&>(other).i_ ||
                   j_ != static_cast<const CubedSphereIterator&>(other).j_ ||
                   t_ != static_cast<const CubedSphereIterator&>(other).t_;
        }

        // Clone the grid iterator
        virtual std::unique_ptr<Base> clone() const {
            auto result    = new CubedSphereIterator(grid_, false);
            result->i_     = i_;
            result->j_     = j_;
            result->t_     = t_;
            result->point_ = point_;
            result->size_  = size_;
            result->n_     = n_;
            return std::unique_ptr<Base>(result);
        }

        const CubedSphere& grid_;
        idx_t i_;
        idx_t j_;
        idx_t t_;
        idx_t size_;
        idx_t n_;
        typename Base::value_type point_;
        ComputePoint compute_point;
    };

    // -----------------------------------------------------------------------------------------------

public:
    // Iterators for returning xy or lonlat
    using IteratorXY     = CubedSphereIterator<Grid::IteratorXY, ComputePointXY>;
    using IteratorLonLat = CubedSphereIterator<Grid::IteratorLonLat, ComputePointLonLat>;

    class IteratorTIJ_Base : public IteratorT<IteratorTIJ_Base, PointTIJ> {};
    using IteratorTIJ = CubedSphereIterator<IteratorTIJ_Base, ComputePointTIJ>;

    static std::string static_type();

    // Constructors
    CubedSphere(const std::string&, int, Projection, const std::string& stagger);
    CubedSphere(int, Projection, const std::string& stagger);
    CubedSphere(const CubedSphere&);

    // Destructor
    virtual ~CubedSphere() override;

    // Return total grid size
    virtual idx_t size() const override { return accumulate(npts_.begin(), npts_.end(), 0); }

    // Return information about the grid
    virtual Spec spec() const override;
    virtual std::string name() const override;
    virtual std::string type() const override;

    // Return N_, where (N_ * N_) is the number of cells on a tile
    inline idx_t N() const { return N_; }

    // Access to the tile class
    inline atlas::grid::CubedSphereTiles tiles() const { return tiles_; }

    // Tile specific access to x and y locations
    // -----------------------------------------

    inline double xsPlusIndex(idx_t idx, idx_t t) const {
        return static_cast<double>(xs_[t]) + static_cast<double>(idx);
    }

    inline double xsrMinusIndex(idx_t idx, idx_t t) const {
        return static_cast<double>(xsr_[t]) - static_cast<double>(idx);
    }

    inline double ysPlusIndex(idx_t idx, idx_t t) const {
        return static_cast<double>(ys_[t]) + static_cast<double>(idx);
    }

    inline double ysrMinusIndex(idx_t idx, idx_t t) const {
        return static_cast<double>(ysr_[t]) - static_cast<double>(idx);
    }

    // Lambdas for access to appropriate functions for tile
    // ----------------------------------------------------

    std::vector<std::function<double(int, int, int)>> xtile;
    std::vector<std::function<double(int, int, int)>> ytile;

    // Functions for returning xy
    // --------------------------

    inline void xyt(idx_t i, idx_t j, idx_t t, double crd[]) const {
        crd[0] = xtile.at(t)(i, j, t);
        crd[1] = ytile.at(t)(i, j, t);
        crd[2] = static_cast<double>(t);
    }

    PointXY xyt(idx_t i, idx_t j, idx_t t) const { return PointXY(xtile.at(t)(i, j, t), ytile.at(t)(i, j, t)); }

    inline void xy(idx_t i, idx_t j, idx_t t, double xy[]) const {
        double crd[3];
        this->xyt(i, j, t, crd);
        this->xyt2xy(crd, xy);
    }

    PointXY xy(idx_t i, idx_t j, idx_t t) const {
        double crd[2];
        this->xy(i, j, t, crd);
        return PointXY(crd[0], crd[1]);
    }

    // Functions for returning lonlat, either as array or PointLonLat
    // --------------------------------------------------------------

    void lonlat(idx_t i, idx_t j, idx_t t, double lonlat[]) const {
        this->xy(i, j, t, lonlat);      // outputing xy in lonlat
        projection_.xy2lonlat(lonlat);  // converting xy to lonlat
    }

    PointLonLat lonlat(idx_t i, idx_t j, idx_t t) const {
        double lonlat[2];
        this->lonlat(i, j, t, lonlat);
        return PointLonLat(lonlat[LON], lonlat[LAT]);
    }

    // Check whether i, j, t is in grid
    // --------------------------------
    inline bool inGrid(idx_t i, idx_t j, idx_t t) const {
        constexpr idx_t tmax = 5;
        if (t >= 0 && t <= tmax) {
            if (j >= jmin_[t] && j <= jmax_[t]) {
                if (i >= imin_[t][j] && i <= imax_[t][j]) {
                    return true;
                }
            }
        }
        return false;
    }

    // Check on whether the final element
    // ----------------------------------

    bool finalElement(idx_t i, idx_t j, idx_t t) const {
        constexpr idx_t tmax = 5;
        if (t == tmax) {
            idx_t jmax = jmax_[tmax];
            if (j == jmax) {
                if (i == imax_[tmax][jmax]) {
                    return true;
                }
            }
        }
        return false;
    }

    // Move to next grid element in an iterator
    // ----------------------------------------

    // Note that i is the fastest index, followed by j, followed by t
    std::unique_ptr<int[]> nextElement(const idx_t i, const idx_t j, const idx_t t) const {
        auto ijt = std::make_unique<int[]>(3);

        ijt[0] = i;
        ijt[1] = j;
        ijt[2] = t;

        if (i < imax_[t][j]) {
            ijt[0] = i + 1;
            ijt[1] = j;
            ijt[2] = t;
            return ijt;
        }

        if (i == imax_[t][j]) {
            if (j < jmax_[t]) {
                // move to next column
                ijt[0] = 0;
                ijt[1] = j + 1;
                ijt[2] = t;
                return ijt;
            }

            if (j == jmax_[t]) {
                if (t < nTiles_ - 1) {
                    // move to next tile
                    ijt[0] = 0;
                    ijt[1] = 0;
                    ijt[2] = t + 1;
                    return ijt;
                }

                if (t == nTiles_ - 1) {
                    // We are at the final point so we go to
                    // to a point that defines the "end()" of the
                    // iterator i.e. it is not a point on the grid
                    // For now it is set at (N_, N_, nTiles -1)
                    ijt[0] = N_;
                    ijt[1] = N_;
                    ijt[2] = nTiles_ - 1;
                    return ijt;
                }
            }
        }

        return ijt;
    }

    // Iterator start/end positions
    // ----------------------------

    virtual std::unique_ptr<Grid::IteratorXY> xy_begin() const override {
        return std::make_unique<IteratorXY>(*this);
    }
    virtual std::unique_ptr<Grid::IteratorXY> xy_end() const override {
        return std::make_unique<IteratorXY>(*this, false);
    }
    virtual std::unique_ptr<Grid::IteratorLonLat> lonlat_begin() const override {
        return std::make_unique<IteratorLonLat>(*this);
    }
    virtual std::unique_ptr<Grid::IteratorLonLat> lonlat_end() const override {
        return std::make_unique<IteratorLonLat>(*this, false);
    }
    virtual std::unique_ptr<IteratorTIJ> tij_begin() const {
        return std::make_unique<IteratorTIJ>(*this);
    }
    virtual std::unique_ptr<IteratorTIJ> tij_end() const {
        return std::make_unique<IteratorTIJ>(*this, false);
    }

    // Default configurations

    Config meshgenerator() const override;

    Config partitioner() const override;

    const std::string& stagger() const { return stagger_; }

protected:
    virtual void print(std::ostream&) const override;

    virtual void hash(eckit::Hash&) const override;

    virtual RectangularLonLatDomain lonlatBoundingBox() const override;

    Domain computeDomain() const;

private:
    void xy2xyt(const double xy[], double xyt[]) const;  // note: unused!

    void xyt2xy(const double xyt[], double xy[]) const;

protected:
    // (N_ * N_) = number of cells on a tile
    idx_t N_;

    // Number of tiles
    static const idx_t nTiles_ = 6;

    // Start points in x,y direction
    double xs_[nTiles_];
    double ys_[nTiles_];
    double xsr_[nTiles_];  // x order reversed
    double ysr_[nTiles_];  // y order reversed (for FV3 panels 4, 5, 6)

    // Number of unique points on each tile
    std::vector<int> npts_;

    std::string tileType_;

    std::array<idx_t, 6> jmin_;
    std::array<idx_t, 6> jmax_;
    std::vector<std::vector<idx_t>> imin_;
    std::vector<std::vector<idx_t>> imax_;

    std::string stagger_;

private:
    std::string name_ = {"cubedsphere"};
    CubedSphereProjectionBase* cs_projection_;  // store pointer to dynamic_cast for convenience
    atlas::grid::CubedSphereTiles tiles_;
    std::array<std::array<double, 6>, 2> tiles_offsets_xy2ab_;
    std::array<std::array<double, 6>, 2> tiles_offsets_ab2xy_;
};


}  // namespace grid
}  // namespace detail
}  // namespace grid
}  // namespace atlas
