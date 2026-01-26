#pragma once

#include "atlas/grid/detail/grid/Grid.h"
#include "atlas/runtime/Exception.h"
#include "atlas/util/Config.h"
#include "atlas/util/Point.h"

namespace atlas {
namespace grid {
namespace detail {
namespace grid {

class CubedSphere2 : public Grid {
private:

    // Get the lonlat and return as PointLonLat object
    struct ComputePointXY {
        ComputePointXY(const CubedSphere2& grid): grid_(grid) {}
        void operator()(idx_t n, PointXY& point) { grid_.xy(n, point); }
        const CubedSphere2& grid_;
    };

    // Get the lonlat and return as PointLonLat object
    struct ComputePointLonLat {
        ComputePointLonLat(const CubedSphere2& grid): grid_(grid) {}
        void operator()(idx_t n, PointLonLat& point) { grid_.lonlat(n, point); }
        const CubedSphere2& grid_;
    };

    template <typename Base, typename ComputePoint>
    class CubedSphere2Iterator : public Base {
    public:
        CubedSphere2Iterator(const CubedSphere2& grid, bool begin = true):
            grid_(grid),
            size_(grid.size()),
            n_(begin ? 0 : size_),
            compute_point{grid_}
            {
                if (n_ < size_) {
                    compute_point(n_, point_);
                }
            }

        // Return the point and move iterator to the next location
        bool next(typename Base::value_type& point) {
            if (n_ < size_) {
                compute_point(n_, point);
                ++n_;
                return true;
            }
            return false;
        }

        // * operator    
        const typename Base::reference operator*() const {return point_; }

        // ++ operator, move to next element and return point
        const Base& operator++() {
            ++n_;
            if (n_ < size_) {
                compute_point(n_, point_);
            }
            return *this;
        }

        // +=  operator, move some distance through the iterator and return point
        const Base& operator+=(typename Base::difference_type distance) {
            n_ += distance;
            if (n_ >= 0 && n_ < size_) {
                compute_point(n_, point_);
            }
            return *this;
        }

        // == operator, check two iterators are on the same index
        bool operator==(const Base& other) const {
            return n_ == static_cast<const CubedSphere2Iterator&>(other).n_;
        }

        // != operator, check two iterators are not on the same index
        bool operator!=(const Base& other) const {
            return n_ != static_cast<const CubedSphere2Iterator&>(other).n_;
        }

        // Return the number of points between two iterators
        typename Base::difference_type distance(const Base& other) const {
            return static_cast<const CubedSphere2Iterator&>(other).n_ - n_;
        }

        // Clone the iterator in its current position
        std::unique_ptr<Base> clone() const {
            auto result = new CubedSphere2Iterator(grid_, false);
            result->n_ = n_;
            result->compute_point(n_, result->point_);
            return std::unique_ptr<Base>(result);
        }

        const CubedSphere2& grid_;
        idx_t size_;
        idx_t n_;
        typename Base::value_type point_;
        ComputePoint compute_point;
    };

public:
    using IteratorXY = CubedSphere2Iterator<Grid::IteratorXY, ComputePointXY>;
    using IteratorLonLat = CubedSphere2Iterator<Grid::IteratorLonLat, ComputePointLonLat>;

    using Spec = atlas::util::Config;

    CubedSphere2(idx_t resolution);
    CubedSphere2(idx_t resolution, Projection projection);

    std::string name() const override;
    std::string type() const override;
    static std::string static_type();
    idx_t N() const {return N_;}

    void hash(eckit::Hash&) const override;
    RectangularLonLatDomain lonlatBoundingBox() const override;

    idx_t size() const override;
    Spec spec() const override;

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

    void xy(idx_t n, Point2& point) const;
    Point2 xy(idx_t n) const;

    void lonlat(idx_t n, Point2& point) const;
    Point2 lonlat(idx_t n) const;

protected:
    void print(std::ostream&) const override;

private:
    using CSIndices = std::array<idx_t, 3>;

    CSIndices get_cs_indices(gidx_t n) const;
    PointXY ij_to_tangent_coord(idx_t i, idx_t j) const;
    PointXYZ tangent_to_xyz_coord(const PointXY& tan_coord, idx_t tile) const;

protected:
    // (N_ * N_) = number of cells on a tile
    idx_t N_;

    // Number of tiles
    static constexpr idx_t nTiles_ = 6;

private:
    using Matrix = std::array<std::array<double, 3>, 3>;

    /*
    // currently unused
    std::array<Matrix, 6> lfric_rotations_ = {{
        {{ {0, 1, 0}, {0, 0, -1}, {1, 0, 0} }},
        {{ {-1, 0, 0}, {0, 0, -1}, {0, 1, 0} }},
        {{ {0, -1, 0}, {0, 0, -1}, {-1, 0, 0} }},
        {{ {1, 0, 0}, {0, 0, -1}, {0, -1, 0} }},
        {{ {-1, 0, 0}, {0, 1, 0}, {0, 0, 1} }},
        {{ {-1, 0, 0}, {0, -1, 0}, {0, 0, -1} }}
    }};
    */

    std::array<Matrix, 6> lfric_rotations_transposed_ = {{
        {{ {0, 0, 1}, {1, 0, 0}, {0, -1, 0} }},
        {{ {-1, 0, 0}, {0, 0, 1}, {0, -1, 0} }},
        {{ {0, 0, -1}, {-1, 0, 0}, {0, -1, 0} }},
        {{ {1, 0, 0}, {0, 0, -1}, {0, -1, 0} }},
        {{ {-1, 0, 0}, {0, 1, 0}, {0, 0, 1} }},
        {{ {-1, 0, 0}, {0, -1, 0}, {0, 0, -1} }}
    }};
};

}  // namespace grid
}  // namespace detail
}  // namespace grid
}  // namespace atlas
