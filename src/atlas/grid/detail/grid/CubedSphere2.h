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
public:
    using Spec = atlas::util::Config;

    CubedSphere2(idx_t resolution);

    std::string name() const override;
    std::string type() const override;

    void hash(eckit::Hash&) const override;
    RectangularLonLatDomain lonlatBoundingBox() const override;

    idx_t size() const override;
    Spec spec() const override;

    std::unique_ptr<IteratorXY> xy_begin() const override {
        ATLAS_NOTIMPLEMENTED;
    }
    std::unique_ptr<IteratorXY> xy_end() const override { ATLAS_NOTIMPLEMENTED; }
    std::unique_ptr<IteratorLonLat> lonlat_begin() const override {
        ATLAS_NOTIMPLEMENTED;
    }
    std::unique_ptr<IteratorLonLat> lonlat_end() const override {
        ATLAS_NOTIMPLEMENTED;
    }

    void lonlat(idx_t n, Point2& point) const;
    Point2 lonlat(idx_t n) const;

    void xy(idx_t n, Point2& point) const {lonlat(n, point);}
    Point2 xy(idx_t n) const {return lonlat(n);}

protected:
    void print(std::ostream&) const override;

private:
    int get_tile(idx_t n) const ;
    int get_tij(idx_t n) const ;
    int get_ti(idx_t n) const ;
    int get_tj(idx_t n) const ;
    double index_to_curvilinear(idx_t n) const ;

protected:
    // (N_ * N_) = number of cells on a tile
    idx_t N_;

    // Number of tiles
    static const idx_t nTiles_ = 6;

private:
    std::string type_ = {"cubedsphere2"};
    static constexpr int lfric_rotations_[6][9] = {
        {  0,  0,  1,  1,  0,  0,  0, -1,  0},
        { -1,  0,  0,  0,  0,  1,  0, -1,  0},
        {  0,  0, -1, -1,  0,  0,  0, -1,  0},
        {  1,  0,  0,  0,  0, -1,  0, -1,  0},
        { -1,  0,  0,  0,  1,  0,  0,  0,  1},
        { -1,  0,  0,  0, -1,  0,  0,  0, -1}
    };
    static constexpr double rad_to_deg_ = 180 / M_PI;
};

}  // namespace grid
}  // namespace detail
}  // namespace grid
}  // namespace atlas
