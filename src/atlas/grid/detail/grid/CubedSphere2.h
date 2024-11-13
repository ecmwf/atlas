#pragma once

#include "atlas/grid/detail/grid/Grid.h"
#include "atlas/runtime/Exception.h"
#include "atlas/util/Config.h"

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

protected:
    void print(std::ostream&) const override;

protected:
    // (N_ * N_) = number of cells on a tile
    idx_t N_;

    // Number of tiles
    static const idx_t nTiles_ = 6;

private:
    std::string type_ = {"cubedsphere2"};
};

}  // namespace grid
}  // namespace detail
}  // namespace grid
}  // namespace atlas
