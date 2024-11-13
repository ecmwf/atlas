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

    std::string name() const override { ATLAS_NOTIMPLEMENTED; }
    std::string type() const override { ATLAS_NOTIMPLEMENTED; }

    void hash(eckit::Hash&) const override { ATLAS_NOTIMPLEMENTED; }
    RectangularLonLatDomain lonlatBoundingBox() const override {
        ATLAS_NOTIMPLEMENTED;
    }

    idx_t size() const override { ATLAS_NOTIMPLEMENTED; }
    Spec spec() const override { ATLAS_NOTIMPLEMENTED; }
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
    void print(std::ostream&) const override { ATLAS_NOTIMPLEMENTED; }
};

}  // namespace grid
}  // namespace detail
}  // namespace grid
}  // namespace atlas
