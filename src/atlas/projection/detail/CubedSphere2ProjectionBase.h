#pragma once

#include "atlas/domain/Domain.h"
#include "atlas/projection/Jacobian.h"
#include "atlas/projection/detail/ProjectionImpl.h"
#include "atlas/util/Config.h"

namespace atlas {
namespace projection {
namespace detail {

class CubedSphere2ProjectionBase final : public ProjectionImpl {
public:
    CubedSphere2ProjectionBase(const eckit::Parametrisation&);

    static std::string static_type() { return "cubedsphere2"; }
    std::string type() const override { return static_type(); }
    void xy2lonlat(double crd[]) const {}
    void lonlat2xy(double crd[]) const {}

    virtual Jacobian jacobian(const PointLonLat&) const override;
    bool strictlyRegional() const override { return false; }
    
    virtual RectangularLonLatDomain lonlatBoundingBox(const Domain&) const override {
        return GlobalDomain();
    }

    using Spec = atlas::util::Config;

    virtual Spec spec() const override;

    virtual std::string units() const override { return "degrees"; }

    virtual void hash(eckit::Hash&) const override;

};

}  // namespace detail
}  // namespace projection
}  // namespace atlas
