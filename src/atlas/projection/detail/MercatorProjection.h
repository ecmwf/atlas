#pragma once

#include "atlas/projection/detail/ProjectionImpl.h"

namespace atlas {
namespace projection {
namespace detail {

template <typename Rotation>
class MercatorProjectionT: public ProjectionImpl {

public:

    // constructor
    MercatorProjectionT(const eckit::Parametrisation& p);

    // class name
    static std::string static_type() { return Rotation::typePrefix()+"mercator"; }
    virtual std::string type() const override { return static_type(); }

    // projection and inverse projection
    virtual void xy2lonlat(double crd[]) const override;
    virtual void lonlat2xy(double crd[]) const override;
    
    virtual bool strictlyRegional() const override { return true; }  // Mercator projection cannot be used for global grids

    // specification
    virtual eckit::Properties spec() const override;

    virtual std::string units() const override { return "meters"; }
    
    virtual void hash( eckit::MD5& ) const override;

protected:

    double lon0_;            // central longitude
    double radius_;          // sphere radius
    double inv_radius_;      // 1/(sphere radius)

    void setup(const eckit::Parametrisation& p);

private:

  Rotation rotation_;

};

typedef MercatorProjectionT<NotRotated> MercatorProjection;
typedef MercatorProjectionT<Rotated>    RotatedMercatorProjection;


}  // namespace detail
}  // namespace projection
}  // namespace atlas
