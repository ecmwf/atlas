#pragma once

#include "atlas/grid/detail/projection/Projection.h"
#include "atlas/grid/detail/projection/Rotation.h"

namespace atlas {
namespace grid {
namespace projection {

template <typename Rotation>
class SchmidtProjectionT: public Projection {

public:

    // constructor
    SchmidtProjectionT(const eckit::Parametrisation& p);
    SchmidtProjectionT() {}

    // class name
    static std::string static_type() { return Rotation::typePrefix()+"schmidt"; }
    virtual std::string type() const override { return static_type(); }

    // projection and inverse projection
    virtual void xy2lonlat(double crd[]) const override;
    virtual void lonlat2xy(double crd[]) const override;

    virtual bool strictlyRegional() const override { return false; }  // schmidt is global grid

    // specification
    virtual eckit::Properties spec() const override;

    virtual std::string units() const override { return "degrees"; }

private:

    double c_;    // stretching factor
    Rotation rotation_;

};

typedef SchmidtProjectionT<NotRotated> SchmidtProjection;
typedef SchmidtProjectionT<Rotated>    RotatedSchmidtProjection;

}  // namespace projection
}  // namespace grid
}  // namespace atlas
