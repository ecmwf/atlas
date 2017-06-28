#pragma once

#include "atlas/projection/detail/ProjectionImpl.h"

namespace atlas {
namespace projection {
namespace detail {

template <typename Rotation>
class SchmidtProjectionT: public ProjectionImpl {

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
    virtual Spec spec() const override;

    virtual std::string units() const override { return "degrees"; }

    virtual void hash( eckit::Hash& ) const override;

private:

    double c_;    // stretching factor
    Rotation rotation_;

};

typedef SchmidtProjectionT<NotRotated> SchmidtProjection;
typedef SchmidtProjectionT<Rotated>    RotatedSchmidtProjection;

}  // namespace detail
}  // namespace projection
}  // namespace atlas
