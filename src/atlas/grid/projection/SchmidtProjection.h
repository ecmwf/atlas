#pragma once

#include "atlas/grid/projection/Projection.h"
#include "atlas/grid/projection/Rotation.h"

namespace atlas {
namespace grid {
namespace projection {

template <typename Rotation>
class SchmidtProjectionT: public Projection {

public:

    // constructor
    SchmidtProjectionT(const eckit::Parametrisation& p);
    SchmidtProjectionT() {}

    // copy constructor
    SchmidtProjectionT( const SchmidtProjectionT& rhs );

    // clone method
    virtual Projection* clone() const ;

    // class name
    static std::string className() { return "atlas."+Rotation::classNamePrefix()+"SchmidtProjection"; }
    static std::string projection_type_str() { return Rotation::typePrefix()+"schmidt"; }
    virtual std::string virtual_projection_type_str() const { return Rotation::typePrefix()+"schmidt"; }

    // projection and inverse projection
    virtual eckit::geometry::LLPoint2 coords2lonlat(eckit::geometry::Point2) const;
    virtual eckit::geometry::Point2 lonlat2coords(eckit::geometry::LLPoint2) const;

    // purely regional? - no!
    bool isRegional() { return false; }  // schmidt is global grid

    // specification
    virtual eckit::Properties spec() const;

private:

    double c_;    // stretching factor
    Rotation rotation_;

};

typedef SchmidtProjectionT<NotRotated> SchmidtProjection;
typedef SchmidtProjectionT<Rotated>    RotatedSchmidtProjection;

}  // namespace projection
}  // namespace grid
}  // namespace atlas
