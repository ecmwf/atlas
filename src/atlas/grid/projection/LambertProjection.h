#pragma once

#include "atlas/grid/projection/Projection.h"

namespace atlas {
namespace grid {
namespace projection {

class LambertProjection: public Projection {

public:

    // constructor
    LambertProjection(const eckit::Parametrisation& p);

    // copy constructor
    LambertProjection( const LambertProjection& rhs );

    // clone method
    virtual Projection* clone() const ;

    // destructor
    ~LambertProjection() {}

    // class name
    static std::string className() { return "atlas.LambertProjection"; }
    static std::string projection_type_str() {return "lambert";}
    virtual std::string virtual_projection_type_str() const {return "lambert";}

    // projection and inverse projection
    virtual void xy2lonlat(double crd[]) const;
    virtual void lonlat2xy(double crd[]) const;

    virtual bool isStrictlyRegional() const { return true; }  // lambert projection cannot be used for global grids

    // specification
    virtual eckit::Properties spec() const;

    virtual std::string units() const { return "meters"; }

private:
  
    double lat1_, lat2_;     // First and second latitude at which the secant cone cuts the sphere
    bool is_tangent_;        // If the first and second latitude are equal, then the projection is on a tangent cone
    double lon0_;            // central longitude
    double radius_;          // sphere radius
    double n_, inv_n_, F_, rho0_, sign_;    // projection constants

    void setup();
};

}  // namespace projection
}  // namespace grid
}  // namespace atlas
