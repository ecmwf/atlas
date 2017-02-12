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
    virtual void coords2lonlat(double crd[]) const;
    virtual void lonlat2coords(double crd[]) const;

    bool isRegional() { return true; }  // lambert projection cannot be used for global grids

    // specification
    virtual eckit::Properties spec() const;

private:
  
    double lat1_, lat2_;     // secant latitudes
    bool isTangent_;
    double lon0_;            // central longitude
    double radius_;          // sphere radius
    double n_, inv_n_, F_, rho0_;    // projection constants

    void setup();
};

}  // namespace projection
}  // namespace grid
}  // namespace atlas
