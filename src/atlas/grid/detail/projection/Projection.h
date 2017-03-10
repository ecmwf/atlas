#pragma once

#include <array>
#include "atlas/util/Point.h"
#include "eckit/config/Parametrisation.h"
#include "eckit/value/Properties.h"
#include "eckit/memory/Builder.h"
#include "eckit/memory/Owned.h"

namespace atlas {
namespace grid {
namespace projection {

class Projection : public eckit::Owned {

public:

    using ARG1       = const eckit::Parametrisation&;
    using builder_t  = eckit::BuilderT1<Projection>;
    static std::string className() {return "atlas.Projection";}

public:

    static Projection* create(); // creates the LonLatProjection
    static Projection* create(const eckit::Parametrisation& p);

    Projection() {}
    Projection( const Projection& rhs ) {}   // copy constructor
    virtual Projection* clone() const =0;    // clone method acting like virtual copy constructor
    virtual ~Projection() {}                 // destructor should be virtual when using a virtual copy constructor

    virtual std::string type() const =0;

    virtual void xy2lonlat(double crd[]) const =0;
    virtual void lonlat2xy(double crd[]) const =0;

    PointLonLat lonlat( const PointXY& ) const;
    PointXY xy( const PointLonLat& ) const;

    virtual bool strictlyRegional() const =0;

    virtual eckit::Properties spec() const =0;

    virtual std::string units() const =0;

    virtual operator bool() const { return true; }

};

inline PointLonLat Projection::lonlat( const PointXY& xy ) const {
  PointLonLat lonlat(xy);
  xy2lonlat(lonlat.data());
  return lonlat;
}

inline PointXY Projection::xy( const PointLonLat& lonlat ) const {
  PointXY xy(lonlat);
  lonlat2xy(xy.data());
  return xy;
}

}  // namespace projection
}  // namespace grid
}  // namespace atlas
