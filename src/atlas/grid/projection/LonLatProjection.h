#pragma once

#include "atlas/grid/projection/Projection.h"
#include "atlas/grid/projection/Rotation.h"

namespace atlas {
namespace grid {
namespace projection {

template <typename Rotation>
class LonLatProjectionT : public Projection {

public:

  // constructor
  LonLatProjectionT( const eckit::Parametrisation& );

  // copy constructor
  LonLatProjectionT( const LonLatProjectionT& );

  // clone method
  virtual Projection* clone() const ;

  // destructor
  ~LonLatProjectionT() {}

  // class name
  static std::string className() { return "atlas."+Rotation::classNamePrefix()+"LonLatProjection"; }
  static std::string projection_type_str() { return Rotation::typePrefix()+"lonlat"; }
  virtual std::string virtual_projection_type_str() const { return Rotation::typePrefix()+"lonlat"; }

  // projection and inverse projection
  virtual eckit::geometry::LLPoint2 coords2lonlat(eckit::geometry::Point2 xy) const {
    eckit::geometry::LLPoint2 lonlat(xy[0],xy[1]);
    rotation_.rotate(lonlat);
    return lonlat;
  }

  virtual eckit::geometry::Point2 lonlat2coords(eckit::geometry::LLPoint2 ll) const {
    eckit::geometry::LLPoint2 coords(ll.lon(),ll.lat());
    rotation_.unrotate(coords);
    return coords;
  }

  // purely regional? - no!
  bool isRegional() { return false; }

  // specification
  virtual eckit::Properties spec() const;

private:

  Rotation rotation_;

};

typedef LonLatProjectionT<NotRotated> LonLatProjection;
typedef LonLatProjectionT<Rotated>    RotatedLonLatProjection;

}  // namespace projection
}  // namespace grid
}  // namespace atlas
