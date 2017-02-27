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
  virtual void xy2lonlat(double crd[]) const { rotation_.rotate(crd);   }
  virtual void lonlat2xy(double crd[]) const { rotation_.unrotate(crd); }

  virtual bool isStrictlyRegional() const { return false; }

  // specification
  virtual eckit::Properties spec() const;

  virtual std::string units() const { return "degrees"; }
  
  virtual operator bool() const { return rotation_.rotated(); }

private:

  Rotation rotation_;

};

typedef LonLatProjectionT<NotRotated> LonLatProjection;
typedef LonLatProjectionT<Rotated>    RotatedLonLatProjection;

}  // namespace projection
}  // namespace grid
}  // namespace atlas
