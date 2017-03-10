#pragma once

#include "atlas/grid/detail/projection/Projection.h"
#include "atlas/grid/detail/projection/Rotation.h"

#include "atlas/runtime/Log.h"

namespace atlas {
namespace grid {
namespace projection {

template <typename Rotation>
class LonLatProjectionT : public Projection {

private:
  friend class Projection;
  LonLatProjectionT(): Projection() {}

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
  static std::string static_type() { return Rotation::typePrefix()+"lonlat"; }
  virtual std::string type() const { return static_type(); }

  // projection and inverse projection
  virtual void xy2lonlat(double crd[]) const { rotation_.rotate(crd);   }
  virtual void lonlat2xy(double crd[]) const { rotation_.unrotate(crd); }

  virtual bool strictlyRegional() const { return false; }

  // specification
  virtual eckit::Properties spec() const;

  virtual std::string units() const { return "degrees"; }

  virtual operator bool() const override { return rotation_.rotated(); Log::info() << "rotated = " << rotation_.rotated();  }

private:

  Rotation rotation_;

};

typedef LonLatProjectionT<NotRotated> LonLatProjection;
typedef LonLatProjectionT<Rotated>    RotatedLonLatProjection;

}  // namespace projection
}  // namespace grid
}  // namespace atlas
