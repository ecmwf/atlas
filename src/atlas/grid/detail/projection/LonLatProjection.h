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

  // destructor
  ~LonLatProjectionT() {}

  // class name
  static std::string static_type() { return Rotation::typePrefix()+"lonlat"; }
  virtual std::string type() const override { return static_type(); }

  // projection and inverse projection
  virtual void xy2lonlat(double crd[]) const override { rotation_.rotate(crd);   }
  virtual void lonlat2xy(double crd[]) const override { rotation_.unrotate(crd); }

  virtual bool strictlyRegional() const override { return false; }

  // specification
  virtual eckit::Properties spec() const override;

  virtual std::string units() const override { return "degrees"; }

  virtual operator bool() const override { return rotation_.rotated(); }

  virtual void hash( eckit::MD5& ) const override;

private:

  Rotation rotation_;

};

typedef LonLatProjectionT<NotRotated> LonLatProjection;
typedef LonLatProjectionT<Rotated>    RotatedLonLatProjection;

// --------------------------------------------------------------------------------------------------------------------

class ShiftedLonLatProjection : public Projection {

private:
  friend class Projection;
  ShiftedLonLatProjection(): Projection() {}

public:

  // constructor
  ShiftedLonLatProjection( double lon, double lat );
  ShiftedLonLatProjection( const eckit::Parametrisation& );

  // destructor
  ~ShiftedLonLatProjection() {}

  // class name
  static std::string static_type() { return "shifted_lonlat"; }
  virtual std::string type() const override { return static_type(); }

  // projection and inverse projection
  virtual void xy2lonlat(double crd[]) const override { crd[0] += lon_; crd[1] += lat_; }
  virtual void lonlat2xy(double crd[]) const override { crd[0] -= lon_; crd[1] -= lat_; }

  virtual bool strictlyRegional() const override { return false; }

  // specification
  virtual eckit::Properties spec() const override;

  virtual std::string units() const override { return "degrees"; }

  virtual operator bool() const override { return lon_!=0. || lat_!=0.; }

  virtual void hash( eckit::MD5& ) const override;

private:

  double lon_;
  double lat_;

};

// --------------------------------------------------------------------------------------------------------------------


}  // namespace projection
}  // namespace grid
}  // namespace atlas
