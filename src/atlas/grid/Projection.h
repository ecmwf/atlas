#pragma once

#include <string>
#include "eckit/config/Parametrisation.h"
#include "eckit/value/Properties.h"
#include "eckit/memory/SharedPtr.h"

#include "atlas/util/Point.h"
#include "atlas/grid/detail/projection/Projection.h"

//---------------------------------------------------------------------------------------------------------------------

// Forward declarations
namespace eckit {
class Parametrisation;
class MD5;
}

//---------------------------------------------------------------------------------------------------------------------

namespace atlas {
namespace grid {

//---------------------------------------------------------------------------------------------------------------------

class Projection {

public:

  Projection();
  Projection( const Projection& );
  Projection( const atlas::grid::projection::Projection* );
  Projection( const eckit::Parametrisation& );

  operator atlas::grid::projection::Projection*() { return projection_.get(); }
  operator const atlas::grid::projection::Projection&() const { return *projection_.get(); }

  void xy2lonlat(double crd[]) const;
  void lonlat2xy(double crd[]) const;

  PointLonLat lonlat( const PointXY& ) const;
  PointXY xy( const PointLonLat& ) const;

  bool strictlyRegional() const;

  eckit::Properties spec() const;

  std::string units() const;

  operator bool() const;

  std::string type() const { return projection_->type(); }
  
  void hash( eckit::MD5& ) const;

private:

  eckit::SharedPtr<atlas::grid::projection::Projection> projection_;
};

//---------------------------------------------------------------------------------------------------------------------

class ShiftedLonLatProjection : public Projection {

public:

  ShiftedLonLatProjection( double lon, double lat );

};

//---------------------------------------------------------------------------------------------------------------------

inline void Projection::xy2lonlat(double crd[]) const { return projection_->xy2lonlat(crd); }
inline void Projection::lonlat2xy(double crd[]) const { return projection_->lonlat2xy(crd); }
inline PointLonLat Projection::lonlat( const PointXY& xy) const { return projection_->lonlat(xy); }
inline PointXY Projection::xy( const PointLonLat& lonlat) const { return projection_->xy(lonlat); }
inline bool Projection::strictlyRegional() const { return projection_->strictlyRegional(); }
inline eckit::Properties Projection::spec() const { return projection_->spec(); }
inline std::string Projection::units() const { return projection_->units(); }
inline Projection::operator bool() const { return projection_->operator bool(); }

//---------------------------------------------------------------------------------------------------------------------
}
}
