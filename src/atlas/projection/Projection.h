#pragma once

#include <string>
#include "eckit/config/Parametrisation.h"
#include "eckit/value/Properties.h"
#include "eckit/memory/SharedPtr.h"

#include "atlas/util/Point.h"
#include "atlas/projection/detail/Projection.h"

//---------------------------------------------------------------------------------------------------------------------

// Forward declarations
namespace eckit {
class Parametrisation;
class MD5;
}

//---------------------------------------------------------------------------------------------------------------------

namespace atlas {

//---------------------------------------------------------------------------------------------------------------------

class Projection {
  
public:
  
  using Implementation = projection::detail::Projection;

public:

  Projection();
  Projection( const Projection& );
  Projection( const Implementation* );
  Projection( const eckit::Parametrisation& );

  // operator Implementation*() { return projection_.get(); }
  // operator const Implementation&() const { return *projection_.get(); }

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

  eckit::SharedPtr<Implementation> projection_;
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
