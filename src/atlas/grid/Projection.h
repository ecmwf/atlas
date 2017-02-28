#pragma once

#include <string>
#include <array>
#include "eckit/config/Parametrisation.h"
#include "eckit/value/Properties.h"
#include "eckit/memory/Builder.h"
#include "eckit/memory/Owned.h"

#include "atlas/util/Point.h"
#include "atlas/grid/detail/projection/Projection.h"

//---------------------------------------------------------------------------------------------------------------------

// Forward declarations
namespace eckit {
class Parametrisation;
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

  bool isStrictlyRegional() const;

  eckit::Properties spec() const;
  
  std::string units() const;
  
  operator bool() const;
  
  std::string type() const { return projection_->type(); }
  
private:

  eckit::SharedPtr<atlas::grid::projection::Projection> projection_;
};

//---------------------------------------------------------------------------------------------------------------------

inline void Projection::xy2lonlat(double crd[]) const { return projection_->xy2lonlat(crd); }
inline void Projection::lonlat2xy(double crd[]) const { return projection_->lonlat2xy(crd); }
inline PointLonLat Projection::lonlat( const PointXY& xy) const { return projection_->lonlat(xy); }
inline PointXY Projection::xy( const PointLonLat& lonlat) const { return projection_->xy(lonlat); }
inline bool Projection::isStrictlyRegional() const { return projection_->isStrictlyRegional(); }
inline eckit::Properties Projection::spec() const { return projection_->spec(); }
inline std::string Projection::units() const { return projection_->units(); }
inline Projection::operator bool() const { return projection_->operator bool(); }

//---------------------------------------------------------------------------------------------------------------------
}
}
