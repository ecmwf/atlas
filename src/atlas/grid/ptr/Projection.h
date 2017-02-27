#pragma once

#include <string>
#include "atlas/grid/projection/Projection.h"

//---------------------------------------------------------------------------------------------------------------------

// Forward declarations
namespace eckit {
class Parametrisation;
}

//---------------------------------------------------------------------------------------------------------------------

namespace atlas {
namespace grid {
namespace ptr {

//---------------------------------------------------------------------------------------------------------------------

class Projection {

public:

  Projection();
  Projection( const Projection& );
  Projection( const atlas::grid::projection::Projection* );
  Projection( const eckit::Parametrisation& );

  operator atlas::grid::projection::Projection*() { return projection_.get(); }

  std::string units() const { return projection_.get()->units(); }

private:

  eckit::SharedPtr<atlas::grid::projection::Projection> projection_;
};

//---------------------------------------------------------------------------------------------------------------------

}
}
}
