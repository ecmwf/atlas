#pragma once

#include "eckit/memory/SharedPtr.h"

//---------------------------------------------------------------------------------------------------------------------

// Forward declarations
namespace eckit {
class Parametrisation;
}

namespace atlas {
namespace grid {
namespace domain {
class Domain;
}
}
}

//---------------------------------------------------------------------------------------------------------------------

namespace atlas {
namespace grid {
namespace ptr {

//---------------------------------------------------------------------------------------------------------------------

class Domain {

public:

  Domain();
  Domain( const Domain& );
  Domain( const atlas::grid::domain::Domain* );
  Domain( const eckit::Parametrisation& );

  operator atlas::grid::domain::Domain*() { return domain_.get(); }

private:

  eckit::SharedPtr<atlas::grid::domain::Domain> domain_;
};

//---------------------------------------------------------------------------------------------------------------------

}
}
}
