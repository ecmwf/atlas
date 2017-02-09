
#include "atlas/grid/spacing/Spacing.h"
#include "eckit/config/Parametrisation.h"
#include "eckit/exception/Exceptions.h"

namespace atlas {
namespace grid {
namespace spacing {

Spacing* Spacing::create(const eckit::Parametrisation& params) {
  std::string spacingType;
  if (params.get("spacingType",spacingType)) {
    return eckit::Factory<Spacing>::instance().get(spacingType).create(params);
  }

  // should return error here
  throw eckit::BadParameter("spacingType missing in params",Here());
  return NULL;
}

}  // namespace spacing
}  // namespace grid
}  // namespace atlas
