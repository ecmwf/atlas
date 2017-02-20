
#include "atlas/grid/spacing/Spacing.h"
#include "eckit/config/Parametrisation.h"
#include "eckit/exception/Exceptions.h"

namespace atlas {
namespace grid {
namespace spacing {

Spacing* Spacing::create(const eckit::Parametrisation& params) {
  std::string spacingType;
  if (not params.get("type",spacingType) ) {
    throw eckit::BadParameter("type missing in configuration",Here());
  }
  return eckit::Factory<Spacing>::instance().get(spacingType).create(params);
}

}  // namespace spacing
}  // namespace grid
}  // namespace atlas
