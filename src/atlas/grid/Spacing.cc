#include "atlas/grid/Spacing.h"

namespace atlas {
namespace grid {

Spacing::Spacing():
  spacing_( nullptr ){
}

Spacing::Spacing( const Spacing& other ):
  spacing_( other.spacing_ ) {
}

Spacing::Spacing( const spacing::Spacing *spacing ):
    spacing_( spacing ) {
}

Spacing::Spacing( const eckit::Parametrisation& p ):
    spacing_( atlas::grid::spacing::Spacing::create(p) ) {
}

} // namespace Grid
} // namespace atlas
