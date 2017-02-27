#include "atlas/grid/ptr/Spacing.h"

namespace atlas {
namespace grid {
namespace ptr {

Spacing::Spacing():
    spacing_( nullptr ){
}

Spacing::Spacing( atlas::grid::spacing::Spacing *spacing):
    spacing_( spacing ) {
}

Spacing::Spacing( const eckit::Parametrisation& p ):
    spacing_( atlas::grid::spacing::Spacing::create(p) ) {
}

} // namespace ptr
} // namespace Grid
} // namespace atlas
