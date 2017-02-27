#include "atlas/grid/ptr/Projection.h"

#include "eckit/config/Parametrisation.h"

namespace atlas {
namespace grid {
namespace ptr {

Projection::Projection():
    projection_( atlas::grid::projection::Projection::create() ) {
}

Projection::Projection(const Projection& projection):
    projection_( projection.projection_ ) {
}

Projection::Projection( const projection::Projection *projection ):
    projection_( const_cast<atlas::grid::projection::Projection*>(projection) ) {
}

Projection::Projection( const eckit::Parametrisation& p ):
    projection_( atlas::grid::projection::Projection::create(p) ) {
}

} // namespace ptr
} // namespace Grid
} // namespace atlas
