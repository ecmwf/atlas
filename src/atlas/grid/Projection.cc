#include "atlas/grid/Projection.h"

#include "eckit/config/Parametrisation.h"
#include "atlas/grid/detail/projection/LonLatProjection.h"

namespace atlas {
namespace grid {

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

void Projection::hash( eckit::MD5& md5 ) const {
    return projection_->hash(md5);
}

} // namespace Grid
} // namespace atlas
