#include "atlas/projection/Projection.h"

namespace atlas {

Projection::Projection():
    projection_( Implementation::create() ) {
}

Projection::Projection(const Projection& projection):
    projection_( projection.projection_ ) {
}

Projection::Projection( const Implementation* projection ):
    projection_( const_cast<Implementation*>(projection) ) {
}

Projection::Projection( const eckit::Parametrisation& p ):
    projection_( Implementation::create(p) ) {
}

void Projection::hash( eckit::MD5& md5 ) const {
    return projection_->hash(md5);
}

} // namespace atlas
