#include "atlas/grid/Domain.h"

#include "atlas/grid/detail/domain/Domain.h"

namespace atlas {
namespace grid {

Domain::Domain():
    domain_( nullptr ){
}

Domain::Domain( const Domain& domain):
    domain_( domain.domain_ ) {
}

Domain::Domain( const atlas::grid::domain::Domain *domain):
    domain_( const_cast<atlas::grid::domain::Domain*>(domain) ) {
}

Domain::Domain( const eckit::Parametrisation& p ):
    domain_( atlas::grid::domain::Domain::create(p) ) {
}

} // namespace Grid
} // namespace atlas
