#include "eckit/exception/Exceptions.h"

#include "atlas/domain/detail/Domain.h"
#include "atlas/projection/Projection.h"

namespace atlas {
namespace domain {

Domain* Domain::create() {
    // default: global domain
    util::Config projParams;
    projParams.set( "type", "global" );
    return Domain::create( projParams );
}

Domain* Domain::create( const eckit::Parametrisation& p ) {
    std::string domain_type;
    if ( p.get( "type", domain_type ) ) { return eckit::Factory<Domain>::instance().get( domain_type ).create( p ); }

    // should return error here
    throw eckit::BadParameter( "type missing in Params", Here() );
}

}  // namespace domain
}  // namespace atlas
