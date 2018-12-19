#include "eckit/exception/Exceptions.h"
#include "eckit/utils/Hash.h"

#include "atlas/domain/detail/DomainFactory.h"
#include "atlas/domain/detail/EmptyDomain.h"

namespace atlas {
namespace domain {

EmptyDomain::EmptyDomain() {}

EmptyDomain::EmptyDomain( const eckit::Parametrisation& p ) {}

EmptyDomain::Spec EmptyDomain::spec() const {
    Spec domain_spec;
    domain_spec.set( "type", type() );
    return domain_spec;
}

void EmptyDomain::print( std::ostream& os ) const {
    os << "EmptyDomain";
}

void EmptyDomain::hash( eckit::Hash& h ) const {
    h.add( type() );
}

std::string EmptyDomain::units() const {
    NOTIMP;
}

namespace {
static DomainBuilder<EmptyDomain> register_builder( EmptyDomain::static_type() );
}

}  // namespace domain
}  // namespace atlas
