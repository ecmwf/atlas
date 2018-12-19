
#include <ostream>

#include "eckit/utils/Hash.h"

#include "atlas/domain/detail/DomainFactory.h"
#include "atlas/domain/detail/GlobalDomain.h"

namespace atlas {
namespace domain {

namespace {
constexpr std::array<double, 2> yrange() {
    return {-90., 90.};
}
}  // namespace

GlobalDomain::GlobalDomain( const double west ) : ZonalBandDomain( yrange(), west ) {}

GlobalDomain::GlobalDomain() : ZonalBandDomain( yrange() ) {}

GlobalDomain::GlobalDomain( const eckit::Parametrisation& ) : GlobalDomain() {}

GlobalDomain::Spec GlobalDomain::spec() const {
    Spec domain_spec;
    domain_spec.set( "type", type() );
    return domain_spec;
}

void GlobalDomain::hash( eckit::Hash& h ) const {
    h.add( type() );
}

void GlobalDomain::print( std::ostream& os ) const {
    os << "GlobalDomain";
}

namespace {
static DomainBuilder<GlobalDomain> register_builder( GlobalDomain::static_type() );
}
}  // namespace domain
}  // namespace atlas
