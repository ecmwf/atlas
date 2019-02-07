
#include <ostream>

#include "atlas/domain/detail/DomainFactory.h"
#include "atlas/domain/detail/ZonalBandDomain.h"
#include "atlas/runtime/Exception.h"

namespace atlas {
namespace domain {

namespace {

static bool _is_global( double ymin, double ymax ) {
    const double eps = 1.e-12;
    return std::abs( ( ymax - ymin ) - 180. ) < eps;
}

static std::array<double, 2> get_interval_y( const eckit::Parametrisation& params ) {
    double ymin, ymax;

    if ( !params.get( "ymin", ymin ) ) throw_Exception( "ymin missing in Params", Here() );

    if ( !params.get( "ymax", ymax ) ) throw_Exception( "ymax missing in Params", Here() );

    return {ymin, ymax};
}

/*
constexpr std::array<double, 2> interval_x() {
    return {0., 360.};
}
*/
}  // namespace

constexpr char ZonalBandDomain::units_[];

bool ZonalBandDomain::is_global( const Interval& y ) {
    const double eps = 1.e-12;
    return std::abs( std::abs( y[1] - y[0] ) - 180. ) < eps;
}

ZonalBandDomain::ZonalBandDomain( const eckit::Parametrisation& params ) :
    ZonalBandDomain( get_interval_y( params ) ) {}

ZonalBandDomain::ZonalBandDomain( const Interval& interval_y ) : ZonalBandDomain( interval_y, /*west*/ 0. ) {}

ZonalBandDomain::ZonalBandDomain( const Interval& interval_y, const double west ) :
    RectangularDomain( {west, west + 360.}, interval_y, units_ ) {
    global_   = _is_global( ymin(), ymax() );
    ymin_tol_ = ymin() - 1.e-6;
    ymax_tol_ = ymax() + 1.e-6;
}

bool ZonalBandDomain::contains( double x, double y ) const {
    return contains_y( y );
}

ZonalBandDomain::Spec ZonalBandDomain::spec() const {
    Spec domain_spec;
    domain_spec.set( "type", type() );
    domain_spec.set( "ymin", ymin() );
    domain_spec.set( "ymax", ymax() );
    return domain_spec;
}

void ZonalBandDomain::hash( eckit::Hash& h ) const {
    spec().hash( h );
}

void ZonalBandDomain::print( std::ostream& os ) const {
    os << "ZonalBandDomain["
       << "ymin=" << ymin() << ",ymax=" << ymax() << "]";
}

bool ZonalBandDomain::containsNorthPole() const {
    return ymax_tol_ >= 90.;
}

bool ZonalBandDomain::containsSouthPole() const {
    return ymin_tol_ <= -90.;
}

namespace {
static DomainBuilder<ZonalBandDomain> register_builder( ZonalBandDomain::static_type() );
}

}  // namespace domain
}  // namespace atlas
