/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <ostream>
#include <utility>

#include "eckit/utils/Hash.h"

#include "atlas/domain/detail/DomainFactory.h"
#include "atlas/domain/detail/RectangularDomain.h"
#include "atlas/runtime/Exception.h"

namespace atlas {
namespace domain {

using Interval = RectangularDomain::Interval;

namespace {

static std::array<double, 2> get_interval_x( const eckit::Parametrisation& params ) {
    double xmin, xmax;

    if ( !params.get( "xmin", xmin ) ) {
        throw_Exception( "xmin missing in Params", Here() );
    }

    if ( !params.get( "xmax", xmax ) ) {
        throw_Exception( "xmax missing in Params", Here() );
    }

    return {xmin, xmax};
}

static std::array<double, 2> get_interval_y( const eckit::Parametrisation& params ) {
    double ymin, ymax;

    if ( !params.get( "ymin", ymin ) ) {
        throw_Exception( "ymin missing in Params", Here() );
    }

    if ( !params.get( "ymax", ymax ) ) {
        throw_Exception( "ymax missing in Params", Here() );
    }

    return {ymin, ymax};
}

static std::string get_units( const eckit::Parametrisation& params ) {
    std::string units;
    if ( !params.get( "units", units ) ) {
        throw_Exception( "units missing in Params", Here() );
    }
    return units;
}

}  // namespace

bool RectangularDomain::is_global( const Interval& x, const Interval& y, const std::string& units ) {
    if ( units != "degrees" ) {
        return false;
    }

    const double eps = 1.e-12;
    return std::abs( ( x[1] - x[0] ) - 360. ) < eps && std::abs( std::abs( y[1] - y[0] ) - 180. ) < eps;
}

bool RectangularDomain::is_zonal_band( const Interval& x, const std::string& units ) {
    if ( units != "degrees" ) {
        return false;
    }

    const double eps = 1.e-12;
    return std::abs( ( x[1] - x[0] ) - 360. ) < eps;
}

RectangularDomain::RectangularDomain( const eckit::Parametrisation& params ) :
    RectangularDomain( get_interval_x( params ), get_interval_y( params ), get_units( params ) ) {}

RectangularDomain::RectangularDomain( const Interval& x, const Interval& y, const std::string& units ) :
    xmin_( x[0] ), xmax_( x[1] ), ymin_( y[0] ), ymax_( y[1] ), units_( units ) {
    unit_degrees_ = ( units_ == "degrees" ) ? true : false;

    // Make sure xmax>=xmin and ymax>=ymin
    if ( xmin_ > xmax_ ) {
        std::swap( xmin_, xmax_ );
    }
    if ( ymin_ > ymax_ ) {
        std::swap( ymin_, ymax_ );
    }
    global_ = is_global( {xmin_, xmax_}, {ymin_, ymax_}, units_ );

    const double tol = 1.e-6;
    xmin_tol_        = xmin_ - tol;
    ymin_tol_        = ymin_ - tol;
    xmax_tol_        = xmax_ + tol;
    ymax_tol_        = ymax_ + tol;
}

bool RectangularDomain::contains( double x, double y ) const {
    return contains_x( x ) and contains_y( y );
}

RectangularDomain::Spec RectangularDomain::spec() const {
    Spec domain_spec;
    domain_spec.set( "type", type() );
    domain_spec.set( "xmin", xmin() );
    domain_spec.set( "xmax", xmax() );
    domain_spec.set( "ymin", ymin() );
    domain_spec.set( "ymax", ymax() );
    domain_spec.set( "units", units() );
    return domain_spec;
}

void RectangularDomain::print( std::ostream& os ) const {
    os << "RectangularDomain["
       << "xmin=" << xmin() << ",xmax=" << xmax() << ",ymin=" << ymin() << ",ymax=" << ymax() << ",units=" << units()
       << "]";
}

void RectangularDomain::hash( eckit::Hash& h ) const {
    double multiplier = units() == "meters" ? 1e2 : 1e8;
    auto add_double   = [&]( const double& x ) { h.add( std::round( x * multiplier ) ); };
    h.add( type() );
    h.add( units() );
    add_double( xmin() );
    add_double( xmax() );
    add_double( ymin() );
    add_double( ymax() );
}

bool RectangularDomain::containsNorthPole() const {
    return unit_degrees_ && ymax_tol_ >= 90.;
}

bool RectangularDomain::containsSouthPole() const {
    return unit_degrees_ && ymin_tol_ <= -90.;
}

namespace {
static DomainBuilder<RectangularDomain> register_builder( RectangularDomain::static_type() );
}

extern "C" {
double atlas__LonLatRectangularDomain__north( const RectangularDomain* This ) {
    return This->ymax();
}
double atlas__LonLatRectangularDomain__west( const RectangularDomain* This ) {
    return This->xmin();
}
double atlas__LonLatRectangularDomain__south( const RectangularDomain* This ) {
    return This->ymin();
}
double atlas__LonLatRectangularDomain__east( const RectangularDomain* This ) {
    return This->xmax();
}
}

}  // namespace domain
}  // namespace atlas
