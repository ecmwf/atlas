/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "atlas/domain/Domain.h"
#include "atlas/domain/detail/Domain.h"
#include "atlas/domain/detail/EmptyDomain.h"
#include "atlas/domain/detail/GlobalDomain.h"
#include "atlas/domain/detail/RectangularDomain.h"
#include "atlas/domain/detail/ZonalBandDomain.h"

using RD = atlas::domain::RectangularDomain;
using ZD = atlas::domain::ZonalBandDomain;

namespace atlas {

Domain::Domain( const eckit::Parametrisation& p ) : Handle( atlas::domain::Domain::create( p ) ) {}

RectangularDomain::RectangularDomain( const Interval& x, const Interval& y, const std::string& units ) :
    Domain( ( RD::is_global( x, y, units ) )
                ? new atlas::domain::GlobalDomain( x[0] )
                : ( RD::is_zonal_band( x, units ) ? new atlas::domain::ZonalBandDomain( y, x[0] )
                                                  : new atlas::domain::RectangularDomain( x, y, units ) ) ),
    domain_( dynamic_cast<const atlas::domain::RectangularDomain*>( get() ) ) {}

RectangularDomain::RectangularDomain( const Domain& domain ) :
    Domain( domain ), domain_( dynamic_cast<const atlas::domain::RectangularDomain*>( get() ) ) {}

bool RectangularDomain::contains_x( double x ) const {
    return domain_->contains_x( x );
}

bool RectangularDomain::contains_y( double y ) const {
    return domain_->contains_y( y );
}

bool RectangularDomain::zonal_band() const {
    return domain_->zonal_band();
}

double RectangularDomain::xmin() const {
    return domain_->xmin();
}

double RectangularDomain::xmax() const {
    return domain_->xmax();
}

double RectangularDomain::ymin() const {
    return domain_->ymin();
}

double RectangularDomain::ymax() const {
    return domain_->ymax();
}

ZonalBandDomain::ZonalBandDomain( const Interval& y ) :
    RectangularLonLatDomain( ( ZD::is_global( y ) ) ? new atlas::domain::GlobalDomain()
                                                    : new atlas::domain::ZonalBandDomain( y ) ),
    domain_( dynamic_cast<const atlas::domain::ZonalBandDomain*>( get() ) ) {}

ZonalBandDomain::ZonalBandDomain( const Interval& y, const double& west ) :
    RectangularLonLatDomain( ( ZD::is_global( y ) ) ? new atlas::domain::GlobalDomain( west )
                                                    : new atlas::domain::ZonalBandDomain( y, west ) ),
    domain_( dynamic_cast<const atlas::domain::ZonalBandDomain*>( get() ) ) {}

ZonalBandDomain::ZonalBandDomain( const Domain& domain ) :
    RectangularLonLatDomain( domain ), domain_( dynamic_cast<const atlas::domain::ZonalBandDomain*>( get() ) ) {}

GlobalDomain::GlobalDomain( const double& west ) :
    ZonalBandDomain( new atlas::domain::GlobalDomain( west ) ),
    domain_( dynamic_cast<const atlas::domain::GlobalDomain*>( get() ) ) {}

GlobalDomain::GlobalDomain() :
    ZonalBandDomain( new atlas::domain::GlobalDomain() ),
    domain_( dynamic_cast<const atlas::domain::GlobalDomain*>( get() ) ) {}

GlobalDomain::GlobalDomain( const Domain& domain ) :
    ZonalBandDomain( domain ), domain_( dynamic_cast<const atlas::domain::GlobalDomain*>( get() ) ) {}


std::string atlas::Domain::type() const {
    return get()->type();
}

bool Domain::contains( double x, double y ) const {
    return get()->contains( x, y );
}

bool Domain::contains( const PointXY& p ) const {
    return get()->contains( p );
}

std::string Domain::units() const {
    return get()->units();
}

Domain::Spec Domain::spec() const {
    return get()->spec();
}

bool Domain::global() const {
    return get()->global();
}

bool Domain::empty() const {
    return get()->empty();
}

void Domain::hash( eckit::Hash& h ) const {
    get()->hash( h );
}

bool Domain::containsNorthPole() const {
    return get()->containsNorthPole();
}

bool Domain::containsSouthPole() const {
    return get()->containsSouthPole();
}

void Domain::print( std::ostream& os ) const {
    return get()->print( os );
}

std::ostream& operator<<( std::ostream& os, const Domain& d ) {
    d.print( os );
    return os;
}

}  // namespace atlas
