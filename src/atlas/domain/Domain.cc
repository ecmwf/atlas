/*
 * (C) Copyright 1996-2017 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "atlas/domain/Domain.h"

#include "atlas/domain/detail/Domain.h"
#include "atlas/domain/detail/EmptyDomain.h"
#include "atlas/domain/detail/RectangularDomain.h"
#include "atlas/domain/detail/ZonalBandDomain.h"
#include "atlas/domain/detail/GlobalDomain.h"

using RD = atlas::domain::RectangularDomain;

namespace atlas {

Domain::Domain():
    domain_( new domain::EmptyDomain() ) {
}

Domain::Domain( const Domain& domain):
    domain_( domain.domain_ ) {
}

Domain::Domain( const Implementation* domain):
    domain_( domain ) {
}

Domain::Domain( const eckit::Parametrisation& p ):
    domain_( atlas::domain::Domain::create(p) ) {
}

RectangularDomain::RectangularDomain( const Interval& x, const Interval& y, const std::string& units ) :
  Domain( 
    ( RD::is_global(x,y,units) ) ? new atlas::domain::GlobalDomain() :
    ( RD::is_zonal_band(x,units) ? new atlas::domain::ZonalBandDomain(y) :
    new atlas::domain::RectangularDomain(x,y,units) ) ) {
}

} // namespace atlas
