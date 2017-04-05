/*
 * (C) Copyright 1996-2017 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "atlas/grid/Domain.h"

#include "atlas/grid/detail/domain/Domain.h"
#include "atlas/grid/detail/domain/EmptyDomain.h"
#include "atlas/grid/detail/domain/RectangularDomain.h"
#include "atlas/grid/detail/domain/ZonalBandDomain.h"
#include "atlas/grid/detail/domain/GlobalDomain.h"

using RD = atlas::grid::domain::RectangularDomain;

namespace atlas {
namespace grid {

Domain::Domain():
    domain_( new domain::EmptyDomain() ) {
}

Domain::Domain( const Domain& domain):
    domain_( domain.domain_ ) {
}

Domain::Domain( const domain_t* domain):
    domain_( domain ) {
}

Domain::Domain( const eckit::Parametrisation& p ):
    domain_( atlas::grid::domain::Domain::create(p) ) {
}

RectangularDomain::RectangularDomain( const Interval& x, const Interval& y, const std::string& units ) :
  Domain( 
    ( RD::is_global(x,y,units) ) ? new atlas::grid::domain::GlobalDomain() :
    ( RD::is_zonal_band(x,units) ? new atlas::grid::domain::ZonalBandDomain(y) :
    new atlas::grid::domain::RectangularDomain(x,y,units) ) ) {
}

} // namespace Grid
} // namespace atlas
