/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "eckit/utils/Hash.h"

#include "atlas/option/Options.h"
#include "atlas/projection/Projection.h"
#include "atlas/projection/detail/LonLatProjection.h"
#include "atlas/projection/detail/ProjectionImpl.h"
#include "atlas/util/Config.h"

namespace atlas {

Projection::Projection() : Handle( new projection::detail::LonLatProjection() ) {}

Projection::Projection( const eckit::Parametrisation& p ) : Handle( Implementation::create( p ) ) {}

Projection::Projection( const std::string& type, const eckit::Parametrisation& p ) :
    Handle( Implementation::create( type, p ) ) {}

atlas::Projection::operator bool() const {
    return get()->operator bool();
}

void Projection::hash( eckit::Hash& h ) const {
    return get()->hash( h );
}

void atlas::Projection::xy2lonlat( double crd[] ) const {
    return get()->xy2lonlat( crd );
}

void atlas::Projection::xy2lonlat( Point2& point ) const {
    return get()->xy2lonlat( point );
}


void atlas::Projection::lonlat2xy( double crd[] ) const {
    return get()->lonlat2xy( crd );
}

void atlas::Projection::lonlat2xy( Point2& point ) const {
    return get()->lonlat2xy( point );
}

atlas::Projection::Jacobian atlas::Projection::jacobian( const PointLonLat& p ) const {
    return get()->jacobian( p );
}

PointLonLat atlas::Projection::lonlat( const PointXY& xy ) const {
    return get()->lonlat( xy );
}

PointXY atlas::Projection::xy( const PointLonLat& lonlat ) const {
    return get()->xy( lonlat );
}

bool atlas::Projection::strictlyRegional() const {
    return get()->strictlyRegional();
}

RectangularLonLatDomain Projection::lonlatBoundingBox( const Domain& domain ) const {
    return get()->lonlatBoundingBox( domain );
}

Projection::Spec atlas::Projection::spec() const {
    return get()->spec();
}

std::string atlas::Projection::units() const {
    return get()->units();
}

std::string Projection::type() const {
    return get()->type();
}

}  // namespace atlas
