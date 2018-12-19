/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "eckit/exception/Exceptions.h"
#include "eckit/utils/Hash.h"

#include "atlas/projection/detail/LonLatProjection.h"
#include "atlas/projection/detail/ProjectionFactory.h"
#include "atlas/projection/detail/ProjectionImpl.h"
#include "atlas/util/Config.h"


namespace atlas {
namespace projection {
namespace detail {

const ProjectionImpl* ProjectionImpl::create( const eckit::Parametrisation& p ) {
    std::string projectionType;
    if ( p.get( "type", projectionType ) ) { return ProjectionFactory::build( projectionType, p ); }

    // should return error here
    throw eckit::BadParameter( "type missing in Params", Here() );
}

Rotated::Rotated( const PointLonLat& south_pole, double rotation_angle ) :
    util::Rotation( south_pole, rotation_angle ) {}

Rotated::Rotated( const eckit::Parametrisation& p ) : util::Rotation( p ) {}

void Rotated::spec( Spec& s ) const {
    std::vector<double> npole{northPole().lon(), northPole().lat()};
    std::vector<double> spole{southPole().lon(), southPole().lat()};
    s.set( "north_pole", npole );
    s.set( "south_pole", spole );
    s.set( "rotation_angle", rotationAngle() );
}

void Rotated::hash( eckit::Hash& hsh ) const {
    hsh.add( "rotated" );
    hsh.add( southPole().lon() );
    hsh.add( southPole().lat() );
    hsh.add( rotationAngle() );
}

}  // namespace detail
}  // namespace projection
}  // namespace atlas
