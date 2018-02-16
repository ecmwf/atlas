/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "atlas/projection/detail/LonLatProjection.h"

namespace atlas {
namespace projection {
namespace detail {

template <typename Rotation>
LonLatProjectionT<Rotation>::LonLatProjectionT( const eckit::Parametrisation& config ) :
    ProjectionImpl(),
    rotation_( config ) {}

template <typename Rotation>
typename LonLatProjectionT<Rotation>::Spec LonLatProjectionT<Rotation>::spec() const {
    Spec proj_spec;
    proj_spec.set( "type", static_type() );
    rotation_.spec( proj_spec );
    return proj_spec;
}

template <typename Rotation>
void LonLatProjectionT<Rotation>::hash( eckit::Hash& hsh ) const {
    hsh.add( static_type() );
    rotation_.hash( hsh );
}

register_BuilderT1( ProjectionImpl, LonLatProjection, LonLatProjection::static_type() );
register_BuilderT1( ProjectionImpl, RotatedLonLatProjection, RotatedLonLatProjection::static_type() );

}  // namespace detail
}  // namespace projection
}  // namespace atlas
