/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "atlas/projection/Projection.h"

namespace atlas {

Projection::Projection() : projection_( Implementation::create() ) {}

Projection::Projection( const Projection& projection ) : projection_( projection.projection_ ) {}

Projection::Projection( const Implementation* projection ) : projection_( const_cast<Implementation*>( projection ) ) {}

Projection::Projection( const eckit::Parametrisation& p ) : projection_( Implementation::create( p ) ) {}

void Projection::hash( eckit::Hash& h ) const {
    return projection_->hash( h );
}

}  // namespace atlas
