/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor
 * does it submit to any jurisdiction.
 */

#include "atlas/mesh/Mesh.h"

namespace atlas {

//----------------------------------------------------------------------------------------------------------------------

Mesh::Mesh() : impl_( new Implementation() ) {}

Mesh::Mesh( const Mesh& mesh ) : impl_( mesh.impl_ ) {}

Mesh::Mesh( const Implementation* impl ) : impl_( const_cast<Implementation*>( impl ) ) {}

Mesh::Mesh( eckit::Stream& stream ) : impl_( new Implementation( stream ) ) {}

//----------------------------------------------------------------------------------------------------------------------

}  // namespace atlas
