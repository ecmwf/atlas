/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "atlas/mesh/actions/BuildXYZField.h"

#include "atlas/array/MakeView.h"
#include "atlas/field/Field.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/mesh/Nodes.h"
#include "atlas/runtime/Trace.h"
#include "atlas/util/Earth.h"
#include "atlas/util/Point.h"

namespace atlas {
namespace mesh {
namespace actions {

//----------------------------------------------------------------------------------------------------------------------

BuildXYZField::BuildXYZField( const std::string& name, bool force_recompute ) :
    name_( name ), force_recompute_( force_recompute ) {}

Field& BuildXYZField::operator()( Mesh& mesh ) const {
    return operator()( mesh.nodes() );
}

Field& BuildXYZField::operator()( mesh::Nodes& nodes ) const {
    bool recompute = force_recompute_;
    if ( !nodes.has_field( name_ ) ) {
        nodes.add( Field( name_, array::make_datatype<double>(), array::make_shape( nodes.size(), 3 ) ) );
        recompute = true;
    }
    if ( recompute ) {
        ATLAS_TRACE( "BuildXYZField" );
        array::ArrayView<double, 2> lonlat = array::make_view<double, 2>( nodes.lonlat() );
        array::ArrayView<double, 2> xyz    = array::make_view<double, 2>( nodes.field( name_ ) );

        PointXYZ p2;
        for ( idx_t n = 0; n < nodes.size(); ++n ) {
            const PointLonLat p1( lonlat( n, 0 ), lonlat( n, 1 ) );
            util::Earth::convertSphericalToCartesian( p1, p2 );
            xyz( n, 0 ) = p2.x();
            xyz( n, 1 ) = p2.y();
            xyz( n, 2 ) = p2.z();
        }
    }
    return nodes.field( name_ );
}

//----------------------------------------------------------------------------------------------------------------------

}  // namespace actions
}  // namespace mesh
}  // namespace atlas
