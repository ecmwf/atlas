
/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "eckit/config/Resource.h"

#include "atlas/array/ArrayView.h"
#include "atlas/array/MakeView.h"
#include "atlas/interpolation/method/PointIndex3.h"
#include "atlas/mesh/HybridElements.h"
#include "atlas/runtime/Trace.h"
#include "atlas/util/Topology.h"

namespace atlas {
namespace interpolation {
namespace method {

ElemIndex3* create_element_kdtree( const Mesh& mesh, const Field& field_centres ) {
    ATLAS_TRACE();
    const array::ArrayView<const double, 2> centres = array::make_view<double, 2>( field_centres );
    const array::ArrayView<const int, 1> flags      = array::make_view<int, 1>( mesh.cells().flags() );
    auto include_element                            = [&]( unsigned int e ) {
        using util::Topology;
        return not Topology::view( flags( e ) ).check( Topology::INVALID );
    };

    static bool fastBuildKDTrees = eckit::Resource<bool>( "$ATLAS_FAST_BUILD_KDTREES", true );

    ElemIndex3* tree      = new ElemIndex3();
    const size_t nb_cells = centres.shape( 0 );

    if ( fastBuildKDTrees ) {
        std::vector<ElemIndex3::Value> p;
        p.reserve( nb_cells );

        for ( unsigned int j = 0; j < nb_cells; ++j ) {
            if ( include_element( j ) ) {
                p.emplace_back( ElemIndex3::Point( centres( j, XX ), centres( j, YY ), centres( j, ZZ ) ),
                                ElemIndex3::Payload( j ) );
            }
        }

        tree->build( p.begin(), p.end() );
    }
    else {
        for ( unsigned int j = 0; j < nb_cells; ++j ) {
            if ( include_element( j ) ) {
                tree->insert(
                    ElemIndex3::Value( ElemIndex3::Point( centres( j, XX ), centres( j, YY ), centres( j, ZZ ) ),
                                       ElemIndex3::Payload( j ) ) );
            }
        }
    }
    return tree;
}

ElemIndex3* create_element_centre_index( const Mesh& mesh ) {
    return create_element_kdtree( mesh, mesh.cells().field( "centre" ) );
}

}  // namespace method
}  // namespace interpolation
}  // namespace atlas
