/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <cmath>

#include "atlas/grid.h"
#include "atlas/mesh.h"
#include "atlas/meshgenerator.h"

#include "atlas/runtime/Log.h"
#include "atlas/util/CoordinateEnums.h"
#include "atlas/mesh/actions/ReorderHilbert.h"

#include "tests/AtlasTestEnvironment.h"

namespace atlas {
namespace test {

//-----------------------------------------------------------------------------

CASE( "test_hilbert_reordering" ) {
    Mesh mesh = StructuredMeshGenerator().generate( Grid("O64" ) );
    mesh::actions::ReorderHilbert{mesh}();
    auto lonlat = array::make_view<double,2>( mesh.nodes().lonlat() );
    for( idx_t n=0; n< std::min(100,lonlat.shape(0)); ++n ) {
        Log::info() << n << '\t' << PointLonLat{ lonlat(n,0), lonlat(n,1 ) } << std::endl;
    }
}

//-----------------------------------------------------------------------------

}  // namespace test
}  // namespace atlas

int main( int argc, char** argv ) {
    return atlas::test::run( argc, argv );
}
