/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <algorithm>
#include "atlas/functionspace/NodeColumns.h"
#include "atlas/functionspace/EdgeColumns.h"
#include "atlas/grid/Grid.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/mesh/Nodes.h"
#include "atlas/mesh/HybridElements.h"
#include "atlas/meshgenerator/StructuredMeshGenerator.h"
#include "atlas/parallel/mpi/mpi.h"
#include "tests/AtlasTestEnvironment.h"

using namespace atlas::functionspace;
using namespace atlas::grid;
using namespace atlas::meshgenerator;

namespace atlas {
namespace test {

template< typename Container >
Container reversed( const Container& a ) {
    Container a_reversed = a;
    std::reverse( a_reversed.begin(), a_reversed.end() );
    return a_reversed;
}

//-----------------------------------------------------------------------------

CASE( "halo nodes" ) {

    Grid grid( "O8" );
    std::vector<int> halos {0, 1, 2, 3, 4};
    std::vector<int> nodes {560, 592, 624, 656, 688};

    for( bool reduce : {false, true} ) {

        SECTION( std::string( reduce ? "reduced" : "increased") ) {

            Mesh mesh = StructuredMeshGenerator().generate( grid );
            EXPECT( mesh.nodes().size() == nodes[0] );

            for( int h : ( reduce ? reversed(halos) : halos ) ) {
                NodeColumns fs( mesh, option::halo(h) );
                if( mpi::comm().size() == 1 ) {
                    EXPECT( fs.nb_nodes() == nodes[h] );
                }
            }
        }
    }
}

//-----------------------------------------------------------------------------

}  // namespace test
}  // namespace atlas

int main( int argc, char** argv ) {
    return atlas::test::run( argc, argv );
}
