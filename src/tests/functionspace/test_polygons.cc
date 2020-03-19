/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "atlas/array/ArrayView.h"
#include "atlas/array/MakeView.h"
#include "atlas/field/Field.h"
#include "atlas/functionspace/StructuredColumns.h"
#include "atlas/grid/Partitioner.h"
#include "atlas/grid/StructuredGrid.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/util/LonLatPolygon.h"

#include "tests/AtlasTestEnvironment.h"

using namespace eckit;
using namespace atlas::functionspace;
using namespace atlas::util;

namespace atlas {
namespace test {

//-----------------------------------------------------------------------------

CASE( "test_polygons" ) {
    std::string gridname = eckit::Resource<std::string>( "--grid", "O32" );
    Grid grid( gridname );
    functionspace::StructuredColumns fs( grid );

    auto polygons = util::LonLatPolygons( fs.polygons() );

    std::vector<int> sizes( mpi::size() );
    std::vector<int> simplified_sizes( mpi::size() );
    for ( idx_t i = 0; i < mpi::size(); ++i ) {
        sizes[i]            = fs.polygons()[i].size();
        simplified_sizes[i] = polygons[i].size();
    }

    // Test iterator:
    for ( auto& polygon : fs.polygons() ) {
        Log::info() << "size of polygon = " << polygon.size() << std::endl;
    }

    for ( auto& polygon : polygons ) {
        Log::info() << "size of lonlatpolygon = " << polygon.size() << std::endl;
    }

    auto points = std::vector<PointLonLat>{
        {45., 75.},  {90., 75.},  {135., 75.},  {180., 75.},  {225., 75.},  {270., 75.},  {315., 75.},
        {45., 15.},  {90., 15.},  {135., 15.},  {180., 15.},  {225., 15.},  {270., 15.},  {315., 15.},
        {45., -75.}, {90., -75.}, {135., -75.}, {180., -75.}, {225., -75.}, {270., -75.}, {315., -75.},
    };

    std::vector<int> part( points.size() );
    for ( size_t n = 0; n < points.size(); ++n ) {
        Log::debug() << n << "  " << points[n];

        // A brute force approach.
        // This could be enhanced by a kd-tree search to nearest polygon centroid
        for ( idx_t p = 0; p < polygons.size(); ++p ) {
            if ( polygons[p].contains( points[n] ) ) {
                Log::debug() << " : " << p;
                part[n] = p;
            }
        }
        Log::debug() << std::endl;
    }

    if ( mpi::size() == 1 ) {
        auto expected_part  = std::vector<int>{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
        auto expected_sizes = std::vector<int>{132};
        auto expected_simplified_sizes = std::vector<int>{5};

        EXPECT( part == expected_part );
        EXPECT( sizes == expected_sizes );
        EXPECT( simplified_sizes == expected_simplified_sizes );
    }
    if ( mpi::size() == 4 ) {
        auto expected_part  = std::vector<int>{0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3};
        auto expected_sizes = std::vector<int>{50, 84, 84, 50};
        auto expected_simplified_sizes = std::vector<int>{7, 43, 43, 7};

        EXPECT( part == expected_part );
        EXPECT( sizes == expected_sizes );
        EXPECT( simplified_sizes == expected_simplified_sizes );
    }
}

//-----------------------------------------------------------------------------

}  // namespace test
}  // namespace atlas

int main( int argc, char** argv ) {
    return atlas::test::run( argc, argv );
}
