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
    size_t root          = 0;
    std::string gridname = eckit::Resource<std::string>( "--grid", "O32" );
    Grid grid( gridname );
    functionspace::StructuredColumns fs( grid );

    auto lonlat_polygons = []( const functionspace::StructuredColumns& fs ) {
        const auto& polygons = fs.polygons();
        std::vector<std::unique_ptr<util::LonLatPolygon>> ll_poly( polygons.size() );
        for ( idx_t p = 0; p < polygons.size(); ++p ) {
            ll_poly[p].reset( new util::LonLatPolygon( *polygons[p] ) );
        }
        return ll_poly;
    };

    auto polygons = lonlat_polygons( fs );

    auto points = std::vector<PointLonLat>{
        {45., 75.},  {90., 75.},  {135., 75.},  {180., 75.},  {225., 75.},  {270., 75.},  {315., 75.},
        {45., 15.},  {90., 15.},  {135., 15.},  {180., 15.},  {225., 15.},  {270., 15.},  {315., 15.},
        {45., -75.}, {90., -75.}, {135., -75.}, {180., -75.}, {225., -75.}, {270., -75.}, {315., -75.},
    };

    std::vector<int> part( points.size() );
    for ( idx_t n = 0; n < points.size(); ++n ) {
        Log::debug() << n << "  " << points[n];
        for ( idx_t p = 0; p < polygons.size(); ++p ) {
            if ( polygons[p]->contains( points[n] ) ) {
                Log::debug() << " : " << p;
                part[n] = p;
            }
        }
        Log::debug() << std::endl;
    }

    if ( mpi::size() == 4 ) {
        auto expected_part = std::vector<int>{0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3};
        EXPECT( part == expected_part );
    }
}

//-----------------------------------------------------------------------------

}  // namespace test
}  // namespace atlas

int main( int argc, char** argv ) {
    return atlas::test::run( argc, argv );
}
