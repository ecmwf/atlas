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
#include "atlas/functionspace/NodeColumns.h"
#include "atlas/functionspace/StructuredColumns.h"
#include "atlas/grid/Partitioner.h"
#include "atlas/grid/StructuredGrid.h"
#include "atlas/meshgenerator.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/util/PolygonLocator.h"
#include "atlas/util/PolygonXY.h"

#include "tests/AtlasTestEnvironment.h"

using namespace atlas::functionspace;
using namespace atlas::util;

namespace atlas {
namespace test {

namespace {
std::string grid_name() {
    static std::string _gridname = eckit::Resource<std::string>( "--grid", "O32" );
    return _gridname;
}

std::string functionspace_name() {
    static std::string _name = eckit::Resource<std::string>( "--functionspace", "StructuredColumns" );
    return _name;
}

std::string configuration() {
    std::stringstream out;
    out << "mpi::size()=" << mpi::size() << ",grid_name()=" << grid_name()
        << ",functionspace_name()=" << functionspace_name();
    return out.str();
};


FunctionSpace structured_columns() {
    static functionspace::StructuredColumns _fs = []() {
        Grid grid( grid_name() );
        return functionspace::StructuredColumns{grid};
    }();
    return _fs;
}

FunctionSpace node_columns() {
    static functionspace::NodeColumns _fs = []() {
        Grid grid( grid_name() );
        Mesh mesh = StructuredMeshGenerator().generate( grid );
        return functionspace::NodeColumns{mesh};
    }();
    return _fs;
}


FunctionSpace functionspace() {
    static FunctionSpace _fs = functionspace_name() == "NodeColumns" ? node_columns() : structured_columns();
    return _fs;
}


std::vector<PointLonLat>& points() {
    static std::vector<PointLonLat> _points{
        {45., 75.},  {90., 75.},  {135., 75.},  {180., 75.},  {225., 75.},  {270., 75.},  {315., 75.},
        {45., 15.},  {90., 15.},  {135., 15.},  {180., 15.},  {225., 15.},  {270., 15.},  {315., 15.},
        {45., -75.}, {90., -75.}, {135., -75.}, {180., -75.}, {225., -75.}, {270., -75.}, {315., -75.},
    };
    return _points;
}

void check_part( const std::vector<int>& vec ) {
    Log::debug() << "part = " << vec << std::endl;
    std::vector<int> expected;
    if ( mpi::size() == 1 && grid_name() == "O32" ) {
        expected = std::vector<int>{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    }
    if ( mpi::size() == 4 && grid_name() == "O32" && functionspace_name() == "StructuredColumns" ) {
        expected = std::vector<int>{0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3};
    }
    if ( mpi::size() == 4 && grid_name() == "O32" && functionspace_name() == "NodeColumns" ) {
        expected = std::vector<int>{0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3};
    }
    if ( expected.size() ) {
        EXPECT( vec == expected );
    }
    else {
        Log::warning() << "Check for part not implemented for configuration " << configuration() << std::endl;
    }
}

void check_sizes( const std::vector<int>& vec ) {
    Log::debug() << "sizes = " << vec << std::endl;
    std::vector<int> expected;
    if ( mpi::size() == 1 && grid_name() == "O32" && functionspace_name() == "StructuredColumns" ) {
        expected = std::vector<int>{5};
    }
    if ( mpi::size() == 4 && grid_name() == "O32" && functionspace_name() == "StructuredColumns" ) {
        expected = std::vector<int>{7, 43, 43, 7};
    }
    if ( mpi::size() == 1 && grid_name() == "O32" && functionspace_name() == "NodeColumns" ) {
        expected = std::vector<int>{167};
    }
    if ( mpi::size() == 4 && grid_name() == "O32" && functionspace_name() == "NodeColumns" ) {
        expected = std::vector<int>{169, 147, 149, 165};
    }
    if ( expected.size() ) {
        EXPECT( vec == expected );
    }
    else {
        Log::warning() << "Check for sizes not implemented for configuration " << configuration() << std::endl;
    }
}

void check_simplified_sizes( const std::vector<int>& vec ) {
    Log::debug() << "simplified_sizes = " << vec << std::endl;
    std::vector<int> expected;
    if ( mpi::size() == 1 && grid_name() == "O32" && functionspace_name() == "StructuredColumns" ) {
        expected = std::vector<int>{5};
    }
    if ( mpi::size() == 4 && grid_name() == "O32" && functionspace_name() == "StructuredColumns" ) {
        expected = std::vector<int>{7, 43, 43, 7};
    }
    if ( mpi::size() == 1 && grid_name() == "O32" && functionspace_name() == "NodeColumns" ) {
        expected = std::vector<int>{5};
    }
    if ( mpi::size() == 4 && grid_name() == "O32" && functionspace_name() == "NodeColumns" ) {
        expected = std::vector<int>{8, 5, 8, 7};
    }
    if ( expected.size() ) {
        EXPECT( vec == expected );
    }
    else {
        Log::warning() << "Check for simplified_sizes not implemented for configuration " << configuration()
                       << std::endl;
    }
}


}  // namespace

//-----------------------------------------------------------------------------

CASE( "info" ) {
    functionspace().polygon().outputPythonScript( "polygon.py" );
}

CASE( "test_polygons" ) {
    auto fs = functionspace();

    ATLAS_TRACE( "computations after setup" );
    auto polygons = ListPolygonXY( fs.polygons() );

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
        Log::info() << "size of PolygonXY = " << polygon.size() << std::endl;
    }

    std::vector<int> part( points().size() );
    for ( size_t n = 0; n < points().size(); ++n ) {
        Log::debug() << n << "  " << points()[n];

        // A brute force approach.
        // Use PolygonLocator class for optimized result
        for ( idx_t p = 0; p < polygons.size(); ++p ) {
            if ( polygons[p].contains( points()[n] ) ) {
                Log::debug() << " : " << p;
                part[n] = p;
            }
        }
        Log::debug() << std::endl;
    }

    check_part( part );
    check_sizes( sizes );
    check_simplified_sizes( simplified_sizes );

    PolygonLocator find_partition( polygons );
    for ( size_t n = 0; n < points().size(); ++n ) {
        EXPECT_EQ( find_partition( points()[n] ), part[n] );
    }
    Log::info() << std::endl;
}

CASE( "test_polygon_locator_from_const_ref_polygons" ) {
    auto polygons = ListPolygonXY{functionspace().polygons()};
    PolygonLocator find_partition( polygons );
    EXPECT_EQ( find_partition( PointLonLat{0., 90.} ), 0 );
    EXPECT_EQ( find_partition( PointLonLat{0., -90.} ), mpi::size() - 1 );
}

CASE( "test_polygon_locator_from_move" ) {
    PolygonLocator find_partition( ListPolygonXY{functionspace().polygons()} );
    EXPECT_EQ( find_partition( PointLonLat{0., 90.} ), 0 );
    EXPECT_EQ( find_partition( PointLonLat{0., -90.} ), mpi::size() - 1 );
}

CASE( "test_polygon_locator_from_shared" ) {
    auto polygons = std::make_shared<ListPolygonXY>( functionspace().polygons() );
    PolygonLocator find_partition( polygons );
    EXPECT_EQ( find_partition( PointLonLat{0., 90.} ), 0 );
    EXPECT_EQ( find_partition( PointLonLat{0., -90.} ), mpi::size() - 1 );
}

//-----------------------------------------------------------------------------

}  // namespace test
}  // namespace atlas

int main( int argc, char** argv ) {
    return atlas::test::run( argc, argv );
}
