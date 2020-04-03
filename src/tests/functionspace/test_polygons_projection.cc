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
#include "atlas/output/Gmsh.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/util/LonLatPolygon.h"
#include "atlas/util/PolygonLocator.h"

#include "tests/AtlasTestEnvironment.h"

using namespace atlas::functionspace;
using namespace atlas::util;

namespace atlas {
namespace test {

namespace {
std::string grid_name() {
    static std::string _gridname = eckit::Resource<std::string>( "--grid", "L90" );
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

Grid grid() {
    static Grid _grid{
        option::name( grid_name() ) | util::Config( "projection", util::Config( "type", "rotated_lonlat" )(
                                                                      "south_pole", std::vector<double>{45., -60.} ) ),
        RectangularDomain{{0., 90}, {-45, 45.}}};
    return _grid;
}

Mesh mesh() {
    static Mesh _mesh{StructuredMeshGenerator{util::Config( "partitioner", "checkerboard" )}.generate( grid() )};
    return _mesh;
}


FunctionSpace structured_columns() {
    static functionspace::StructuredColumns _fs{grid(), grid::Partitioner( "checkerboard" )};
    return _fs;
}

FunctionSpace node_columns() {
    static functionspace::NodeColumns _fs{mesh()};
    return _fs;
}

FunctionSpace functionspace() {
    static FunctionSpace _fs = functionspace_name() == "NodeColumns" ? node_columns() : structured_columns();
    return _fs;
}


std::vector<PointLonLat>& points() {
    static std::vector<PointLonLat> _points{{90., 45.}, {135., 45.}, {60., 0.}, {90., 0.}};
    return _points;
}

void check_part( const std::vector<int>& vec ) {
    Log::debug() << "part = " << vec << std::endl;
    std::vector<int> expected;
    if ( mpi::size() == 1 && grid_name() == "L90" ) {
        expected = std::vector<int>{0, 0, 0, 0};
    }
    if ( mpi::size() == 4 && grid_name() == "L90" && functionspace_name() == "StructuredColumns" ) {
        expected = std::vector<int>{0, 1, 2, 3};
    }
    if ( mpi::size() == 4 && grid_name() == "L90" && functionspace_name() == "NodeColumns" ) {
        expected = std::vector<int>{0, 1, 2, 3};
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
    if ( mpi::size() == 1 && grid_name() == "L90" && functionspace_name() == "StructuredColumns" ) {
        expected = std::vector<int>{5};
    }
    if ( mpi::size() == 4 && grid_name() == "L90" && functionspace_name() == "StructuredColumns" ) {
        expected = std::vector<int>{7, 9, 5, 5};
    }
    if ( mpi::size() == 1 && grid_name() == "L90" && functionspace_name() == "NodeColumns" ) {
        expected = std::vector<int>{361};
    }
    if ( mpi::size() == 4 && grid_name() == "L90" && functionspace_name() == "NodeColumns" ) {
        expected = std::vector<int>{185, 183, 181, 179};
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
    if ( mpi::size() == 1 && grid_name() == "L90" && functionspace_name() == "StructuredColumns" ) {
        expected = std::vector<int>{5};
    }
    if ( mpi::size() == 4 && grid_name() == "L90" && functionspace_name() == "StructuredColumns" ) {
        expected = std::vector<int>{7, 9, 5, 5};
    }
    if ( mpi::size() == 1 && grid_name() == "L90" && functionspace_name() == "NodeColumns" ) {
        expected = std::vector<int>{5};
    }
    if ( mpi::size() == 4 && grid_name() == "L90" && functionspace_name() == "NodeColumns" ) {
        expected = std::vector<int>{7, 9, 5, 5};
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

CASE( "info" ) {
    Log::info() << grid().spec() << std::endl;
    functionspace().polygon().outputPythonScript( "rotated_polygon.py" );
    output::Gmsh{"mesh.msh", util::Config( "coordinates", "lonlat" )}.write( mesh() );
}

//-----------------------------------------------------------------------------

CASE( "test_polygons" ) {
    auto fs = functionspace();

    ATLAS_TRACE( "computations after setup" );
    auto polygons = util::LonLatPolygons( fs.polygons() );

    auto projection = fs.projection();
    EXPECT( projection );

    std::vector<int> sizes( mpi::size() );
    std::vector<int> simplified_sizes( mpi::size() );
    for ( idx_t i = 0; i < mpi::size(); ++i ) {
        sizes[i]            = fs.polygons()[i].size();
        simplified_sizes[i] = polygons[i].size();
    }

    // Test iterator:
    Log::info() << "Iterating PartitionPolygons" << std::endl;
    for ( auto& polygon : fs.polygons() ) {
        Log::info() << "  size of polygon = " << polygon.size() << std::endl;
    }

    Log::info() << "Iterating LonLatPolygons" << std::endl;
    for ( auto& polygon : polygons ) {
        Log::info() << "  size of lonlatpolygon = " << polygon.size() << std::endl;
    }

    std::vector<int> part( points().size() );
    for ( size_t n = 0; n < points().size(); ++n ) {
        Log::debug() << n << "  " << points()[n];

        // A brute force approach.
        // Use PolygonLocator class for optimized result
        for ( idx_t p = 0; p < polygons.size(); ++p ) {
            if ( polygons[p].contains( projection.xy( points()[n] ) ) ) {
                Log::debug() << " : " << p;
                part[n] = p;
            }
        }
        Log::debug() << std::endl;
    }

    check_part( part );
    check_sizes( sizes );
    check_simplified_sizes( simplified_sizes );

    PolygonLocator find_partition( polygons, projection );
    for ( size_t n = 0; n < points().size(); ++n ) {
        int found_partition;
        EXPECT_NO_THROW( found_partition = find_partition( points()[n] ) );
        EXPECT_EQ( found_partition, part[n] );
    }
}

//-----------------------------------------------------------------------------

}  // namespace test
}  // namespace atlas

int main( int argc, char** argv ) {
    return atlas::test::run( argc, argv );
}
