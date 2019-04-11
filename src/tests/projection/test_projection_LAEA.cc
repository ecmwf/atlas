/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <vector>

#include "atlas/projection/Projection.h"
#include "atlas/util/Config.h"
#include "atlas/util/Point.h"
#include "atlas/grid.h"

#include "tests/AtlasTestEnvironment.h"

using atlas::util::Config;
using ListPointXY = std::vector<atlas::PointXY>;
namespace atlas {
namespace test {

static const double km = 1000.;

//-----------------------------------------------------------------------------

CASE( "test_projection_LAEA" ) {
    auto points_xy = ListPointXY { // in Meters!!!
      {0.,0.},
      {100.*km,0.},
      {0.,100.*km}
    };
    auto projection = Projection{ Config
            ("type","lambert_azimuthal_equal_area")
            ("central_longitude", -67.)
            ("standard_parallel", 50.) };
    EXPECT( projection.type() == "lambert_azimuthal_equal_area" );
    EXPECT( projection.units() == "meters" );
    EXPECT( projection.strictlyRegional() );
    
    for( auto pxy : points_xy ) {
        double tolerance_micrometers = 1.e-6;
        auto p_identity = projection.xy( projection.lonlat( pxy ) );
        EXPECT( is_approximately_equal( p_identity.x(), pxy.x(), tolerance_micrometers ) );
        EXPECT( is_approximately_equal( p_identity.y(), pxy.y(), tolerance_micrometers ) );
    }
}

//-----------------------------------------------------------------------------

CASE( "test_grid_creation_from_GRIB" ) {

    auto mir_parametrisation = []() -> Config {
        Config parametrisation;
        parametrisation.set("gridType",                             "lambert_azimuthal_equal_area");
        parametrisation.set("standardParallelInDegrees",             50. );
        parametrisation.set("centralLongitudeInDegrees",            -67. );
        parametrisation.set("latitudeOfFirstGridPointInDegrees",     75. );
        parametrisation.set("longitudeOfFirstGridPointInDegrees",   -90. );
        parametrisation.set("xDirectionGridLengthInMetres",          100.*km );
        parametrisation.set("yDirectionGridLengthInMetres",          100.*km );
        parametrisation.set("numberOfPointsAlongXAxis",              100 );
        parametrisation.set("numberOfPointsAlongYAxis",              50  );
        parametrisation.set("radius", util::Earth::radius() );
        return parametrisation;
    };

    auto from_GRIB = []( const eckit::Parametrisation& parametrisation ) -> Grid::Spec {
        // Convert parametrisation to Grid::Spec
        std::vector<double> firstPointLonLat(2);
        double standard_parallel, central_longitude, dx, dy;
        long nx, ny;
        double radius = util::Earth::radius();
        
        parametrisation.get("standardParallelInDegrees", standard_parallel);
        parametrisation.get("centralLongitudeInDegrees", central_longitude);
        parametrisation.get("radius", radius);
        parametrisation.get("longitudeOfFirstGridPointInDegrees", firstPointLonLat[0]);
        parametrisation.get("latitudeOfFirstGridPointInDegrees", firstPointLonLat[1]);
        parametrisation.get("numberOfPointsAlongXAxis", nx);
        parametrisation.get("numberOfPointsAlongYAxis", ny);
        parametrisation.get("xDirectionGridLengthInMetres", dx);
        parametrisation.get("yDirectionGridLengthInMetres", dy);
        
        Grid::Spec gridspec; // could be JSON or YAML
        gridspec.set("type","regional"); // --> indicates following values will be parsed
        gridspec.set("nx",nx);
        gridspec.set("ny",ny);
        gridspec.set("dx",dx);
        gridspec.set("dy",dy);
        gridspec.set("lonlat(xmin,ymax)",firstPointLonLat);
        gridspec.set("projection", [&]() {
            Grid::Spec projection;
            projection.set("type","lambert_azimuthal_equal_area");
            projection.set("standard_parallel",standard_parallel);
            projection.set("central_longitude",central_longitude);
            projection.set("radius",radius);
            return projection;
        }());
        return gridspec;
    };
    
    auto g = Grid( from_GRIB( mir_parametrisation() ) );
    auto domain = RectangularDomain( g.domain() );
    auto firstPointXY = PointXY{ domain.xmin(), domain.ymax() };
    auto firstPointLonLat = g.projection().lonlat( firstPointXY );
    EXPECT( is_approximately_equal( firstPointLonLat.lon(), mir_parametrisation().getDouble("longitudeOfFirstGridPointInDegrees") ) );
    EXPECT( is_approximately_equal( firstPointLonLat.lat(), mir_parametrisation().getDouble("latitudeOfFirstGridPointInDegrees") ) );

    
    // Get circumscribed LonLat bounding box ( Could be turned into memberfunction of StructuredGrid )
    {
        double n, s, w, e;

        double maxlat = -90.;
        double minlat =  90.;
        double minlon =  360.;
        double maxlon = -360.;
        auto structured = StructuredGrid(g);
        for( idx_t j : {0,structured.ny()-1} ) {
            for( idx_t i=0; i<structured.nx(j); ++i ) {
                double lat = structured.lonlat(i,j).lat();
                minlat = std::min( minlat, lat );
                maxlat = std::max( maxlat, lat );
            }
        }
        for( idx_t j=0; j<structured.ny(); ++j ) {
            for( idx_t i : {0,structured.nx(j)-1} ) {
                double lon = structured.lonlat(i,j).lon();
                minlon = std::min( minlon, lon );
                maxlon = std::max( maxlon, lon );
            }
        }
        n = maxlat;
        s = minlat;
        w = minlon;
        e = maxlon;

        Log::info() << "n = " << n << std::endl;
        Log::info() << "s = " << s << std::endl;
        Log::info() << "w = " << w << std::endl;
        Log::info() << "e = " << e << std::endl;
    }
}

//-----------------------------------------------------------------------------

}  // namespace test
}  // namespace atlas

int main( int argc, char** argv ) {
    return atlas::test::run( argc, argv );
}
