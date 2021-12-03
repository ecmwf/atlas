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
#include <iomanip>
#include <sstream>

#include "atlas/grid/Iterator.h"
#include "atlas/grid/StructuredGrid.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/meshgenerator.h"
#include "atlas/output/Gmsh.h"
#include "atlas/runtime/Log.h"
#include "atlas/util/Config.h"

#include "tests/AtlasTestEnvironment.h"

using Config = atlas::util::Config;

namespace atlas {
namespace test {

//-----------------------------------------------------------------------------

namespace {
struct Increments {
    double dx;
    double dy;
    Increments(double _dx, double _dy): dx(_dx), dy(_dy) {}
};

static Grid RegularLL(const Increments& increments, const Domain& _domain) {
    // As MIR defines a "RegularLL" grid:
    //   - If the domain is global or zonal band, the grid is periodic, with the domain's West defining the West boundary
    //     and the East-most point is at (East-dx) where (East = West+360)
    //   - If the domain is not a zonal band or global, the East-most point is at (East)
    double dx = increments.dx;
    double dy = increments.dy;

    RectangularDomain domain{_domain};
    double n = domain.ymax();
    double s = domain.ymin();
    double w = domain.xmin();
    double e = domain.xmax();
    long nj  = long((n - s) / dy) + 1;
    long ni  = long((e - w) / dx) + long(domain.zonal_band() ? 0 : 1);

    using atlas::grid::LinearSpacing;
    StructuredGrid::XSpace xspace(LinearSpacing(w, e, ni, not domain.zonal_band()));
    StructuredGrid::YSpace yspace(LinearSpacing(n, s, nj));

    return StructuredGrid(xspace, yspace, StructuredGrid::Projection(), domain);
}

static Grid RegularLL(const Increments& increments) {
    return RegularLL(increments, RectangularDomain{{0., 360.}, {-90., 90.}});
}

static Grid RegularLL11(const Domain& domain) {
    return RegularLL(Increments{1., 1.}, domain);
}

static Grid RegularLL11() {
    return RegularLL(Increments{1., 1.}, RectangularDomain{{0., 360.}, {-90., 90.}});
}

}  // namespace

CASE("Test number of points for global grids") {
    // --- Global grids --- //
    EXPECT(Grid("O16").size() == 1600);

    EXPECT(Grid("O1280").size() == 6599680);

    EXPECT(RegularLL(Increments{1., 1.}).size() == 65160);
}

CASE("Test number of points for non-global grids") {
    Grid global = RegularLL11();

    auto test_regional_size = [&](const Domain& domain, idx_t size) -> bool {
        if (RegularLL11(domain).size() != size) {
            Log::error() << "RegularLL11( domain ).size() != size ---> [ " << RegularLL11(domain).size()
                         << " != " << size << " ]" << std::endl;
            Log::error() << "    with Domain = " << domain << std::endl;
            return false;
        }
        return true;
    };

    auto test_cropping_size = [&](const Domain& domain, idx_t size) -> bool {
        if (Grid(global, domain).size() != size) {
            Log::error() << "Grid( global, domain ).size() != size ---> [ " << Grid(global, domain).size()
                         << " != " << size << " ]" << std::endl;
            Log::error() << "    with Domain = " << domain << std::endl;
            return false;
        }
        return true;
    };

    // -- Parallels
    {
        // North Pole latitude
        EXPECT(test_regional_size(RectangularDomain{{0., 360.}, {90., 90.}}, 360));
        EXPECT(test_cropping_size(RectangularDomain{{0., 360.}, {90., 90.}}, 360));

        // North Pole latitude plus one extra
        EXPECT(test_regional_size(RectangularDomain{{0., 360.}, {89., 90.}}, 720));
        EXPECT(test_cropping_size(RectangularDomain{{0., 360.}, {89., 90.}}, 720));

        // Equator latitude
        EXPECT(test_regional_size(RectangularDomain{{0., 360.}, {0., 0.}}, 360));
        EXPECT(test_cropping_size(RectangularDomain{{0., 360.}, {0., 0.}}, 360));

        // South Pole latitude
        EXPECT(test_regional_size(RectangularDomain{{0., 360.}, {-90., -90.}}, 360));
        EXPECT(test_cropping_size(RectangularDomain{{0., 360.}, {-90., -90.}}, 360));

        // South Pole latitude plus one extra
        EXPECT(test_regional_size(RectangularDomain{{0., 360.}, {-90., -89.}}, 720));
        EXPECT(test_cropping_size(RectangularDomain{{0., 360.}, {-90., -89.}}, 720));
    }

    // -- Meridians
    {
        // Greenwich Meridian
        EXPECT(test_regional_size(RectangularDomain{{0., 0.}, {-90., 90.}}, 181));
        EXPECT(test_cropping_size(RectangularDomain{{0., 0.}, {-90., 90.}}, 181));

        // Greenwich Meridian + one to East
        EXPECT(test_regional_size(RectangularDomain{{0., 1.}, {-90., 90.}}, 2 * 181));
        EXPECT(test_cropping_size(RectangularDomain{{0., 1.}, {-90., 90.}}, 2 * 181));

        // Greenwich Meridian + one to West
        EXPECT(test_regional_size(RectangularDomain{{-1., 0.}, {-90., 90.}}, 2 * 181));
        EXPECT(test_cropping_size(RectangularDomain{{-1., 0.}, {-90., 90.}}, 2 * 181));

        // Meridian, one degree to the East of Greenwich
        EXPECT(test_regional_size(RectangularDomain{{-1., -1.}, {-90., 90.}}, 181));
        EXPECT(test_cropping_size(RectangularDomain{{-1., -1.}, {-90., 90.}}, 181));

        // Date line at 180 degrees East
        EXPECT(test_regional_size(RectangularDomain{{180., 180.}, {-90., 90.}}, 181));
        EXPECT(test_cropping_size(RectangularDomain{{180., 180.}, {-90., 90.}}, 181));

        // Date line at 180 degrees West
        EXPECT(test_regional_size(RectangularDomain{{-180., -180.}, {-90., 90.}}, 181));
        EXPECT(test_cropping_size(RectangularDomain{{-180., -180.}, {-90., 90.}}, 181));
    }

    // Crop to single point {0.,0.}
    EXPECT(test_cropping_size(RectangularDomain{{-0.5, 0.5}, {-0.5, 0.5}}, 1));
    EXPECT(test_cropping_size(RectangularDomain{{0., 0.}, {0., 0.}}, 1));

    // Regional grid with same domain leads to bounds given by domain --> 4 points
    EXPECT(test_regional_size(RectangularDomain{{-0.5, 0.5}, {-0.5, 0.5}}, 4));
    EXPECT(test_regional_size(RectangularDomain{{0., 0.}, {0., 0.}}, 1));
}

CASE("Test first point of cropped global grids (MIR-374)") {
    auto test_first_point = [](const std::string& name, const Domain& domain, const PointLonLat& p) -> bool {
        double eps = 1.e-6;
        Grid global_grid{name};
        Grid cropped_grid{global_grid, domain};

        PointLonLat first_point = *cropped_grid.lonlat().begin();
        if (not is_approximately_equal(first_point.lon(), p.lon(), eps) ||
            not is_approximately_equal(first_point.lat(), p.lat(), eps)) {
            auto old_precision = Log::error().precision(16);
            Log::error() << "First point doesn't match for Grid " << name << " with domain " << domain << std::endl;
            Log::error() << "    first_point = " << first_point << std::endl;
            Log::error() << "    expected    = " << p << std::endl;
            Log::error().precision(old_precision);
            return false;
        }
        return true;
    };

    // clang-format off
    EXPECT( test_first_point( "O16",  RectangularDomain{ {-180.,180.}, {-90.,90.} }, PointLonLat{-180.,85.76058712044382} ) );
    EXPECT( test_first_point( "O640", RectangularDomain{ {-180.,180.}, {-90.,90.} }, PointLonLat{-180.,89.89239644558999} ) );
    EXPECT( test_first_point( "O16",  RectangularDomain{ {-180.,180.}, {-60.,90.} }, PointLonLat{-180.,85.76058712044382} ) );
    EXPECT( test_first_point( "O640", RectangularDomain{ {-180.,180.}, {-60.,90.} }, PointLonLat{-180.,89.89239644558999} ) );
    EXPECT( test_first_point( "O16",  RectangularDomain{ {-180.,180.}, {-90.,60.} }, PointLonLat{-180.,58.14295404920328} ) );
    EXPECT( test_first_point( "O640", RectangularDomain{ {-180.,180.}, {-90.,60.} }, PointLonLat{-180.,59.953135752239} ) );
    EXPECT( test_first_point( "O16",  RectangularDomain{ {-10. ,10. }, {-10.,10.} }, PointLonLat{-9.473684210526358,8.306702856518804} ) );
    EXPECT( test_first_point( "O640", RectangularDomain{ {-10. ,10. }, {-10.,10.} }, PointLonLat{-9.878048780487802,9.910190568389} ) );
    // clang-format on
}

//-----------------------------------------------------------------------------

}  // namespace test
}  // namespace atlas

int main(int argc, char** argv) {
    return atlas::test::run(argc, argv);
}
