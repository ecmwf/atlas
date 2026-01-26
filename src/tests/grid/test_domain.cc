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

#include "atlas/domain.h"
#include "atlas/grid/StructuredGrid.h"
#include "atlas/runtime/Log.h"
#include "atlas/util/Config.h"

#include "atlas/mesh.h"
#include "atlas/functionspace/NodeColumns.h"

#include "tests/AtlasTestEnvironment.h"

namespace atlas {
namespace test {

//----------------------------------------------------------------------------------------------------------------------

CASE("test_domain_rectangular") {
    Domain domain = RectangularDomain({0, 180}, {-25, 25});
    EXPECT(not domain.global());
    EXPECT(domain.type() == std::string("rectangular"));

    util::Config domain_cfg = domain.spec();
    Domain from_cfg(domain_cfg);
    Log::info() << from_cfg.spec() << std::endl;
    EXPECT(from_cfg.type() == std::string("rectangular"));
}

CASE("test_domain_rectangular_tolerance") {
    using grid::LinearSpacing;

    Domain domain        = RectangularDomain({-27, 45}, {33, 73});
    RectangularDomain rd = domain;

    StructuredGrid::XSpace xspace(LinearSpacing(rd.xmin(), rd.xmax(), 721, true));
    StructuredGrid::YSpace yspace(LinearSpacing(rd.ymin(), rd.ymax(), 401));

    StructuredGrid grid(xspace, yspace, StructuredGrid::Projection(), rd);

    EXPECT(grid.nx(0) == 721);
    EXPECT(grid.ny() == 401);
}

CASE("test_domain_zonal_from_rectangular") {
    Domain domain = RectangularDomain({0, 360}, {-25, 25});
    EXPECT(not domain.global());
    EXPECT(domain.type() == std::string("zonal_band"));

    util::Config domain_cfg = domain.spec();
    Domain from_cfg(domain_cfg);
    Log::info() << from_cfg.spec() << std::endl;
    EXPECT(from_cfg.type() == std::string("zonal_band"));
}

CASE("test_domain_global_from_rectangular") {
    Domain domain = RectangularDomain({-180, 180}, {-90, 90});
    EXPECT(domain.global());
    EXPECT(domain.type() == std::string("global"));

    util::Config domain_cfg = domain.spec();
    Domain from_cfg(domain_cfg);
    Log::info() << from_cfg.spec() << std::endl;
    EXPECT(from_cfg.type() == std::string("global"));

    RectangularDomain rd = domain;
    EXPECT(rd == true);
    EXPECT(rd.xmin() == -180.);
    EXPECT(rd.xmax() == 180.);
    EXPECT(rd.ymin() == -90.);
    EXPECT(rd.ymax() == 90.);

    ZonalBandDomain zd = domain;
    EXPECT(zd == true);
    EXPECT(zd.xmin() == -180.);
    EXPECT(zd.xmax() == 180.);
    EXPECT(zd.ymin() == -90.);
    EXPECT(zd.ymax() == 90.);
}

CASE("test_domain_global_from_zonalband") {
    Domain domain = ZonalBandDomain({-45, 45});
    EXPECT(not domain.global());
    EXPECT(domain.type() == std::string("zonal_band"));

    util::Config domain_cfg = domain.spec();
    Domain from_cfg(domain_cfg);
    Log::info() << from_cfg.spec() << std::endl;
    EXPECT(from_cfg.type() == std::string("zonal_band"));

    RectangularDomain rd = domain;
    EXPECT(rd == true);
    EXPECT(rd.xmin() == 0);
    EXPECT(rd.xmax() == 360.);
    EXPECT(rd.ymin() == -45.);
    EXPECT(rd.ymax() == 45.);

    ZonalBandDomain zd = domain;
    EXPECT(zd == true);
    EXPECT(zd.ymin() == -45.);
    EXPECT(zd.ymax() == 45.);
}

CASE("test github #282") {
    size_t expected_points = 209019;
    auto computed_points = [](Domain domain) -> size_t {
        auto grid = Grid("L576x361", domain);
        auto fs = functionspace::NodeColumns(Mesh{grid}, option::halo(1));
        return fs.size();
    };
    EXPECT_EQ(expected_points, computed_points(Domain{}));
    EXPECT_EQ(expected_points, computed_points(GlobalDomain{}));
    EXPECT_EQ(expected_points, computed_points(GlobalDomain{-180}));
    EXPECT_EQ(expected_points, computed_points(ZonalBandDomain{{-90,90}}));
    EXPECT_EQ(expected_points, computed_points(ZonalBandDomain{{-90,90},-180}));
    EXPECT_EQ(expected_points, computed_points(Domain{util::Config
        ("type","zonal_band")
        ("ymin",-90)("ymax",90)}));
    EXPECT_EQ(expected_points, computed_points(Domain{util::Config
        ("type","zonal_band")
        ("west",-180)("ymin",-90)("ymax",90)}));
    EXPECT_EQ(expected_points, computed_points(RectangularDomain{{-180,180},{-90,90}}));
    EXPECT_EQ(expected_points, computed_points(Domain{util::Config
        ("type","rectangular")("units","degrees")
        ("xmin",0)   ("xmax",360)("ymin",-90)("ymax",90)}));
    EXPECT_EQ(expected_points, computed_points(Domain{util::Config
        ("type","rectangular")("units","degrees")
        ("xmin",-180)("xmax",180)("ymin",-90)("ymax",90)}));
}

//-----------------------------------------------------------------------------

}  // namespace test
}  // namespace atlas

int main(int argc, char** argv) {
    return atlas::test::run(argc, argv);
}
