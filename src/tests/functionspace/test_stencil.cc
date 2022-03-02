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
#include "eckit/linalg/Vector.h"
#include "eckit/types/Types.h"

#include "atlas/array.h"
#include "atlas/field/Field.h"
#include "atlas/functionspace/NodeColumns.h"
#include "atlas/functionspace/PointCloud.h"
#include "atlas/functionspace/StructuredColumns.h"
#include "atlas/grid/Partitioner.h"
#include "atlas/grid/Stencil.h"
#include "atlas/grid/StencilComputer.h"
#include "atlas/grid/StructuredGrid.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/meshgenerator.h"
#include "atlas/output/Gmsh.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/util/CoordinateEnums.h"
#include "atlas/util/MicroDeg.h"

#include "tests/AtlasTestEnvironment.h"

using namespace eckit;
using namespace atlas::functionspace;
using namespace atlas::util;
using namespace atlas::grid;


namespace atlas {
namespace test {

//-----------------------------------------------------------------------------

std::vector<double> IFS_vertical_coordinates(idx_t nlev) {
    std::vector<double> zcoord(nlev);
    double dzcoord = 1. / double(nlev);
    for (idx_t jlev = 0; jlev < nlev; ++jlev) {
        zcoord[jlev] = 0.5 * dzcoord + jlev * dzcoord;
    }
    return zcoord;
}

std::vector<double> zrange(idx_t nlev, double min, double max) {
    std::vector<double> zcoord(nlev);
    double dzcoord = (max - min) / double(nlev - 1);
    for (idx_t jlev = 0; jlev < nlev; ++jlev) {
        zcoord[jlev] = min + jlev * dzcoord;
    }
    return zcoord;
}

double cubic(double x, double min, double max) {
    double x0   = min;
    double x1   = 0.5 * (max + min);
    double x2   = max;
    double xmax = 0.5 * (x0 + x1);
    return (x - x0) * (x - x1) * (x - x2) / ((xmax - x0) * (xmax - x1) * (xmax - x2));
}


CASE("test finding of North-West grid point") {
    std::string gridname = eckit::Resource<std::string>("--grid", "O8");

    StructuredGrid grid(gridname);

    constexpr double tol = 0.5e-6;

    constexpr idx_t halo = 2;

    ComputeNorth compute_j_north(grid, halo);
    ComputeWest compute_i_west(grid, halo);

    struct IJ {
        idx_t i;
        idx_t j;
    };
    idx_t ny = grid.ny();

    if (mpi::comm().size() == 1) {
        auto entries = {
            std::make_tuple(PointXY{0. + 0.5 * tol, grid.y(0) + 0.5 * tol}, IJ{0, 0}),
            std::make_tuple(PointXY{0. - 0.5 * tol, grid.y(0) - 0.5 * tol}, IJ{0, 0}),
            std::make_tuple(PointXY{0. + 2.0 * tol, grid.y(0) + 2.0 * tol}, IJ{0, -1}),
            std::make_tuple(PointXY{0. - 2.0 * tol, grid.y(0) - 2.0 * tol}, IJ{-1, 0}),
            std::make_tuple(PointXY{360. + 0.5 * tol, grid.y(ny - 1) + 0.5 * tol}, IJ{grid.nx(ny - 1), ny - 1}),
            std::make_tuple(PointXY{360. - 0.5 * tol, grid.y(ny - 1) - 0.5 * tol}, IJ{grid.nx(ny - 1), ny - 1}),
            std::make_tuple(PointXY{360. + 2.0 * tol, grid.y(ny - 1) + 2.0 * tol}, IJ{grid.nx(ny - 2), ny - 2}),
            std::make_tuple(PointXY{360. - 2.0 * tol, grid.y(ny - 1) - 2.0 * tol}, IJ{grid.nx(ny - 1) - 1, ny - 1}),
        };
        for (auto entry : entries) {
            auto p     = std::get<0>(entry);
            auto index = std::get<1>(entry);
            EXPECT(compute_j_north(p.y()) == index.j);
            EXPECT(compute_i_west(p.x(), index.j) == index.i);
        }
    }
}

CASE("test horizontal stencil") {
    //if ( mpi::comm().size() == 1 ) {
    std::string gridname = eckit::Resource<std::string>("--grid", "O8");

    StructuredGrid grid(gridname);
    int halo = eckit::Resource<int>("--halo", 2);
    util::Config config;
    config.set("halo", halo);
    config.set("levels", 9);
    config.set("periodic_points", true);
    functionspace::StructuredColumns fs(grid, grid::Partitioner("equal_regions"), config);

    HorizontalStencil<4> stencil;

    ComputeHorizontalStencil compute_stencil(grid, stencil.width());

    auto departure_points = {
        PointXY(0., 90.),
        PointXY(0., -90.),
        PointXY(0., 0.),
        PointXY(360., 0.),
    };
    for (auto p : departure_points) {
        Log::info() << p << std::endl;
        compute_stencil(p.x(), p.y(), stencil);
        for (idx_t j = 0; j < stencil.width(); ++j) {
            Log::info() << stencil.i(0, j) << " " << stencil.j(j) << "   --   "
                        << "x,y = " << fs.compute_xy(stencil.i(0, j), stencil.j(j)) << std::endl;
        }
        Log::info() << std::endl;
    }
}

CASE("test horizontal stencil linear") {
    //if ( mpi::comm().size() == 1 ) {
    std::string gridname = eckit::Resource<std::string>("--grid", "O8");

    StructuredGrid grid(gridname);
    //int halo = eckit::Resource<int>( "--halo", 2 );
    util::Config config;
    config.set("levels", 9);
    config.set("periodic_points", true);
    functionspace::StructuredColumns fs(grid, grid::Partitioner("equal_regions"), config);

    HorizontalStencil<2> stencil;

    ComputeHorizontalStencil compute_stencil(grid, stencil.width());

    auto departure_points = {
        PointXY(0., 90.), PointXY(0., -90.), PointXY(0., 0.), PointXY(360., 0.), PointXY(359.9, 0.),
    };
    for (auto p : departure_points) {
        Log::info() << p << std::endl;
        compute_stencil(p.x(), p.y(), stencil);
        for (idx_t j = 0; j < stencil.width(); ++j) {
            for (idx_t i = 0; i < stencil.width(); ++i) {
                Log::info() << stencil.i(i, j) << " " << stencil.j(j) << "   --   "
                            << "x,y = " << fs.compute_xy(stencil.i(i, j), stencil.j(j)) << std::endl;
            }
            Log::info() << std::endl;
        }
    }
}


//-----------------------------------------------------------------------------

CASE("test vertical stencil") {
    SECTION("Initialize ComputeLower from raw data (as e.g. given from IFS)") {
        double z[] = {0.1, 0.3, 0.5, 0.7, 0.9};
        EXPECT_NO_THROW(ComputeLower{Vertical(5, array::make_view<double, 1>(z, 5), std::vector<double>{0., 1.})});
    }


    idx_t nlev     = 10;
    auto zcoord    = IFS_vertical_coordinates(nlev);
    double dzcoord = 1. / double(nlev);

    auto vertical = Vertical(nlev, zcoord, std::vector<double>{0., 1.});

    SECTION("Test compute_lower works as expected ") {
        ComputeLower compute_lower(vertical);

        const double eps = 1.e-14;

        for (idx_t k = 0; k < nlev; ++k) {
            idx_t k_expected = std::max<idx_t>(0, std::min(nlev - 1, k));
            EXPECT(compute_lower(zcoord[k]) == k_expected);
            EXPECT(compute_lower(zcoord[k] - eps) == k_expected);
            EXPECT(compute_lower(zcoord[k] + eps) == k_expected);
            EXPECT(compute_lower(zcoord[k] + 0.5 * dzcoord) == k_expected);
        }
    }

    SECTION("Compute vertical stencil") {
        ComputeVerticalStencil compute_vertical_stencil(vertical, 4);
        std::vector<double> departure_points{0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};
        for (auto p : departure_points) {
            VerticalStencil<4> stencil;
            compute_vertical_stencil(p, stencil);
            Log::info() << p << "   :    ";
            for (idx_t k = 0; k < stencil.width(); ++k) {
                Log::info() << stencil.k(k) << "[" << vertical[stencil.k(k)] << "] ";
            }
            Log::info() << "     interval = " << stencil.k_interval() << std::endl;
            if (p < vertical[0]) {
                EXPECT(stencil.k_interval() == -1);
            }
            else if (p < vertical[1]) {
                EXPECT(stencil.k_interval() == 0);
            }
            else if (p > vertical[nlev - 1]) {
                EXPECT(stencil.k_interval() == 3);
            }
            else if (p > vertical[nlev - 2]) {
                EXPECT(stencil.k_interval() == 2);
            }
            else {
                EXPECT(stencil.k_interval() == 1);
            }
        }
    }
}

//-----------------------------------------------------------------------------

#if 1
CASE("ifs method to find nearest grid point") {
    // see satrad/module/gaussgrid.F90
    std::string gridname = eckit::Resource<std::string>("--grid", "O8");
    StructuredGrid grid(gridname);

    auto p = PointXY{0., grid.y(0)};
    idx_t kgrib_lat, kgrib_lon;
    {
        double x = p.x();
        double y = p.y();
        std::vector<double> pdlat(grid.ny());
        for (idx_t j = 0; j < grid.ny(); ++j) {
            pdlat[j] = std::abs(y - grid.y(j));
        }

        auto iterator = std::min_element(pdlat.begin(), pdlat.end());
        kgrib_lat     = static_cast<idx_t>(iterator - pdlat.begin());

        double zfirstlon = grid.x(0, kgrib_lat);
        double zdlon     = grid.x(1, kgrib_lat) - zfirstlon;
        double zsafelon  = std::fmod(x - zfirstlon + 720., 360.);
        kgrib_lon        = static_cast<idx_t>(std::fmod(std::round(zsafelon / zdlon), grid.nx(kgrib_lat)));
    }
    EXPECT(kgrib_lon == 0);
    EXPECT(kgrib_lat == 0);
}
#endif
//-----------------------------------------------------------------------------

}  // namespace test
}  // namespace atlas

int main(int argc, char** argv) {
    return atlas::test::run(argc, argv);
}
