/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "eckit/utils/MD5.h"

#include "atlas/grid.h"
#include "atlas/projection.h"
#include "atlas/util/Config.h"
#include "atlas/util/Point.h"

#include "tests/AtlasTestEnvironment.h"


namespace atlas {
namespace test {

//-----------------------------------------------------------------------------


Projection::Jacobian getJacobian(const Projection& proj, const PointLonLat& lonlat0) {
    double dlon = 0.0001, dlat = 0.0001;

    double xy0[2] = {lonlat0.lon() + 0.00, lonlat0.lat() + 0.00};
    double xy1[2] = {lonlat0.lon() + dlon, lonlat0.lat() + 0.00};
    double xy2[2] = {lonlat0.lon() + 0.00, lonlat0.lat() + dlat};

    proj.lonlat2xy(xy0);
    proj.lonlat2xy(xy1);
    proj.lonlat2xy(xy2);

    Projection::Jacobian jac;

    jac[0] = {(xy1[0] - xy0[0]) / dlon, (xy2[0] - xy0[0]) / dlat};
    jac[1] = {(xy1[1] - xy0[1]) / dlon, (xy2[1] - xy0[1]) / dlat};

    return jac;
}

void printJacobian(const Projection::Jacobian& jac, const std::string& name) {
    printf(" %s = \n", name.c_str());
    printf("  %12.4f, %12.4f | %12.4f\n", jac[0][0], jac[0][1], sqrt(jac[0][0] * jac[0][0] + jac[0][1] * jac[0][1]));
    printf("  %12.4f, %12.4f | %12.4f\n", jac[1][0], jac[1][1], sqrt(jac[1][0] * jac[1][0] + jac[1][1] * jac[1][1]));
    printf("  %12.4f, %12.4f\n", sqrt(jac[0][0] * jac[0][0] + jac[1][0] * jac[1][0]),
           sqrt(jac[0][1] * jac[0][1] + jac[1][1] * jac[1][1]));
}

template <typename SKIP>
void doTaylorTest(const StructuredGrid& grid, double dmax, SKIP skip) {
    const auto& proj = grid.projection();

    eckit::MD5 hash;

    for (int j = 0, jglo = 0; j < grid.ny(); j++) {
        for (int i = 0; i < grid.nx(j); i++, jglo++) {
            if (!skip(grid, i, j)) {
                double lonlat[2];
                grid.lonlat(i, j, lonlat);

                auto jacA = proj.jacobian(PointLonLat(lonlat));
                auto jacB = getJacobian(proj, PointLonLat(lonlat));
                auto jacC = jacB * jacA.inverse() - Projection::Jacobian::identity();

                hash.add(jacA[0][0]);
                hash.add(jacA[0][1]);
                hash.add(jacA[1][0]);
                hash.add(jacA[1][1]);

                double diff = jacC.norm();

                if (diff > dmax) {
                    printf("------------------\n");
                    printf(" i, j = %8d, %8d\n", i, j);
                    printJacobian(jacA, "jacA");
                    printf("\n");
                    printJacobian(jacB, "jacB");
                    printf("\n");
                    printf("%12.4e\n", diff);
                    printf("\n");
                }

                EXPECT(diff < dmax);
            }
        }
    }
    Log::info() << "Jacobian checksum (not crossplatform) : " << hash.digest() << std::endl;
}


void doTaylorTest(const StructuredGrid& grid, double dmax) {
    auto noskip = [](const StructuredGrid& grid, int i, int j) { return false; };
    doTaylorTest(grid, dmax, noskip);
}

CASE("test_rotated_schmidt") {
    auto projection = [] {
        util::Config config;
        config.set("type", "rotated_schmidt");
        config.set("stretching_factor", 2.4);
        config.set("rotation_angle", 0.0);
        config.set("north_pole", {2.0, 46.7});
        return Projection{config};
    }();

    auto jacobian = projection.jacobian(PointLonLat{0., 0.});
    double tol    = 1.e-10;
    EXPECT_APPROX_EQ(jacobian[0][0], 1.372570713035, tol);
    EXPECT_APPROX_EQ(jacobian[0][1], -0.045140586761, tol);
    EXPECT_APPROX_EQ(jacobian[1][0], 0.045110961287, tol);
    EXPECT_APPROX_EQ(jacobian[1][1], 1.371669903783, tol);

    doTaylorTest(Grid{"N16", projection}, 1e-3,
                 // Avoid large jumps of pseudo-longitude (ie 359.5 -> 0)
                 [](const StructuredGrid& grid, int i, int j) {
                     return (i == 0) || (i == grid.nx(j) - 1) || (i == grid.nx(j) / 2);
                 });
}


CASE("test_lambert_conformal_conic") {
    const int Nx = 64, Ny = 64, Nux = 53, Nuy = 53;
    const double DxInMetres = 50000., DyInMetres = 50000.;
    const double LaDInDegrees = 46.2, Latin1InDegrees = 46.2, Latin2InDegrees = 46.2, LoVInDegrees = 2.0;
    const double XMinInMetres = -Nux / 2 * DxInMetres, YMinInMetres = -Nuy / 2 * DyInMetres;

    auto projection = [&] {
        util::Config config;
        config.set("type", "lambert_conformal_conic");
        config.set("longitude0", LoVInDegrees);
        config.set("latitude0", LaDInDegrees);
        config.set("latitude1", Latin1InDegrees);
        config.set("latitude2", Latin2InDegrees);
        return Projection{config};
    }();

    auto jacobian = projection.jacobian(PointLonLat{0., 0.});
    double tol    = 1.e-5;
    EXPECT_APPROX_EQ(jacobian[0][0], 148530.135993682081, tol);
    EXPECT_APPROX_EQ(jacobian[0][1], 3742.887654051571, tol);
    EXPECT_APPROX_EQ(jacobian[1][0], -3742.887654051572, tol);
    EXPECT_APPROX_EQ(jacobian[1][1], 148530.135993682081, tol);

    util::Config grid;
    grid.set("type", "regional");
    grid.set("dx", DxInMetres);
    grid.set("dy", DyInMetres);
    grid.set("xmin", XMinInMetres);
    grid.set("ymin", YMinInMetres);
    grid.set("nx", Nx);
    grid.set("ny", Ny);
    grid.set("projection", projection.spec());

    doTaylorTest(Grid(grid), 1e-6);
}

CASE("test_lonlat") {
    doTaylorTest(StructuredGrid("N16"), 1e-9);
}

CASE("test_constructors") {
    Projection::Jacobian a = {};
    Projection::Jacobian b = {1., 2., 3., 4.};
    Projection::Jacobian c = {{5., 6.}, {7., 8.}};

    Log::info() << "a =" << std::endl;
    Log::info() << a << std::endl;
    Log::info() << "b =" << std::endl;
    Log::info() << b << std::endl;
    Log::info() << "c =" << std::endl;
    Log::info() << c << std::endl;

    EXPECT_EQ(a[0][0], 0.);
    EXPECT_EQ(a[0][1], 0.);
    EXPECT_EQ(a[1][0], 0.);
    EXPECT_EQ(a[1][1], 0.);

    EXPECT_EQ(b[0][0], 1.);
    EXPECT_EQ(b[0][1], 2.);
    EXPECT_EQ(b[1][0], 3.);
    EXPECT_EQ(b[1][1], 4.);

    EXPECT_EQ(c[0][0], 5.);
    EXPECT_EQ(c[0][1], 6.);
    EXPECT_EQ(c[1][0], 7.);
    EXPECT_EQ(c[1][1], 8.);
}

CASE("test_multiplication") {
    const Projection::Jacobian A = {{0., -1.}, {1., 0.}};
    const Point2 x               = {2., 1.};

    const auto Ax      = A * x;
    const auto AinvAx  = A.inverse() * Ax;
    const auto Atimes2 = A * 2.;

    Log::info() << "A =" << std::endl;
    Log::info() << A << std::endl;
    Log::info() << "2A =" << std::endl;
    Log::info() << Atimes2 << std::endl;
    Log::info() << "x =" << std::endl;
    Log::info() << x << std::endl;
    Log::info() << "Ax =" << std::endl;
    Log::info() << Ax << std::endl;
    Log::info() << "A^-1 Ax =" << std::endl;
    Log::info() << AinvAx << std::endl;

    EXPECT_EQ(Atimes2[0][0], 0.);
    EXPECT_EQ(Atimes2[0][1], -2.);
    EXPECT_EQ(Atimes2[1][0], 2.);
    EXPECT_EQ(Atimes2[1][1], 0.);
    EXPECT_EQ(Ax, Point2(-1., 2.));
    EXPECT_EQ(AinvAx, x);
}


//-----------------------------------------------------------------------------

}  // namespace test
}  // namespace atlas

int main(int argc, char** argv) {
    return atlas::test::run(argc, argv);
}
