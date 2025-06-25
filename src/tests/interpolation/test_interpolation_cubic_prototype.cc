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

#include "atlas/array.h"
#include "atlas/field/Field.h"
#include "atlas/functionspace/NodeColumns.h"
#include "atlas/functionspace/PointCloud.h"
#include "atlas/functionspace/StructuredColumns.h"
#include "atlas/grid/Partitioner.h"
#include "atlas/grid/StructuredGrid.h"
#include "atlas/interpolation.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/meshgenerator.h"
#include "atlas/output/Gmsh.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/util/CoordinateEnums.h"
#include "atlas/util/MicroDeg.h"

#if ATLAS_ECKIT_HAVE_ECKIT_585
#include "eckit/linalg/LinearAlgebraSparse.h"
#else
#include "eckit/linalg/LinearAlgebra.h"
#endif
#include "eckit/linalg/SparseMatrix.h"
#include "eckit/linalg/Vector.h"
#include "eckit/types/Types.h"


#include "CubicInterpolationPrototype.h"
#include "tests/AtlasTestEnvironment.h"

using namespace atlas::functionspace;
using namespace atlas::grid;
using namespace atlas::util;
using namespace atlas::array;

namespace atlas {
namespace test {

//-----------------------------------------------------------------------------

std::vector<double> IFS_full_levels_uniform(idx_t nlev) {
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

CASE("test vertical cubic interpolation") {
    idx_t nlev    = 10;
    auto vertical = Vertical{nlev, IFS_full_levels_uniform(nlev)};
    bool limiter  = true;
    CubicVerticalInterpolation interpolate(vertical, util::Config("limiter", limiter));
    std::vector<double> departure_points{0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 0.3246};
    array::ArrayT<double> array(nlev);
    auto view = array::make_view<double, 1>(array);
    Log::info() << "source:" << std::endl;
    for (idx_t k = 0; k < nlev; ++k) {
        view(k) = cubic(vertical(k), 0., 1.);
        Log::info() << "  " << vertical(k) << " : " << view(k) << std::endl;
    }

    Log::info() << "interpolation:" << std::endl;
    for (auto p : departure_points) {
        Log::info() << "  " << std::setw(6) << p << " :    " << std::setw(12) << interpolate(p, view)
                    << "     expected : " << cubic(p, 0., 1.) << std::endl;
        if (p >= vertical.front() && p <= vertical.back() && !limiter) {
            EXPECT(eckit::types::is_approximately_equal(interpolate(p, view), cubic(p, 0., 1.)));
        }
    }
}

//-----------------------------------------------------------------------------

CASE("test horizontal cubic interpolation") {
    //if ( mpi::comm().size() == 1 ) {
    std::string gridname = eckit::Resource<std::string>("--grid", "O8");

    StructuredGrid grid(gridname);
    int halo = eckit::Resource<int>("--halo", 2);
    Config config;
    config.set("halo", halo);
    config.set("levels", 9);
    config.set("periodic_points", true);
    StructuredColumns fs(grid, Partitioner("equal_regions"), config);

    Field field = fs.createField<double>(option::levels(0));
    auto f      = array::make_view<double, 1>(field);
    auto xy     = array::make_view<double, 2>(fs.xy());

    auto fx  = [](double x) { return cubic(x, 0., 360.); };
    auto fy  = [](double y) { return cubic(y, -90., 90.); };
    auto fxy = [fx, fy](double x, double y) { return fx(x) * fy(y); };

    for (idx_t j = 0; j < fs.size(); ++j) {
        f(j) = fxy(xy(j, XX), xy(j, YY));
    }

    CubicHorizontalInterpolation cubic_interpolation(fs);

    auto departure_points = {
        PointXY(0.13257, 45.6397),
    };
    for (auto p : departure_points) {
        Log::info() << p << "  -->  " << cubic_interpolation(p.x(), p.y(), f) << std::endl;
        EXPECT(is_approximately_equal(cubic_interpolation(p.x(), p.y(), f), fxy(p.x(), p.y())));
    }
}

//-----------------------------------------------------------------------------

bool operator==(const eckit::linalg::Triplet& t1, const eckit::linalg::Triplet& t2) {
    if (t1.row() != t2.row()) {
        return false;
    }
    if (t1.col() != t2.col()) {
        return false;
    }
    if (!is_approximately_equal(t1.value(), t2.value())) {
        return false;
    }

    return true;
}
bool operator!=(const eckit::linalg::Triplet& t1, const eckit::linalg::Triplet& t2) {
    return !(t1 == t2);
}
bool operator==(const std::vector<eckit::linalg::Triplet>& t1, const std::vector<eckit::linalg::Triplet>& t2) {
    if (t1.size() != t2.size()) {
        return false;
    }
    for (size_t i = 0; i < t1.size(); ++i) {
        if (t1[i] != t2[i]) {
            return false;
        }
    }
    return true;
}
//std::ostream& operator<<( std::ostream& out, const LocalView<double, 1>& array ) {
//    out << "{ ";
//    for ( idx_t i = 0; i < array.size(); ++i ) {
//        out << array( i );
//        if ( i < array.size() - 1 ) {
//            out << ",";
//        }
//        out << " ";
//    }
//    out << "}";
//    return out;
//}

//-----------------------------------------------------------------------------

CASE("test horizontal cubic interpolation triplets") {
    //if ( mpi::comm().size() == 1 ) {
    std::string gridname = eckit::Resource<std::string>("--grid", "O8");

    StructuredGrid grid(gridname);
    int halo = eckit::Resource<int>("--halo", 2);
    Config config;
    config.set("halo", halo);
    config.set("levels", 9);
    config.set("periodic_points", true);
    StructuredColumns fs(grid, Partitioner("equal_regions"), config);

    Field field = fs.createField<double>(option::levels(0));
    auto f      = array::make_view<double, 1>(field);
    auto xy     = array::make_view<double, 2>(fs.xy());

    for (idx_t j = 0; j < fs.size(); ++j) {
        f(j) = xy(j, XX);
    }

    CubicHorizontalInterpolation cubic_interpolation(fs);

    auto departure_points = PointCloud{{
        {0.13257, 45.6397},
        {360., -90.},
    }};
    auto departure_lonlat = make_view<double, 2>(departure_points.lonlat());

    CubicHorizontalInterpolation::WorkSpace ws;
    CubicHorizontalInterpolation::Triplets triplets;
    for (idx_t row = 0; row < departure_points.size(); ++row) {
        auto triplets_row =
            cubic_interpolation.compute_triplets(row, departure_lonlat(row, XX), departure_lonlat(row, YY), ws);
        Log::info() << departure_lonlat.slice(row, array::Range::all()) << "  -->  {\n";
        for (auto triplet : triplets_row) {
            Log::info() << "    " << triplet << " ,\n";
        }
        Log::info() << " } " << std::endl;
        std::copy(triplets_row.begin(), triplets_row.end(), std::back_inserter(triplets));
    }

    auto triplets2 = cubic_interpolation.reserve_triplets(departure_points.size());
    {
        idx_t row{0};
        for (auto p : departure_points.iterate().xy()) {
            cubic_interpolation.insert_triplets(row++, p.x(), p.y(), triplets2, ws);
        }
    }


    EXPECT(triplets2 == triplets);

    eckit::linalg::SparseMatrix matrix(departure_points.size(), fs.size(), triplets2);
    Log::info() << matrix << std::endl;

    std::vector<double> tgt(departure_points.size());
    eckit::linalg::Vector v_src(const_cast<double*>(f.data()), f.size());
    eckit::linalg::Vector v_tgt(tgt.data(), tgt.size());
#if ATLAS_ECKIT_HAVE_ECKIT_585
    eckit::linalg::LinearAlgebraSparse::backend().spmv(matrix, v_src, v_tgt);
#else
    eckit::linalg::LinearAlgebra::backend().spmv(matrix, v_src, v_tgt);
#endif
    Log::info() << "output = " << tgt << std::endl;
}

//-----------------------------------------------------------------------------

CASE("test 3d cubic interpolation") {
    const double tolerance = 1.e-15;

    //if ( mpi::comm().size() == 1 ) {
    std::string gridname = eckit::Resource<std::string>("--grid", "O8");
    idx_t nlev           = 11;

    Vertical vertical(nlev, zrange(nlev, 0., 1.));
    Log::info() << zrange(nlev, 0., 1.) << std::endl;
    ;

    StructuredGrid grid(gridname);
    int halo = eckit::Resource<int>("--halo", 2);
    Config config;
    config.set("halo", halo);
    config.set("periodic_points", true);
    StructuredColumns fs(grid, vertical, Partitioner("equal_regions"), config);

    Field input   = fs.createField<double>();
    auto f        = array::make_view<double, 2>(input);
    auto xy       = array::make_view<double, 2>(fs.xy());
    const auto& z = fs.vertical();

    auto fx = [](double x) { return cubic(x, 0., 360.); };
    auto fy = [](double y) { return cubic(y, -90., 90.); };
    auto fz = [](double z) { return cubic(z, 0., 1.); };
    auto fp = [fx, fy, fz](const PointXYZ& p) { return fx(p.x()) * fy(p.y()) * fz(p.z()); };

    for (idx_t n = 0; n < fs.size(); ++n) {
        for (idx_t k = fs.k_begin(); k < fs.k_end(); ++k) {
            PointXYZ p{
                xy(n, XX),
                xy(n, YY),
                z(k),
            };
            f(n, k) = fp(p);
        }
    }
    input.set_dirty(false);  // to avoid halo-exchange


    auto departure_points = PointCloud({
        {90., -45., 0.25}, {0., -45., 0.25},  {180., -45., 0.25}, {360., -45., 0.25}, {90., -90., 0.25},
        {90., 0., 0.25},   {90., 90., 0.25},  {10., -45., 0.25},  {20., -45., 0.25},  {30., -45., 0.25},
        {40., -45., 0.25}, {50., -45., 0.25}, {60., -45., 0.25},  {70., -45., 0.25},  {80., -45., 0.25},
        {90., -45., 0.25}, {10., -10., 0.25}, {20., -20., 0.25},  {30., -30., 0.25},  {40., -40., 0.25},
        {50., -50., 0.25}, {60., -60., 0.25}, {70., -70., 0.25},  {80., -80., 0.25},  {90., -90., 0.25},
        {60., -60., 0.6},  {90., -45., 0.16}, {90., -45., 0.6},   {90., -45., 1.},    {90., -45., 0.},
        {90., -45., 0.1},  {45., 33., 0.7},   {359., 20., 0.1},   {360., 88., 0.9},   {0.1, -89., 1.},

        {90., -45., 0.25}, {0., -45., 0.25},  {180., -45., 0.25}, {360., -45., 0.25}, {90., -90., 0.25},
        {90., 0., 0.25},   {90., 90., 0.25},  {10., -45., 0.25},  {20., -45., 0.25},  {30., -45., 0.25},
        {40., -45., 0.25}, {50., -45., 0.25}, {60., -45., 0.25},  {70., -45., 0.25},  {80., -45., 0.25},
        {90., -45., 0.25}, {10., -10., 0.25}, {20., -20., 0.25},  {30., -30., 0.25},  {40., -40., 0.25},
        {50., -50., 0.25}, {60., -60., 0.25}, {70., -70., 0.25},  {80., -80., 0.25},  {90., -90., 0.25},
        {60., -60., 0.6},  {90., -45., 0.16}, {90., -45., 0.6},   {90., -45., 1.},    {90., -45., 0.},
        {90., -45., 0.1},  {45., 33., 0.7},   {359., 20., 0.1},   {360., 88., 0.9},   {0.1, -89., 1.},

        {90., -45., 0.25}, {0., -45., 0.25},  {180., -45., 0.25}, {360., -45., 0.25}, {90., -90., 0.25},
        {90., 0., 0.25},   {90., 90., 0.25},  {10., -45., 0.25},  {20., -45., 0.25},  {30., -45., 0.25},
        {40., -45., 0.25}, {50., -45., 0.25}, {60., -45., 0.25},  {70., -45., 0.25},  {80., -45., 0.25},
        {90., -45., 0.25}, {10., -10., 0.25}, {20., -20., 0.25},  {30., -30., 0.25},  {40., -40., 0.25},
        {50., -50., 0.25}, {60., -60., 0.25}, {70., -70., 0.25},  {80., -80., 0.25},  {90., -90., 0.25},
        {60., -60., 0.6},  {90., -45., 0.16}, {90., -45., 0.6},   {90., -45., 1.},    {90., -45., 0.},
        {90., -45., 0.1},  {45., 33., 0.7},   {359., 20., 0.1},   {360., 88., 0.9},   {0.1, -89., 1.},
    });

    SECTION("prototype") {
        Cubic3DInterpolation cubic_interpolation(fs);

        for (auto p : departure_points.iterate().xyz()) {
            double interpolated = cubic_interpolation(p, f);
            double exact        = fp(p);
            Log::info() << p << "  -->  " << interpolated << "      [exact] " << exact << std::endl;
            EXPECT(is_approximately_equal(interpolated, exact, tolerance));
        }
    }


    SECTION("official version") {
        auto matrix_free = Config("matrix_free", true);
        Interpolation interpolation(option::type("structured-tricubic") | matrix_free, fs, departure_points);

        Field output = Field("output", make_datatype<double>(), make_shape(departure_points.size()));
        interpolation.execute(input, output);

        auto output_view = array::make_view<double, 1>(output);
        idx_t n{0};
        for (auto p : departure_points.iterate().xyz()) {
            double interpolated = output_view(n++);
            double exact        = fp(p);
            Log::info() << p << "  -->  " << interpolated << "      [exact] " << exact << std::endl;
            EXPECT(is_approximately_equal(interpolated, exact, tolerance));
        }
    }

    SECTION("SL-like") {
        auto matrix_free = Config("matrix_free", true);

        auto dp_field = fs.createField<double>(option::variables(3));

        {
            auto iterator     = departure_points.iterate().xyz().begin();
            auto iterator_end = departure_points.iterate().xyz().end();
            auto dp           = array::make_view<double, 3>(dp_field);
            for (idx_t n = 0; n < dp.shape(0); ++n) {
                for (idx_t k = 0; k < dp.shape(1); ++k) {
                    PointXYZ p{0, 0, 0};
                    if (iterator != iterator_end) {
                        p = *iterator;
                        ++iterator;
                    }
                    dp(n, k, LON) = p.x();
                    dp(n, k, LAT) = p.y();
                    dp(n, k, ZZ)  = p.z();
                }
            }
        }
        Interpolation interpolation(option::type("structured-tricubic") | matrix_free, fs, dp_field);

        Field output = fs.createField<double>();
        interpolation.execute(input, output);

        auto output_view = array::make_view<double, 2>(output);


        auto iterator     = departure_points.iterate().xyz().begin();
        auto iterator_end = departure_points.iterate().xyz().end();

        for (idx_t n = 0; n < output_view.shape(0); ++n) {
            for (idx_t k = 0; k < output_view.shape(1); ++k) {
                PointXYZ p{0, 0, 0};
                if (iterator != iterator_end) {
                    p = *iterator;
                    ++iterator;

                    double interpolated = output_view(n, k);
                    double exact        = fp(p);
                    Log::info() << p << "  -->  " << interpolated << "      [exact] " << exact << std::endl;
                    EXPECT(is_approximately_equal(interpolated, exact, tolerance));
                }
            }
        }
    }
}

}  // namespace test
}  // namespace atlas

int main(int argc, char** argv) {
    return atlas::test::run(argc, argv);
}
