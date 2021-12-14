/**
* (C) Crown copyright 2021, Met Office
*
* This software is licensed under the terms of the Apache Licence Version 2.0
* which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
*/


#include "atlas/grid.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/option.h"
#include "atlas/output/Gmsh.h"
#include "atlas/projection.h"
#include "atlas/util/Config.h"
#include "tests/AtlasTestEnvironment.h"

using namespace atlas::util;
using namespace atlas::grid;

namespace {

///< use vectors of double for testing

const std::vector<double> lon_LAM_str{
    346.9838795150958504, 347.9838795150958504, 348.9838795150958504, 349.8677242998691668, 350.6659835640296592,
    351.3869448275862055, 352.0380931034482614, 352.6892413793103174, 353.3403896551723733, 353.9915379310344861,
    354.642686206896542,  355.2938344827585979, 355.9449827586206538, 356.5961310344827666, 357.2472793103448225,
    357.8984275862068785, 358.5495758620689344, 359.2007241379310472, 359.8518724137931031, 360.503020689655159,
    361.1541689655172149, 361.8053172413793277, 362.4564655172413836, 363.1076137931034395, 363.7587620689654955,
    364.4099103448276082, 365.0610586206896642, 365.7122068965517201, 366.363355172413776,  367.0843164359702087,
    367.8825757001305305, 368.7664204849036764, 369.7664204849036764, 370.7664204849036764};

const std::vector<double> lat_LAM_str{
    -9.41181948490414122,  -8.41181948490414122,  -7.41181948490414122,   -6.527974700130756425, -5.729715435970231141,
    -5.00875417241366172,  -4.357605896551548952, -3.706457620689436183,  -3.055309344827323415, -2.404161068965210646,
    -1.753012793103097877, -1.101864517240985109, -0.4507162413788723399, 0.2004320344832404288, 0.8515803103453531975,
    1.502728586207465966,  2.153876862069578735,  2.805025137931691503,   3.456173413793804272,  4.107321689655917041,
    4.758469965518029809,  5.409618241380142578,  6.060766517242255347,   6.711914793104368115,  7.363063068966480884,
    8.014211344828593653,  8.665359620690706421,  9.386320884247275842,   10.18458014840780024,  11.06842493318118414,
    12.06842493318118414,  13.06842493318118414};


const std::vector<double> lon_LAM_reg = {
    348.13120344827576,   348.78235172413787,   349.4334999999999809, 350.0846482758620368, 350.7357965517240928,
    351.3869448275862055, 352.0380931034482614, 352.6892413793103174, 353.3403896551723733, 353.9915379310344861,
    354.642686206896542,  355.2938344827585979, 355.9449827586206538, 356.5961310344827666, 357.2472793103448225,
    357.8984275862068785, 358.5495758620689344, 359.2007241379310472, 359.8518724137931031, 360.503020689655159,
    361.1541689655172149, 361.8053172413793277, 362.4564655172413836, 363.1076137931034395, 363.7587620689654955,
    364.4099103448276082, 365.0610586206896642, 365.7122068965517201, 366.363355172413776,  367.0145034482758888,
    367.6656517241379447, 368.3168000000000006, 368.9679482758621,    369.6190965517242};

const std::vector<double> lat_LAM_reg = {-8.264495551724226,     -7.613347275862113,    -6.962199,
                                         -6.311050724137887258,  -5.659902448275774489, -5.00875417241366172,
                                         -4.357605896551548952,  -3.706457620689436183, -3.055309344827323415,
                                         -2.404161068965210646,  -1.753012793103097877, -1.101864517240985109,
                                         -0.4507162413788723399, 0.2004320344832404288, 0.8515803103453531975,
                                         1.502728586207465966,   2.153876862069578735,  2.805025137931691503,
                                         3.456173413793804272,   4.107321689655917041,  4.758469965518029809,
                                         5.409618241380142578,   6.060766517242255347,  6.711914793104368115,
                                         7.363063068966480884,   8.014211344828593653,  8.665359620690706421,
                                         9.31650789655281919,    9.967656172414931959,  10.61880444827704473,
                                         11.269952724139157,     11.92110100000127};

const int nx                 = 34;  // How do you find this?
const int ny                 = 32;
const double xrange_outer[2] = {348.13120344827576, 369.6190965517242};
const double yrange_outer[2] = {-8.264495551724226, 11.92110100000127};
const double xrange_inner[2] = {351.386944827586319, 366.363355172413776};
const double yrange_inner[2] = {-5.008754172413662, 8.665359620690706};

const double rim_width   = 4.;
const double delta_outer = 1.;
const double delta_inner = 0.6511482758621128;

auto make_var_ratio_projection = [](double var_ratio) {
    Config conf;
    conf.set("type", "variable_resolution");
    conf.set("var_ratio", var_ratio);
    conf.set("delta_low", delta_outer);
    conf.set("delta_hi", delta_inner);
    conf.set("x_reg_start", xrange_inner[0]);
    conf.set("x_reg_end", xrange_inner[1]);
    conf.set("y_reg_start", yrange_inner[0]);
    conf.set("y_reg_end", yrange_inner[1]);
    conf.set("startx", xrange_outer[0]);
    conf.set("endx", xrange_outer[1]);
    conf.set("starty", yrange_outer[0]);
    conf.set("endy", yrange_outer[1]);
    conf.set("north_pole", {0.0, 90.0});
    conf.set("rim_widthx", rim_width);
    conf.set("rim_widthy", rim_width);
    return atlas::Projection(conf);
};

auto not_equal = [](double a, double b) { return std::abs(b - a) > 1.e-5; };

};  // namespace


namespace atlas {
namespace test {

CASE("Understanding of the above data") {
    const std::vector<double>& v = lon_LAM_str;
    std::vector<double> delta(v.size() - 1);
    for (size_t i = 0; i < delta.size(); ++i) {
        delta[i] = v[i + 1] - v[i];
    }
    Log::info() << "delta = " << delta << std::endl;
    // Outputs:
    // delta = [1,1,0.883845,0.798259,0.720961,0.651148,0.651148,0.651148,0.651148,
    //          0.651148,0.651148,0.651148,0.651148,0.651148,0.651148,0.651148,0.651148,
    //          0.651148,0.651148,0.651148,0.651148,0.651148,0.651148,0.651148,0.651148,
    //          0.651148,0.651148,0.651148,0.720961,0.798259,0.883845,1,1]

    std::vector<double> delta_half(delta.begin() + delta.size() / 2, delta.end());
    Log::info() << "delta_half = " << delta_half << std::endl;
    // Outputs:
    // delta_half = [0.651148,0.651148,0.651148,0.651148,0.651148,0.651148,0.651148,
    //               0.651148,0.651148,0.651148,0.651148,0.651148,0.720961,0.798259,
    //               0.883845,1,1]
    std::vector<double> progression(delta_half.size() - 1);
    for (size_t i = 0; i < progression.size(); ++i) {
        progression[i] = delta_half[i + 1] / delta_half[i];
    }
    Log::info() << "progression = " << progression << std::endl;
    // Outputs:
    // progression = [1,1,1,1,1,1,1,1,1,1,1,1.10722,1.10722,1.10722,1.13142,1]

    int nx           = 34;
    int nx_var       = 6;
    int nx_outer     = 4;
    double var_ratio = 1.13;

    constexpr double epstest =
        std::numeric_limits<float>::epsilon();  ///< correction used to change from double to integer

    int inner_i_begin = -1;
    int inner_i_end   = -1;
    for (int i = 0; i < nx; ++i) {
        if (std::abs(v[i] - xrange_inner[0]) < 1.e-10) {
            inner_i_begin = i;
        }
        if (std::abs(v[i] - xrange_inner[1]) < 1.e-10) {
            inner_i_end = i + 1;
        }
    }

    double inner_size = ((nx - 1) - nx_var - nx_outer) * delta_inner;
    ///< number of regular internal grid points, integer
    int nx_inner = (inner_size + epstest) / delta_inner + 1;

    EXPECT_EQ(nx_inner, inner_i_end - inner_i_begin);
    EXPECT_EQ(nx_inner, inner_i_end - inner_i_begin);
    EXPECT_EQ(inner_i_begin, nx_outer / 2 + nx_var / 2);
    EXPECT_EQ(inner_i_end, nx - nx_outer / 2 - nx_var / 2);
    for (int i = 0; i < nx_outer / 2; ++i) {
        EXPECT_APPROX_EQ(v[i + 1] - v[i], delta_outer, 1.e-10);
    }
    for (int i = nx_outer / 2; i < inner_i_begin; ++i) {
        double d = v[i + 1] - v[i];
        EXPECT(not_equal(d, delta_inner));
        EXPECT(not_equal(d, delta_outer));
    }
    for (int i = inner_i_begin; i < inner_i_end - 1; ++i) {
        EXPECT_APPROX_EQ(v[i + 1] - v[i], delta_inner, 1.e-10);
    }
    for (int i = inner_i_end - 1; i < nx - nx_outer / 2 - 1; ++i) {
        double d = v[i + 1] - v[i];
        EXPECT(not_equal(d, delta_inner));
        EXPECT(not_equal(d, delta_outer));
    }
    for (int i = nx - nx_outer / 2 - 1; i < nx - 1; ++i) {
        EXPECT_APPROX_EQ(v[i + 1] - v[i], delta_outer, 1.e-10);
    }

    double var_ints_f = ((nx - 1) - nx_outer - (nx_inner - 1)) / 2.;
    double logr       = std::log(var_ratio);
    double log_ratio  = (var_ints_f - 0.5) * logr;
    double new_ratio  = std::exp(log_ratio / std::floor(var_ints_f));

    EXPECT_APPROX_EQ(new_ratio, 1.10722, 1.e-5);
}


CASE("var_ratio = 1.13") {
    ///< definition of grid
    auto grid = RegularGrid{grid::LinearSpacing{xrange_outer[0], xrange_outer[1], nx},
                            grid::LinearSpacing{yrange_outer[0], yrange_outer[1], ny}, make_var_ratio_projection(1.13)};

    ///< check over regular grid points stretched using new atlas object and check using look-up table
    for (idx_t j = 0; j < grid.ny(); ++j) {
        for (idx_t i = 0; i < grid.nx(); ++i) {
            auto ll = grid.lonlat(i, j);
            EXPECT_APPROX_EQ(ll.lon(), lon_LAM_str[i], 1.e-10);
            EXPECT_APPROX_EQ(ll.lat(), lat_LAM_str[j], 1.e-10);
        }
    }

    idx_t ymid             = grid.ny() / 2;
    idx_t xmid             = grid.nx() / 2;
    auto expect_equal_dlon = [&](int i, double dlon) {
        EXPECT_APPROX_EQ(grid.lonlat(i + 1, ymid).lon() - grid.lonlat(i, ymid).lon(), dlon, 1.e-10);
    };
    auto expect_equal_dlat = [&](int j, double dlat) {
        EXPECT_APPROX_EQ(grid.lonlat(xmid, j + 1).lat() - grid.lonlat(xmid, j).lat(), dlat, 1.e-10);
    };
    expect_equal_dlon(0, delta_outer);
    expect_equal_dlon(xmid, delta_inner);
    expect_equal_dlat(0, delta_outer);
    expect_equal_dlat(ymid, delta_inner);

    // Check that the spacing in xy coordinates matches "delta_inner"
    for (int i = 0; i < grid.nx() - 1; ++i) {
        EXPECT_APPROX_EQ(grid.xy(i + 1, ymid).x() - grid.xy(i, ymid).x(), delta_inner, 1.e-10);
    }
    for (int j = 0; j < grid.ny() - 1; ++j) {
        EXPECT_APPROX_EQ(grid.xy(xmid, j + 1).y() - grid.xy(xmid, j).y(), delta_inner, 1.e-10);
    }


    ///< Set meshes
    ///< define mesh to write in file
    auto mesh = Mesh{grid};

    /**
     *  Write mesh in gmsh object.
     *  output under:
     *   <build-directory>/atlas/src/tests/grid/
     */

    output::Gmsh{"stretch_mesh_lonlat.msh", Config("coordinates", "lonlat")("info", true)}.write(mesh);
    output::Gmsh{"stretch_mesh_xy.msh", Config("coordinates", "xy")("info", true)}.write(mesh);


    /** Create additional regular grid with same delta_hi
     *  for additional .msh file with lon lat in approximately the same range
     *  as the stretched one.
     */

    int nx_new = 37;
    int ny_new = 35;

    ///< 37 the exact range would have 36.52 points
    double alpha    = ((delta_inner * nx_new) - (xrange_outer[1] - xrange_outer[0])) / 2.;
    double startx_n = xrange_outer[0] - alpha;
    double endx_n   = xrange_outer[1] + alpha;
    double starty_n = yrange_outer[0] - alpha;
    double endy_n   = yrange_outer[1] + alpha;


    ///< create regular grid
    auto grid_reg_approx =
        RegularGrid{LinearSpacing{startx_n, endx_n, nx_new}, LinearSpacing{starty_n, endy_n, ny_new}};


    /**
     *  Write mesh in gmsh object.
     *  output under:
     *   <build-directory>/atlas/src/tests/grid/
     */
    output::Gmsh("grid_reg_approx_lonlat.msh", Config("coordinates", "lonlat")("info", true))
        .write(Mesh{grid_reg_approx});
}

CASE("var_ratio = 1.0") {
    ///< TEST of var_ratio = 1.0 configuration
    ///< definition of stretched grid
    auto grid = RegularGrid{grid::LinearSpacing{xrange_outer[0], xrange_outer[1], nx},
                            grid::LinearSpacing{yrange_outer[0], yrange_outer[1], ny}, make_var_ratio_projection(1.0)};

    ///< Check if bounding box is correct
    {
        auto bb = grid.lonlatBoundingBox();
        EXPECT(RectangularLonLatDomain(bb));
        const double tolerance = 1.e-6;
        EXPECT_APPROX_EQ(bb.west(), xrange_outer[0], tolerance);
        EXPECT_APPROX_EQ(bb.east(), xrange_outer[1], tolerance);
        EXPECT_APPROX_EQ(bb.south(), yrange_outer[0], tolerance);
        EXPECT_APPROX_EQ(bb.north(), yrange_outer[1], tolerance);
        for (PointLonLat p : grid.lonlat()) {
            EXPECT(bb.contains(p));
        }
    }

    ///< check over regular grid points stretched using new atlas object and check using look-up table
    for (idx_t j = 0; j < grid.ny(); ++j) {
        for (idx_t i = 0; i < grid.nx(); ++i) {
            auto ll = grid.lonlat(i, j);
            EXPECT_APPROX_EQ(ll.lon(), lon_LAM_reg[i], 1.e-10);
            EXPECT_APPROX_EQ(ll.lat(), lat_LAM_reg[j], 1.e-10);
        }
    }

    auto mesh = Mesh{grid};
    output::Gmsh{"reg_mesh_lonlat.msh", Config("coordinates", "lonlat")("info", true)}.write(mesh);
    output::Gmsh{"reg_mesh_xy.msh", Config("coordinates", "xy")("info", true)}.write(mesh);
}

}  // namespace test
}  // namespace atlas


int main(int argc, char* argv[]) {
    return atlas::test::run(argc, argv);
}
