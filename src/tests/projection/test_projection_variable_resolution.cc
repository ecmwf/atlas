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
#include "atlas/util/Point.h"
#include "atlas/util/Rotation.h"
#include "tests/AtlasTestEnvironment.h"

using namespace atlas::util;
using namespace atlas::grid;
using atlas::util::Rotation;

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

const int nx = 34;  // These number are here hardcoded, but there is a test using x/yrange_outer and delta_inner
const int ny = 32;
const double xrange_outer[2] = {348.13120344827576, 369.6190965517242};
const double yrange_outer[2] = {-8.264495551724226, 11.92110100000127};
const double xrange_inner[2] = {351.386944827586319, 366.363355172413776};
const double yrange_inner[2] = {-5.008754172413662, 8.665359620690706};

const double rim_width   = 4.;
const double delta_outer = 1.;
const double delta_inner = 0.6511482758621128;

///< correction used
constexpr float epstest = std::numeric_limits<float>::epsilon();


auto make_var_ratio_projection = [](double var_ratio) {
    Config conf;
    conf.set("type", "variable_resolution");
    conf.set("outer.dx", delta_outer);        ///< resolution of the external regular grid (rim)
    conf.set("inner.dx", delta_inner);        ///< resolution of the regional model (regular grid)
    conf.set("progression", var_ratio);       ///< power used for the stretching
    conf.set("inner.xmin", xrange_inner[0]);  ///< xstart of the internal regional grid
    conf.set("inner.ymin", yrange_inner[0]);  ///< ystart of the internal regional grid
    conf.set("inner.xend", xrange_inner[1]);  ///< xend of the regular part of stretched internal grid
    conf.set("inner.yend", yrange_inner[1]);  ///< yend of the regular part of stretched internal grid
    conf.set("outer.xmin", xrange_outer[0]);  ///< original domain startx
    conf.set("outer.xend", xrange_outer[1]);  ///< original domain endx
    conf.set("outer.ymin", yrange_outer[0]);  ///< original domain starty
    conf.set("outer.yend", yrange_outer[1]);  ///< original domain endy
    conf.set("outer.width", rim_width);
    return atlas::Projection(conf);
};

auto make_var_ratio_projection_rot = [](double var_ratio, std::vector<double> north_pole) {
    Config conf;
    conf.set("type", "rotated_variable_resolution");
    conf.set("outer.dx", delta_outer);        ///< resolution of the external regular grid (rim)
    conf.set("inner.dx", delta_inner);        ///< resolution of the regional model (regular grid)
    conf.set("progression", var_ratio);       ///< power used for the stretching
    conf.set("inner.xmin", xrange_inner[0]);  ///< xstart of the internal regional grid
    conf.set("inner.ymin", yrange_inner[0]);  ///< ystart of the internal regional grid
    conf.set("inner.xend", xrange_inner[1]);  ///< xend of the regular part of stretched internal grid
    conf.set("inner.yend", yrange_inner[1]);  ///< yend of the regular part of stretched internal grid
    conf.set("outer.xmin", xrange_outer[0]);  ///< original domain startx
    conf.set("outer.xend", xrange_outer[1]);  ///< original domain endx
    conf.set("outer.ymin", yrange_outer[0]);  ///< original domain starty
    conf.set("outer.yend", yrange_outer[1]);  ///< original domain endy
    conf.set("outer.width", rim_width);
    conf.set("north_pole", north_pole);
    return atlas::Projection(conf);
};


auto not_equal = [](double a, double b) { return std::abs(b - a) > 1.e-5; };

static double new_ratio(int n_stretched, double var_ratio) {
    /**
     *  compute ratio,
     *  change stretching factor so that high and low grids
     *  retain original sizes
     */


    ///< number of variable (stretched) grid points in one side
    double var_ints_f = n_stretched * 1.;
    double logr       = std::log(var_ratio);
    double log_ratio  = (var_ints_f - 0.5) * logr;

    return std::exp(log_ratio / n_stretched);
};

/**
 *  n_stretched and n_rim. Outside the regular grid, only one side,
 *  number of grid points in the stretched grid and in the rim region
*/
auto create_stretched_grid = [](const int& n_points_, const int& n_stretched_, const double& var_ratio_,
                                const int& n_rim, double lamphi, const bool& L_long) {
    double new_ratio_{};
    double range_outer[2]{};

    auto normalised = [L_long](double p) {
        if (L_long) {
            p = (p < 180) ? p + 360.0 : p;
        }
        return p;
    };

    lamphi = normalised(lamphi);

    if (L_long) {
        range_outer[0] = xrange_outer[0];
        range_outer[1] = xrange_outer[1];
    }
    else {
        range_outer[0] = yrange_outer[0];
        range_outer[1] = yrange_outer[1];
    }

    if (var_ratio_ > 1.) {
        //< compute number of points in the original regular grid
        // THIS IS TO NORMALIZE
        int n_idx_ = ((lamphi + epstest - range_outer[0]) / delta_inner) + 1;
        // first point in the internal regular grid
        int nstart_int = n_rim + n_stretched_ + 1;
        // last point in the internal regular grid
        int nend_int = n_points_ - (n_rim + n_stretched_);


        double delta      = delta_inner;
        double delta_last = delta_inner;
        double delta_add{};
        double delta_tot{};

        new_ratio_ = new_ratio(n_stretched_, var_ratio_);

        //< stretched area
        //< n_idx_ start from 1
        if (((n_idx_ < nstart_int) && (n_idx_ > n_rim)) || ((n_idx_ > nend_int) && (n_idx_ < n_points_ - n_rim + 1))) {
            //< number of stretched points for the loop
            int n_st{};
            if (n_idx_ < nstart_int) {
                n_st = nstart_int - n_idx_;
            }
            else {
                n_st = n_idx_ - nend_int;
            }

            delta_tot = 0.;
            for (int ix = 0; ix < n_st; ix += 1) {
                delta_last = delta * new_ratio_;
                delta_add  = delta_last - delta_inner;
                delta      = delta_last;
                delta_tot += delta_add;
            }

            if (n_idx_ < nstart_int) {
                lamphi -= delta_tot;
            }
            else {
                lamphi += delta_tot;
            }
        }

        //< rim region
        if (((n_idx_ < n_rim + 1)) || (n_idx_ > n_points_ - n_rim)) {
            delta_tot = 0.;
            //compute total stretched
            for (int ix = 0; ix < n_stretched_; ix += 1) {
                delta_last = delta * new_ratio_;
                delta_add  = delta_last - delta_inner;
                delta      = delta_last;
                delta_tot += delta_add;
            }


            double drim_size = ((rim_width / 2.) / n_rim) - delta_inner;
            int ndrim{};

            if (n_idx_ < nstart_int) {
                ndrim  = nstart_int - n_stretched_ - n_idx_;
                lamphi = lamphi - delta_tot - (ndrim * drim_size);
            }
            else {
                ndrim  = n_idx_ - (n_points_ - n_rim);
                lamphi = lamphi + delta_tot + ndrim * drim_size;
            }
        }
    }

    //< return the new point in the stretched grid
    return normalised(lamphi);
};

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

//< create stretched grid from regular grid

CASE("var_ratio_create = 1.13") {
    ///< definition of grid
    constexpr float epstest = std::numeric_limits<float>::epsilon();
    int nx_reg              = ((xrange_outer[1] - xrange_outer[0] + epstest) / delta_inner) + 1;
    int ny_reg              = ((yrange_outer[1] - yrange_outer[0] + epstest) / delta_inner) + 1;


    auto grid_st =
        RegularGrid{grid::LinearSpacing{xrange_outer[0], xrange_outer[1], nx_reg},
                    grid::LinearSpacing{yrange_outer[0], yrange_outer[1], ny_reg}, make_var_ratio_projection(1.13)};


    ///< check over regular grid points stretched using new atlas object and check using look-up table
    for (idx_t j = 0; j < grid_st.ny(); ++j) {
        for (idx_t i = 0; i < grid_st.nx(); ++i) {
            auto ll = grid_st.lonlat(i, j);
            EXPECT_APPROX_EQ(ll.lon(), lon_LAM_str[i], 1.e-10);
            EXPECT_APPROX_EQ(ll.lat(), lat_LAM_str[j], 1.e-10);
        }
    }


    ///< check over regular grid points stretched using new atlas object and check using function above
    for (idx_t j = 0; j < grid_st.ny(); ++j) {
        for (idx_t i = 0; i < grid_st.nx(); ++i) {
            auto ll_st_lon = create_stretched_grid(nx_reg, 3, 1.13, 2, lon_LAM_reg[i], true);
            auto ll_st_lat = create_stretched_grid(ny_reg, 3, 1.13, 2, lat_LAM_reg[j], false);
            EXPECT_APPROX_EQ(ll_st_lon, lon_LAM_str[i], 1.e-10);
            EXPECT_APPROX_EQ(ll_st_lat, lat_LAM_str[j], 1.e-10);
        }
    }

    ///< check internal regular grid
    idx_t ymid             = grid_st.ny() / 2;
    idx_t xmid             = grid_st.nx() / 2;
    auto expect_equal_dlon = [&](int i, double dlon) {
        EXPECT_APPROX_EQ(grid_st.lonlat(i + 1, ymid).lon() - grid_st.lonlat(i, ymid).lon(), dlon, 1.e-10);
    };
    auto expect_equal_dlat = [&](int j, double dlat) {
        EXPECT_APPROX_EQ(grid_st.lonlat(xmid, j + 1).lat() - grid_st.lonlat(xmid, j).lat(), dlat, 1.e-10);
    };
    expect_equal_dlon(0, delta_outer);
    expect_equal_dlon(xmid, delta_inner);
    expect_equal_dlat(0, delta_outer);
    expect_equal_dlat(ymid, delta_inner);

    //< Check that the spacing in xy coordinates matches "delta_inner"
    for (int i = 0; i < grid_st.nx() - 1; ++i) {
        EXPECT_APPROX_EQ(grid_st.xy(i + 1, ymid).x() - grid_st.xy(i, ymid).x(), delta_inner, 1.e-10);
    }
    for (int j = 0; j < grid_st.ny() - 1; ++j) {
        EXPECT_APPROX_EQ(grid_st.xy(xmid, j + 1).y() - grid_st.xy(xmid, j).y(), delta_inner, 1.e-10);
    }
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

    //< Check that the spacing in xy coordinates matches "delta_inner"
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


CASE("var_ratio_rot = 1.13") {
    ///< check over regular grid points stretched using new atlas object and check using look-up table
    ///< in this case add a rotation
    Config config;
    std::vector<double> north_pole = {-176., 40.};
    config.set("north_pole", north_pole);

    ///< definition of grid I have to rotate this
    auto proj_st = make_var_ratio_projection_rot(1.13, north_pole);
    auto grid    = RegularGrid{grid::LinearSpacing{xrange_outer[0], xrange_outer[1], nx},
                            grid::LinearSpacing{yrange_outer[0], yrange_outer[1], ny}, proj_st};

    Rotation rotation(config);

    for (idx_t j = 0; j < grid.ny(); ++j) {
        for (idx_t i = 0; i < grid.nx(); ++i) {
            ///<compare rotated stretched grid ll, with the defined array to rotate
            auto ll     = grid.lonlat(i, j);
            auto ll_str = {lon_LAM_str[i], lat_LAM_str[j]};
            EXPECT_APPROX_EQ(ll, rotation.rotate(ll_str), 1.e-8);
        }
    }

    idx_t ymid = grid.ny() / 2;
    idx_t xmid = grid.nx() / 2;

    //< Check that the spacing in xy coordinates matches "delta_inner"
    for (int i = 0; i < grid.nx() - 1; ++i) {
        EXPECT_APPROX_EQ(grid.xy(i + 1, ymid).x() - grid.xy(i, ymid).x(), delta_inner, 1.e-8);
    }
    for (int j = 0; j < grid.ny() - 1; ++j) {
        EXPECT_APPROX_EQ(grid.xy(xmid, j + 1).y() - grid.xy(xmid, j).y(), delta_inner, 1.e-8);
    }
    ///< Set meshes
    ///< define mesh to write in file
    auto mesh = Mesh{grid};

    /**
     *  Write mesh in gmsh object.
     *  output under:
     *   <build-directory>/atlas/src/tests/grid/
     */

    output::Gmsh{"stretch_mesh_lonlat_rot.msh", Config("coordinates", "lonlat")("info", true)}.write(mesh);
    output::Gmsh{"stretch_mesh_xy_rot.msh", Config("coordinates", "xy")("info", true)}.write(mesh);


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

    ///
    ///< create regular grid
    auto grid_reg_approx =
        RegularGrid{LinearSpacing{startx_n, endx_n, nx_new}, LinearSpacing{starty_n, endy_n, ny_new}};

    /**
     *  Write mesh in gmsh object.
     *  output under:
     *   <build-directory>/atlas/src/tests/grid/
     */
    output::Gmsh("grid_regrot_approx_lonlat.msh", Config("coordinates", "lonlat")("info", true))
        .write(Mesh{grid_reg_approx});
}

CASE("var_ratio_rot_inv = 1.13") {
    ///< check over regular grid points stretched using new atlas object and check using look-up table
    ///< in this case add a rotation
    Config config;
    ///< definition of grid I have to rotate this
    auto proj_st_nr = make_var_ratio_projection(1.13);
    auto grid_nr    = RegularGrid{grid::LinearSpacing{xrange_outer[0], xrange_outer[1], nx},
                               grid::LinearSpacing{yrange_outer[0], yrange_outer[1], ny}, proj_st_nr};

    ///< definition of grid I have to rotate this
    auto proj_st = make_var_ratio_projection_rot(1.13, {-176., 40.});
    auto grid    = RegularGrid{grid::LinearSpacing{xrange_outer[0], xrange_outer[1], nx},
                            grid::LinearSpacing{yrange_outer[0], yrange_outer[1], ny}, proj_st};


    ///< Test the inverse
    for (idx_t j = 0; j < grid.ny(); ++j) {
        for (idx_t i = 0; i < grid.nx(); ++i) {
            ///<compare rotated stretched grid ll, with the defined array to rotate
            ///< rotated
            Point2 ll_str = grid.lonlat(i, j);
            ///< unrotated
            Point2 ll_str1 = {lon_LAM_str[i], lat_LAM_str[j]};
            //< from stretched rotated to unrotated regular
            grid.projection().lonlat2xy(ll_str);
            Point2 ll_reg = {lon_LAM_reg[i], lat_LAM_reg[j]};
            EXPECT_APPROX_EQ(ll_str, ll_reg, 1.e-8);
        }
    }
}


}  // namespace test
}  // namespace atlas


int main(int argc, char* argv[]) {
    return atlas::test::run(argc, argv);
}
