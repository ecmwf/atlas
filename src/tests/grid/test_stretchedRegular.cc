/**
* (C) Crown copyright 2021, Met Office
*
* This software is licensed under the terms of the Apache Licence Version 2.0
* which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
*/


#include <algorithm>
#include <vector>

#include "atlas/array.h"
#include "atlas/grid.h"
#include "atlas/option.h"
#include "atlas/util/Config.h"
#include "atlas/functionspace.h"
#include "atlas/functionspace/StructuredColumns.h"
#include "atlas/projection.h"
#include "atlas/projection/detail/StretchRegular.h"
#include "tests/AtlasTestEnvironment.h"
#include "atlas/output/Gmsh.h"
#include "atlas/meshgenerator.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/mesh/Mesh.h"

using namespace atlas::util;
using namespace atlas::grid;

namespace {

///< use vectors of double for testing

const double lon_LAM_str[] = {
                                346.9838795150958504,  347.9838795150958504,
                                348.9838795150958504,
                                349.8677242998691668,
                                350.6659835640296592,
                                351.3869448275862055,
                                352.0380931034482614,
                                352.6892413793103174,
                                353.3403896551723733,
                                353.9915379310344861,
                                354.642686206896542,
                                355.2938344827585979,
                                355.9449827586206538,
                                356.5961310344827666,
                                357.2472793103448225,
                                357.8984275862068785,
                                358.5495758620689344,
                                359.2007241379310472,
                                359.8518724137931031,
                                360.503020689655159,
                                361.1541689655172149,
                                361.8053172413793277,
                                362.4564655172413836,
                                363.1076137931034395,
                                363.7587620689654955,
                                364.4099103448276082,
                                365.0610586206896642,
                                365.7122068965517201,
                                366.363355172413776,
                                367.0843164359702087,
                                367.8825757001305305,
                                368.7664204849036764,
                                369.7664204849036764,  370.7664204849036764
    };

const double lat_LAM_str[] = {  -9.41181948490414122,  -8.41181948490414122,
                                -7.41181948490414122,
                                -6.527974700130756425,
                                -5.729715435970231141,
                                -5.00875417241366172,
                                -4.357605896551548952,
                                -3.706457620689436183,
                                -3.055309344827323415,
                                -2.404161068965210646,
                                -1.753012793103097877,
                                -1.101864517240985109,
                                -0.4507162413788723399,
                                0.2004320344832404288,
                                0.8515803103453531975,
                                1.502728586207465966,
                                2.153876862069578735,
                                2.805025137931691503,
                                3.456173413793804272,
                                4.107321689655917041,
                                4.758469965518029809,
                                5.409618241380142578,
                                6.060766517242255347,
                                6.711914793104368115,
                                7.363063068966480884,
                                8.014211344828593653,
                                8.665359620690706421,
                                9.386320884247275842,
                                10.18458014840780024,
                                11.06842493318118414,
                                12.06842493318118414, 13.06842493318118414
    };


const double lon_LAM_reg[] = {348.13120344827576, 348.78235172413787,
                              349.4334999999999809,
                              350.0846482758620368,
                              350.7357965517240928,
                              351.3869448275862055,
                              352.0380931034482614,
                              352.6892413793103174,
                              353.3403896551723733,
                              353.9915379310344861,
                              354.642686206896542,
                              355.2938344827585979,
                              355.9449827586206538,
                              356.5961310344827666,
                              357.2472793103448225,
                              357.8984275862068785,
                              358.5495758620689344,
                              359.2007241379310472,
                              359.8518724137931031,
                              360.503020689655159,
                              361.1541689655172149,
                              361.8053172413793277,
                              362.4564655172413836,
                              363.1076137931034395,
                              363.7587620689654955,
                              364.4099103448276082,
                              365.0610586206896642,
                              365.7122068965517201,
                              366.363355172413776,
                              367.0145034482758888,
                              367.6656517241379447,
                              368.3168000000000006,
                              368.9679482758621, 369.6190965517242};

const double lat_LAM_reg[] = { -8.264495551724226, -7.613347275862113,
                               -6.962199,
                               -6.311050724137887258,
                               -5.659902448275774489,
                               -5.00875417241366172,
                               -4.357605896551548952,
                               -3.706457620689436183,
                               -3.055309344827323415,
                               -2.404161068965210646,
                               -1.753012793103097877,
                               -1.101864517240985109,
                               -0.4507162413788723399,
                               0.2004320344832404288,
                               0.8515803103453531975,
                               1.502728586207465966,
                               2.153876862069578735,
                               2.805025137931691503,
                               3.456173413793804272,
                               4.107321689655917041,
                               4.758469965518029809,
                               5.409618241380142578,
                               6.060766517242255347,
                               6.711914793104368115,
                               7.363063068966480884,
                               8.014211344828593653,
                               8.665359620690706421,
                               9.31650789655281919,
                               9.967656172414931959,
                               10.61880444827704473,
                              11.269952724139157, 11.92110100000127
                              };

};





namespace atlas {
namespace test {
CASE( "Regularstretch" ) {

    auto proj = Projection( "stretch", Config( "delta_low", 1. ) | Config( "delta_hi", 0.6511482758621128 ) |
                                                 Config( "var_ratio", 1.13 ) | Config( "x_reg_start", 351.386944827586319 ) |
                                                 Config( "y_reg_start", -5.008754172413662 ) | Config( "x_reg_end", 366.363355172413776 ) |
                                                 Config( "y_reg_end", 8.665359620690706 ) | Config( "startx", 348.13120344827576) |
                                                 Config( "starty",  -8.264495551724226) | Config( "endx", 369.6190965517242) | Config( "endy", 11.92110100000127) |
                                                 Config( "north_pole", {0.0, 90.0} ) | Config( "rim_widthx", 4. ) | Config( "rim_widthy", 4. ) );





    atlas::util::Config XSpaceConfig;
    XSpaceConfig.set( "type", "linear" );
    XSpaceConfig.set( "N", 34 );
    XSpaceConfig.set("start", 348.13120344827576 );
    XSpaceConfig.set("end", 369.6190965517242 );

    atlas::grid::detail::grid::Structured::XSpace XS(XSpaceConfig);

    atlas::util::Config YSpaceConfig;
    YSpaceConfig.set( "type", "linear" );
    YSpaceConfig.set( "N", 32 );
    YSpaceConfig.set("start",  -8.264495551724226 );
    YSpaceConfig.set("end", 11.92110100000127 );

    atlas::grid::detail::grid::Structured::YSpace YS(YSpaceConfig);

    ///< definition of stretched grid
    auto grid_st = StructuredGrid(XS, YS, proj );

    ///< create regular grid
    atlas::util::Config reg_grid_config;
    reg_grid_config.set( "type", "structured" );
    reg_grid_config.set ("xspace", []() {
      atlas::util::Config config;
      config.set( "type", "linear" );
      config.set( "N", 34 );
      config.set( "start", 348.13120344827576 );
      config.set( "end", 369.6190965517242 );
      return config;
    }() );

    reg_grid_config.set ("yspace", []() {
      atlas::util::Config config;
      config.set( "type", "linear" );
      config.set( "N", 32 );
      config.set( "start", -8.264495551724226 );
      config.set( "end", 11.92110100000127);
      return config;
    }() );

    atlas::StructuredGrid reg_grid(reg_grid_config);

    atlas::functionspace::StructuredColumns nodes_reg(
           reg_grid, atlas::grid::Partitioner( "equal_regions" ), reg_grid_config );

    int sizei = nodes_reg.i_end(0) - nodes_reg.i_begin(0);
    int sizej = nodes_reg.j_end() - nodes_reg.j_begin();

    std::vector<double> lon_p_arr_st(sizei);
    std::vector<double> lon_p_arr(sizei);
    std::vector<double> lat_p_arr_st(sizej);
    std::vector<double> lat_p_arr(sizej);

    ///< check over regular grid points stretched using new atlas object and check using look-up table
    for (atlas::idx_t j = nodes_reg.j_begin(); j < nodes_reg.j_end(); ++j) {
        for (atlas::idx_t i = nodes_reg.i_begin(j); i < nodes_reg.i_end(j); ++i) {
            auto ll1lon = lon_LAM_str[i];
            auto ll1lat = lat_LAM_str[j];
            auto ll2 = grid_st.lonlat( i, j );
            auto ll2lon = ll2.lon();
            auto ll2lat = ll2.lat();
            EXPECT_APPROX_EQ( ll1lon , ll2lon, 1.e-10 );
            EXPECT_APPROX_EQ( ll1lat, ll2lat, 1.e-10 );
        }
    }

    ///< Set meshes
    ///< define mesh to write in file
    auto meshGen_st = atlas::MeshGenerator( "structured" );
    Mesh mesh_st = StructuredMeshGenerator().generate( grid_st );
    auto meshGen_reg = atlas::MeshGenerator( "regular" );
    auto mesh_reg    = meshGen_reg.generate( reg_grid );

    ///< Set gmsh config.
    ///< ghost is the haloe, in a general container
    auto gmshConfig_reg = atlas::util::Config( "coordinates", "xy" ) | atlas::util::Config( "ghost", false );
    gmshConfig_reg.set( "info", true );

    auto gmshConfig_stretch = atlas::util::Config( "coordinates", "lonlat" ) | atlas::util::Config( "ghost", false );
    gmshConfig_stretch.set( "info", true );

    ///< Set source gmsh object.
    const auto gmsh_reg     = atlas::output::Gmsh( "xy_mesh.msh", gmshConfig_reg );
    output::Gmsh gmsh_stretch( "stretch_mesh.msh" , gmshConfig_stretch );

    /**
     *  Write mesh in gmsh object.
     *  output under:
     *   <build-directory>/atlas/src/tests/grid/
     */


    gmsh_reg.write( mesh_reg );
    gmsh_stretch.write( mesh_st );

    /** Create additional regular grid
     *  for additional .msh file with lon lat in approximately the same range
     *  as the stretched one.
     */

    double startx = 348.13120344827576;
    double endx = 369.6190965517242;
    double starty = -8.264495551724226 ;
    double endy = 11.92110100000127;
    double alpha, startx_n, endx_n, starty_n, endy_n;
    int nx_new = 37, ny_new = 35;
    double h_res = 0.6511482758621128;

    ///< 37 the exact range would have 36.52 points
    alpha = ((h_res * nx_new) - (endx - startx))/2 ;
    startx_n = startx - alpha;
    endx_n = endx + alpha;
    starty_n = starty - alpha;
    endy_n = endy + alpha;


    ///< create regular grid
    atlas::util::Config reg_grid_config_ex;
    reg_grid_config_ex.set( "type", "structured" );

    reg_grid_config_ex.set ("xspace", [startx_n, endx_n, nx_new]() {
      atlas::util::Config config;
      config.set( "type", "linear" );
      config.set( "N", nx_new );
      config.set( "start", startx_n);
      config.set( "end", endx_n );
      return config;
    }() );

    reg_grid_config_ex.set ("yspace", [starty_n, endy_n, ny_new]() {
      atlas::util::Config config;
      config.set( "type", "linear" );
      config.set( "N", ny_new );
      config.set( "start", starty_n );
      config.set( "end", endy_n);
      return config;
    }() );

    atlas::StructuredGrid reg_grid_ex(reg_grid_config_ex);

    ///< Set meshes
    auto mesh_reg_ex    = meshGen_reg.generate( reg_grid_ex );

    ///< Set gmsh config.
    ///< ghost is the haloe, in a general container
    auto gmshConfig_reg_ex = atlas::util::Config( "coordinates", "xy" ) | atlas::util::Config( "ghost", false );
    gmshConfig_reg.set( "info", true );

    ///< Set source gmsh object.
    const auto gmsh_reg_ex     = atlas::output::Gmsh( "reg_mesh.msh", gmshConfig_reg );

    /**
     *  Write mesh in gmsh object.
     *  output under:
     *   <build-directory>/atlas/src/tests/grid/
     */
    gmsh_reg_ex.write( mesh_reg_ex );



    ///< TEST of var_ratio = 1.0 configuration

    auto proj_reg = Projection( "stretch", Config( "delta_low", 1. ) | Config( "delta_hi", 0.6511482758621128 ) |
                                                 Config( "var_ratio", 1.0 ) | Config( "x_reg_start", 351.386944827586319 ) |
                                                 Config( "y_reg_start", -5.008754172413662 ) | Config( "x_reg_end", 366.363355172413776 ) |
                                                 Config( "y_reg_end", 8.665359620690706 ) | Config( "startx", 348.13120344827576) |
                                                 Config( "starty",  -8.264495551724226) | Config( "endx", 369.6190965517242) | Config( "endy", 11.92110100000127) |
                                                 Config( "north_pole", {0.0, 90.0} ) | Config( "rim_widthx", 4. ) | Config( "rim_widthy", 4. ) );



    ///< definition of stretched grid
    auto grid_reg = StructuredGrid(XS, YS, proj_reg );
    ///< Check if bounding box is correct
        {
            RectangularLonLatDomain bb{grid_reg.lonlatBoundingBox()};
            EXPECT( RectangularLonLatDomain( bb ) );
            const double tolerance = 1.e-6;
            EXPECT_APPROX_EQ( bb.west(), 348.13120344827576, tolerance );
            EXPECT_APPROX_EQ( bb.east(), 369.6190965517242, tolerance );
            EXPECT_APPROX_EQ( bb.south(), -8.264495551724226, tolerance );
            EXPECT_APPROX_EQ( bb.north(), 11.92110100000127, tolerance );
            for ( PointLonLat p : grid_reg.lonlat() ) {
                EXPECT( bb.contains( p ) );
            }
          }




    ///< check over regular grid points stretched using new atlas object and check using look-up table
    for (atlas::idx_t j = nodes_reg.j_begin(); j < nodes_reg.j_end(); ++j) {
        for (atlas::idx_t i = nodes_reg.i_begin(j); i < nodes_reg.i_end(j); ++i) {
            auto ll1lon = lon_LAM_reg[i];
            auto ll1lat = lat_LAM_reg[j];
            auto ll2 = grid_reg.lonlat( i, j );
            auto ll2lon = ll2.lon();
            auto ll2lat = ll2.lat();
            EXPECT_APPROX_EQ( ll1lon , ll2lon, 1.e-10 );
            EXPECT_APPROX_EQ( ll1lat, ll2lat, 1.e-10 );
        }
    }

}




}  ///< namespace test
}  ///< namespace atlas


int main( int argc, char* argv[] ) {
    return atlas::test::run( argc, argv );
}
