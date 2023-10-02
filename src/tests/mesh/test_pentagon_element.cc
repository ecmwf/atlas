/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "atlas/functionspace.h"
#include "atlas/interpolation.h"
#include "atlas/mesh.h"
#include "atlas/option.h"
#include "atlas/output/Gmsh.h"
#include "atlas/output/Output.h"

#include "tests/AtlasTestEnvironment.h"
#include "tests/TestMeshes.h"

namespace atlas {
namespace test {

//-----------------------------------------------------------------------------

CASE("test_pentagon_like_healpix") {
    Mesh mesh;

    // 1      2      3      4      5
    // +------+------+------+------+
    // |   0  |  1   |  2   |  3   |
    // +6     +7     +8     +9     +10
    //  \   /  \   /  \   /  \   /
    //    +      +      +      +
    //   11     12     13     14

    auto points = std::vector<PointLonLat>{
        {0., 90.},    // 1
        {90, 90.},    // 2
        {180., 90.},  // 3
        {270., 90.},  // 4
        {360., 90.},  // 5
        {0., 80.},    // 6
        {90, 80.},    // 7
        {180., 80.},  // 8
        {270., 80.},  // 9
        {360., 80.},  // 10
        {45., 70.},   // 11
        {135., 70.},  // 12
        {225., 70.},  // 13
        {315., 70.},  // 14
    };


    // nodes
    {
        mesh.nodes().resize(14);
        auto xy            = array::make_view<double, 2>(mesh.nodes().xy());
        auto lonlat        = array::make_view<double, 2>(mesh.nodes().lonlat());
        auto nodes_glb_idx = array::make_view<gidx_t, 1>(mesh.nodes().global_index());
        for (idx_t i = 0; i < mesh.nodes().size(); ++i) {
            lonlat(i, LON)   = points[i][LON];
            lonlat(i, LAT)   = points[i][LAT];
            nodes_glb_idx(i) = i + 1;
        }
        xy.assign(lonlat);
    }

    // elements
    {
        mesh.cells().add(mesh::ElementType::create("Pentagon"), 4);
        auto connect = [&](idx_t cell, std::array<idx_t, 5> nodes) {
            mesh.cells().node_connectivity().block(0).set(cell, nodes.data());
        };
        connect(0, {0, 5, 10, 6, 1});
        connect(1, {1, 6, 11, 7, 2});
        connect(2, {2, 7, 12, 8, 3});
        connect(3, {3, 8, 13, 9, 4});

        auto cell_glb_idx = array::make_view<gidx_t, 1>(mesh.cells().global_index());
        for (idx_t i = 0; i < mesh.cells().size(); ++i) {
            cell_glb_idx(i) = i + 1;
        }
    }

    {
        output::Gmsh gmsh("test_pentagon_xyz.msh", util::Config("coordinates", "xyz"));
        gmsh.write(mesh);
    }

    {
        output::Gmsh gmsh("test_pentagon_lonlat.msh", util::Config("coordinates", "lonlat"));
        gmsh.write(mesh);
    }

    auto src_fs = functionspace::NodeColumns(mesh);
    auto tgt_fs = functionspace::PointCloud({{45, 80}, {0, 89}, {180, 89}, {170, 89}});

    Interpolation interpolation(option::type("finite-element"), src_fs, tgt_fs);

    auto src_field = src_fs.createField<double>();
    auto tgt_field = tgt_fs.createField<double>();

    auto src_lonlat = array::make_view<double, 2>(src_fs.lonlat());
    auto src        = array::make_view<double, 1>(src_field);
    for (auto i = 0; i < src_fs.size(); ++i) {
        src(i) = std::cos(src_lonlat(i, LAT)) * std::sin(src_lonlat(i, LON));
    }

    interpolation.execute(src_field, tgt_field);

    auto tgt_lonlat = array::make_view<double, 2>(tgt_fs.lonlat());
    auto tgt        = array::make_view<double, 1>(tgt_field);
    for (auto i = 0; i < tgt_fs.size(); ++i) {
        double analytical = std::cos(tgt_lonlat(i, LAT)) * std::sin(tgt_lonlat(i, LON));
        Log::info() << i << " " << PointLonLat{tgt_lonlat(i, LON), tgt_lonlat(i, LAT)} << " : interpolated=" << tgt(i)
                    << "; expected=" << analytical << std::endl;
    }
}


//-----------------------------------------------------------------------------

}  // namespace test
}  // namespace atlas

int main(int argc, char** argv) {
    return atlas::test::run(argc, argv);
}
