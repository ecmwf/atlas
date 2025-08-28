/*
 * (C) Copyright 2023 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <array>
#include <numeric>
#include <vector>

#include "atlas/array/MakeView.h"
#include "atlas/grid.h"
#include "atlas/mesh/Elements.h"
#include "atlas/mesh/HybridElements.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/mesh/MeshBuilder.h"
#include "atlas/mesh/Nodes.h"
#include "atlas/util/Config.h"

#include "tests/AtlasTestEnvironment.h"
#include "tests/mesh/helper_mesh_builder.h"

using namespace atlas::mesh;

//#include "atlas/output/Gmsh.h"
//using namespace atlas::output;

namespace atlas {
namespace test {

//-----------------------------------------------------------------------------

template <typename T>
mdspan<T,dims<1>> make_mdspan(std::vector<T>& v) {
    return mdspan<T,dims<1>>{v.data(), v.size()};
}
template <typename T, size_t N>
mdspan<T,extents<size_t,dynamic_extent,N>> make_mdspan(std::vector<std::array<T,N>>& v) {
    return mdspan<T,extents<size_t,dynamic_extent,N>>{reinterpret_cast<T*>(v.data()), v.size()};
}

//-----------------------------------------------------------------------------

CASE("test_tiny_mesh") {
    // small regional grid whose cell-centers are connected as (global nodes and cells):
    //
    //   1 - 5 ----- 6
    //   |  3  \ 1 /2|
    //   2 ----- 3 - 4
    //
    constexpr size_t nb_nodes = 6;
    std::vector<double> lons{{0.0, 0.0, 10.0, 15.0, 5.0, 15.0}};
    std::vector<double> lats{{5.0, 0.0, 0.0, 0.0, 5.0, 5.0}};

    std::vector<int> ghosts(nb_nodes, 0);  // all points owned
    std::vector<gidx_t> global_indices(nb_nodes);
    std::iota(global_indices.begin(), global_indices.end(), 1);  // 1-based numbering
    const idx_t remote_index_base = 0;                           // 0-based numbering
    std::vector<idx_t> remote_indices(nb_nodes);
    std::iota(remote_indices.begin(), remote_indices.end(), remote_index_base);
    std::vector<int> partitions(nb_nodes, 0);  // all points on proc 0

    // triangles
    std::vector<std::array<gidx_t, 3>> tri_boundary_nodes = {{{3, 6, 5}}, {{3, 4, 6}}};
    std::vector<gidx_t> tri_global_indices                = {1, 2};

    // quads
    std::vector<std::array<gidx_t, 4>> quad_boundary_nodes = {{{1, 2, 3, 5}}};
    std::vector<gidx_t> quad_global_indices                = {3};

    gidx_t global_index_base = 1;
    const MeshBuilder mesh_builder{};

    // const Mesh mesh = mesh_builder(lons, lats, ghosts, global_indices, remote_indices, remote_index_base, partitions,
    //                                tri_boundary_nodes, tri_global_indices, quad_boundary_nodes, quad_global_indices);
    const Mesh mesh = mesh_builder(
        make_mdspan(global_indices),
        make_mdspan(lons), make_mdspan(lats),
        make_mdspan(lons), make_mdspan(lats),
        make_mdspan(ghosts), make_mdspan(partitions), make_mdspan(remote_indices), remote_index_base,
        make_mdspan(tri_global_indices), make_mdspan(tri_boundary_nodes),
        make_mdspan(quad_global_indices), make_mdspan(quad_boundary_nodes),
        global_index_base);

    helper::check_mesh_nodes_and_cells(mesh, lons, lats, ghosts, global_indices, remote_indices, remote_index_base,
                                       partitions, tri_boundary_nodes, tri_global_indices, quad_boundary_nodes,
                                       quad_global_indices);

    //Gmsh gmsh("out.msh", util::Config("coordinates", "xyz"));
    //gmsh.write(mesh);
}

CASE("test_cs_c2_mesh_serial") {
    // coordinates of C2 lfric cubed-sphere grid: grid("CS-LFR-2");
    std::vector<double> lons = {337.5, 22.5,  337.5, 22.5,   // +x
                                67.5,  112.5, 67.5,  112.5,  // +y
                                202.5, 202.5, 157.5, 157.5,  // -x
                                292.5, 292.5, 247.5, 247.5,  // -y
                                315,   45,    225,   135,    // +z
                                315,   225,   45,    135};   // -z

    std::vector<double> lats = {-20.941,  -20.941,  20.941,   20.941,     // +x
                                -20.941,  -20.941,  20.941,   20.941,     // +y
                                -20.941,  20.941,   -20.941,  20.941,     // -x
                                -20.941,  20.941,   -20.941,  20.941,     // -y
                                59.6388,  59.6388,  59.6388,  59.6388,    // +z
                                -59.6388, -59.6388, -59.6388, -59.6388};  // -z

    gidx_t global_index_base = 1;
    std::vector<int> ghosts(24, 0);
    std::vector<gidx_t> global_indices(24);
    std::iota(global_indices.begin(), global_indices.end(), global_index_base);
    const idx_t remote_index_base = 1;  // test with 1-based numbering
    std::vector<idx_t> remote_indices(24);
    std::iota(remote_indices.begin(), remote_indices.end(), remote_index_base);
    std::vector<int> partitions(24, 0);

    // triangles
    std::vector<std::array<gidx_t, 3>> tri_boundary_nodes = {//corners
                                                             {{17, 14, 3}},
                                                             {{18, 4, 7}},
                                                             {{20, 8, 12}},
                                                             {{19, 10, 16}},
                                                             {{21, 1, 13}},
                                                             {{23, 5, 2}},
                                                             {{24, 11, 6}},
                                                             {{22, 15, 9}}};
    std::vector<gidx_t> tri_global_indices(8);
    std::iota(tri_global_indices.begin(), tri_global_indices.end(), global_index_base);

    // quads
    std::vector<std::array<gidx_t, 4>> quad_boundary_nodes = {// faces
                                                              {{1, 2, 4, 3}},
                                                              {{5, 6, 8, 7}},
                                                              {{11, 9, 10, 12}},
                                                              {{15, 13, 14, 16}},
                                                              {{17, 18, 20, 19}},
                                                              {{21, 22, 24, 23}},
                                                              // edges between faces
                                                              {{2, 5, 7, 4}},
                                                              {{6, 11, 12, 8}},
                                                              {{9, 15, 16, 10}},
                                                              {{13, 1, 3, 14}},
                                                              {{7, 8, 20, 18}},
                                                              {{12, 10, 19, 20}},
                                                              {{16, 14, 17, 19}},
                                                              {{3, 4, 18, 17}},
                                                              {{23, 24, 6, 5}},
                                                              {{24, 22, 9, 11}},
                                                              {{22, 21, 13, 15}},
                                                              {{21, 23, 2, 1}}};
    std::vector<gidx_t> quad_global_indices(18);
    std::iota(quad_global_indices.begin(), quad_global_indices.end(), global_index_base + tri_global_indices.size());

    const MeshBuilder mesh_builder{};

    SECTION("Build Mesh without a Grid") {
        const Mesh mesh = mesh_builder(
            make_mdspan(global_indices),
            make_mdspan(lons), make_mdspan(lats),
            make_mdspan(lons), make_mdspan(lats),
            make_mdspan(ghosts), make_mdspan(partitions), make_mdspan(remote_indices), remote_index_base,
            make_mdspan(tri_global_indices), make_mdspan(tri_boundary_nodes),
            make_mdspan(quad_global_indices), make_mdspan(quad_boundary_nodes),
            global_index_base);

        helper::check_mesh_nodes_and_cells(mesh, lons, lats, ghosts, global_indices, remote_indices, remote_index_base,
                                           partitions, tri_boundary_nodes, tri_global_indices, quad_boundary_nodes,
                                           quad_global_indices);
        EXPECT(!mesh.grid());

        //Gmsh gmsh("out.msh", util::Config("coordinates", "xyz"));
        //gmsh.write(mesh);
    }

    SECTION("Build Mesh with an UnstructuredGrid from config") {
        std::vector<double> lonlats(2 * lons.size());
        int counter = 0;
        for (size_t i = 0; i < lons.size(); ++i) {
            lonlats[counter] = lons[i];
            counter++;
            lonlats[counter] = lats[i];
            counter++;
        }

        util::Config config{};
        config.set("grid.type", "unstructured");
        config.set("grid.xy", lonlats);
        config.set("validate", true);

        const Mesh mesh = mesh_builder(
            make_mdspan(global_indices),
            make_mdspan(lons), make_mdspan(lats),
            make_mdspan(lons), make_mdspan(lats),
            make_mdspan(ghosts), make_mdspan(partitions), make_mdspan(remote_indices), remote_index_base,
            make_mdspan(tri_global_indices), make_mdspan(tri_boundary_nodes),
            make_mdspan(quad_global_indices), make_mdspan(quad_boundary_nodes),
            global_index_base,
            config);

        helper::check_mesh_nodes_and_cells(mesh, lons, lats, ghosts, global_indices, remote_indices, remote_index_base,
                                           partitions, tri_boundary_nodes, tri_global_indices, quad_boundary_nodes,
                                           quad_global_indices);

        EXPECT(mesh.grid());
        EXPECT(mesh.grid().type() == "unstructured");
    }

    SECTION("Build Mesh with an UnstructuredGrid, with Grid assembled from MeshBuilder arguments") {
        util::Config config{};
        config.set("grid.type", "unstructured");

        const Mesh mesh = mesh_builder(
            make_mdspan(global_indices),
            make_mdspan(lons), make_mdspan(lats),
            make_mdspan(lons), make_mdspan(lats),
            make_mdspan(ghosts), make_mdspan(partitions), make_mdspan(remote_indices), remote_index_base,
            make_mdspan(tri_global_indices), make_mdspan(tri_boundary_nodes),
            make_mdspan(quad_global_indices), make_mdspan(quad_boundary_nodes),
            global_index_base,
            config);



        helper::check_mesh_nodes_and_cells(mesh, lons, lats, ghosts, global_indices, remote_indices, remote_index_base,
                                           partitions, tri_boundary_nodes, tri_global_indices, quad_boundary_nodes,
                                           quad_global_indices);

        EXPECT(mesh.grid());
        EXPECT(mesh.grid().type() == "unstructured");
    }

    SECTION("Build Mesh with a CubedSphereGrid from config") {  // but can't validate this
        util::Config config{};
        config.set("grid.name", "CS-LFR-2");

        const Mesh mesh = mesh_builder(
            make_mdspan(global_indices),
            make_mdspan(lons), make_mdspan(lats),
            make_mdspan(lons), make_mdspan(lats),
            make_mdspan(ghosts), make_mdspan(partitions), make_mdspan(remote_indices), remote_index_base,
            make_mdspan(tri_global_indices), make_mdspan(tri_boundary_nodes),
            make_mdspan(quad_global_indices), make_mdspan(quad_boundary_nodes),
            global_index_base,
            config);

        helper::check_mesh_nodes_and_cells(mesh, lons, lats, ghosts, global_indices, remote_indices, remote_index_base,
                                           partitions, tri_boundary_nodes, tri_global_indices, quad_boundary_nodes,
                                           quad_global_indices);

        EXPECT(mesh.grid());
        EXPECT(mesh.grid().type() == "cubedsphere");
    }
}

//-----------------------------------------------------------------------------

}  // namespace test
}  // namespace atlas

int main(int argc, char** argv) {
    return atlas::test::run(argc, argv);
}
