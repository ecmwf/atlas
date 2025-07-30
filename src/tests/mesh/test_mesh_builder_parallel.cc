/*
 * (C) Copyright 2023 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <algorithm>
#include <array>
#include <iterator>
#include <vector>

#include "atlas/array/MakeView.h"
#include "atlas/mesh/Elements.h"
#include "atlas/mesh/HybridElements.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/mesh/MeshBuilder.h"
#include "atlas/mesh/Nodes.h"
#include "atlas/parallel/mpi/mpi.h"

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

CASE("test_cs_c2_mesh_parallel") {
    ATLAS_ASSERT(mpi::comm().size() == 6);
    const int rank = mpi::comm().rank();

    constexpr gidx_t global_index_base = 1;

    // Coordinates of the C2 LFRic cubed-sphere grid: grid("CS-LFR-2");
    const std::vector<double> global_lons = {337.5, 22.5,  337.5, 22.5,   // +x
                                             67.5,  112.5, 67.5,  112.5,  // +y
                                             202.5, 202.5, 157.5, 157.5,  // -x
                                             292.5, 292.5, 247.5, 247.5,  // -y
                                             315,   45,    225,   135,    // +z
                                             315,   225,   45,    135};   // -z

    const std::vector<double> global_lats = {-20.941,  -20.941,  20.941,   20.941,     // +x
                                             -20.941,  -20.941,  20.941,   20.941,     // +y
                                             -20.941,  20.941,   -20.941,  20.941,     // -x
                                             -20.941,  20.941,   -20.941,  20.941,     // -y
                                             59.6388,  59.6388,  59.6388,  59.6388,    // +z
                                             -59.6388, -59.6388, -59.6388, -59.6388};  // -z

    const std::vector<int> global_partitions = {0, 0, 0, 0,   // +x
                                                1, 1, 1, 1,   // +y
                                                2, 2, 2, 2,   // -x
                                                3, 3, 3, 3,   // -y
                                                4, 4, 4, 4,   // +z
                                                5, 5, 5, 5};  // -z

    const std::vector<int> local_indices = {0, 1, 2, 3,   // +x
                                            0, 1, 2, 3,   // +y
                                            0, 1, 2, 3,   // -x
                                            0, 1, 2, 3,   // -y
                                            0, 1, 2, 3,   // +z
                                            0, 1, 2, 3};  // -z

    std::vector<gidx_t> global_indices;
    std::vector<std::array<gidx_t, 3>> tri_boundary_nodes{};
    std::vector<gidx_t> tri_global_indices{};
    std::vector<std::array<gidx_t, 4>> quad_boundary_nodes{};
    std::vector<gidx_t> quad_global_indices{};

    if (rank == 0) {
        // global indices for points and cells are 1-based
        global_indices = {1, 2, 3, 4, 5, 7, 17, 18};
        tri_boundary_nodes.push_back({{18, 4, 7}});
        tri_global_indices = {2};
        quad_boundary_nodes.push_back({{1, 2, 4, 3}});
        quad_boundary_nodes.push_back({{2, 5, 7, 4}});
        quad_boundary_nodes.push_back({{3, 4, 18, 17}});
        quad_global_indices = {9, 15, 22};
    }
    else if (rank == 1) {
        global_indices = {5, 6, 7, 8, 11, 12, 18, 20};
        tri_boundary_nodes.push_back({{20, 8, 12}});
        tri_global_indices = {3};
        quad_boundary_nodes.push_back({{5, 6, 8, 7}});
        quad_boundary_nodes.push_back({{6, 11, 12, 8}});
        quad_boundary_nodes.push_back({{7, 8, 20, 18}});
        quad_global_indices = {10, 16, 19};
    }
    else if (rank == 2) {
        global_indices = {9, 10, 11, 12, 15, 16, 22, 24};
        tri_boundary_nodes.push_back({{22, 15, 9}});
        tri_global_indices = {8};
        quad_boundary_nodes.push_back({{11, 9, 10, 12}});
        quad_boundary_nodes.push_back({{9, 15, 16, 10}});
        quad_boundary_nodes.push_back({{24, 22, 9, 11}});
        quad_global_indices = {11, 17, 24};
    }
    else if (rank == 3) {
        global_indices = {1, 3, 13, 14, 15, 16, 21, 22};
        tri_boundary_nodes.push_back({{21, 1, 13}});
        tri_global_indices = {5};
        quad_boundary_nodes.push_back({{15, 13, 14, 16}});
        quad_boundary_nodes.push_back({{13, 1, 3, 14}});
        quad_boundary_nodes.push_back({{22, 21, 13, 15}});
        quad_global_indices = {12, 18, 25};
    }
    else if (rank == 4) {
        global_indices = {3, 10, 12, 14, 16, 17, 18, 19, 20};
        tri_boundary_nodes.push_back({{17, 14, 3}});
        tri_boundary_nodes.push_back({{19, 10, 16}});
        tri_global_indices = {1, 4};
        quad_boundary_nodes.push_back({{17, 18, 20, 19}});
        quad_boundary_nodes.push_back({{12, 10, 19, 20}});
        quad_boundary_nodes.push_back({{16, 14, 17, 19}});
        quad_global_indices = {13, 20, 21};
    }
    else {  // rank == 5
        global_indices = {1, 2, 5, 6, 11, 21, 22, 23, 24};
        tri_boundary_nodes.push_back({{23, 5, 2}});
        tri_boundary_nodes.push_back({{24, 11, 6}});
        tri_global_indices = {6, 7};
        quad_boundary_nodes.push_back({{21, 22, 24, 23}});
        quad_boundary_nodes.push_back({{23, 24, 6, 5}});
        quad_boundary_nodes.push_back({{21, 23, 2, 1}});
        quad_global_indices = {14, 23, 26};
    }

    // Compute (local subset of) {lons,lats,ghosts,partitions} from (local subset of) global indices
    std::vector<double> lons;
    std::transform(global_indices.begin(), global_indices.end(), std::back_inserter(lons),
                   [&global_lons](const gidx_t idx) { return global_lons[idx - 1]; });

    std::vector<double> lats;
    std::transform(global_indices.begin(), global_indices.end(), std::back_inserter(lats),
                   [&global_lats](const gidx_t idx) { return global_lats[idx - 1]; });

    std::vector<int> ghosts;
    std::transform(
        global_indices.begin(), global_indices.end(), std::back_inserter(ghosts),
        [&global_partitions, &rank](const gidx_t idx) { return static_cast<int>(global_partitions[idx - 1] != rank); });

    const idx_t remote_index_base = 0;  // 0-based indexing used in local_indices above
    std::vector<idx_t> remote_indices;
    std::transform(global_indices.begin(), global_indices.end(), std::back_inserter(remote_indices),
                   [&local_indices](const gidx_t idx) { return static_cast<idx_t>(local_indices[idx - 1]); });

    std::vector<int> partitions;
    std::transform(global_indices.begin(), global_indices.end(), std::back_inserter(partitions),
                   [&global_partitions](const gidx_t idx) { return global_partitions[idx - 1]; });

    const MeshBuilder mesh_builder{};

    SECTION("Build Mesh without a Grid") {
        const Mesh mesh = mesh_builder(
            make_mdspan(global_indices),
            make_mdspan(lons), make_mdspan(lats), make_mdspan(lons), make_mdspan(lats),
            make_mdspan(ghosts), make_mdspan(partitions), make_mdspan(remote_indices), remote_index_base,
            make_mdspan(tri_global_indices), make_mdspan(tri_boundary_nodes),
            make_mdspan(quad_global_indices), make_mdspan(quad_boundary_nodes),
            global_index_base);
        // const Mesh mesh = mesh_builder(lons, lats, ghosts, global_indices, remote_indices, remote_index_base, partitions,
        //                                tri_boundary_nodes, tri_global_indices, quad_boundary_nodes, quad_global_indices);

        helper::check_mesh_nodes_and_cells(mesh, lons, lats, ghosts, global_indices, remote_indices, remote_index_base,
                                           partitions, tri_boundary_nodes, tri_global_indices, quad_boundary_nodes,
                                           quad_global_indices);
        EXPECT(!mesh.grid());

        //Gmsh gmsh("out.msh", util::Config("coordinates", "xyz"));
        //gmsh.write(mesh);
    }

    SECTION("Build Mesh with an UnstructuredGrid from config") {
        const size_t nb_owned_nodes = std::count(ghosts.begin(), ghosts.end(), 0);
        std::vector<double> owned_lonlats(2 * nb_owned_nodes);
        int counter = 0;
        for (size_t i = 0; i < lons.size(); ++i) {
            if (ghosts[i] == 0) {
                owned_lonlats[counter] = lons[i];
                counter++;
                owned_lonlats[counter] = lats[i];
                counter++;
            }
        }
        size_t nb_nodes_global = 0;
        mpi::comm().allReduce(nb_owned_nodes, nb_nodes_global, eckit::mpi::sum());
        std::vector<double> global_lonlats(2 * nb_nodes_global);
        eckit::mpi::Buffer<double> buffer(mpi::comm().size());
        mpi::comm().allGatherv(owned_lonlats.begin(), owned_lonlats.end(), buffer);
        global_lonlats = std::move(buffer.buffer);

        util::Config config{};
        config.set("grid.type", "unstructured");
        config.set("grid.xy", global_lonlats);
        config.set("validate", true);

        // const Mesh mesh = mesh_builder(lons, lats, ghosts, global_indices, remote_indices, remote_index_base, partitions,
        //                                tri_boundary_nodes, tri_global_indices, quad_boundary_nodes, quad_global_indices,
        //                                config);
        const Mesh mesh = mesh_builder(
            make_mdspan(global_indices),
            make_mdspan(lons), make_mdspan(lats), make_mdspan(lons), make_mdspan(lats),
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

        // const Mesh mesh = mesh_builder(lons, lats, ghosts, global_indices, remote_indices, remote_index_base, partitions,
        //                                tri_boundary_nodes, tri_global_indices, quad_boundary_nodes, quad_global_indices,
        //                                config);
        const Mesh mesh = mesh_builder(
            make_mdspan(global_indices),
            make_mdspan(lons), make_mdspan(lats), make_mdspan(lons), make_mdspan(lats),
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
}

//-----------------------------------------------------------------------------

}  // namespace test
}  // namespace atlas

int main(int argc, char** argv) {
    return atlas::test::run(argc, argv);
}
