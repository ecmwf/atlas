/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "atlas/grid/Grid.h"
#include "atlas/library/config.h"
#include "atlas/mesh/HybridElements.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/mesh/actions/BuildEdges.h"
#include "atlas/mesh/detail/AccumulateFacets.h"
#include "atlas/meshgenerator.h"
#include "atlas/option.h"
#include "atlas/util/Unique.h"

#include "tests/AtlasTestEnvironment.h"

using namespace atlas::mesh;
using namespace atlas::util;

namespace atlas {
namespace test {

//-----------------------------------------------------------------------------

CASE("test_accumulate_facets") {
    Grid grid("O2");
    StructuredMeshGenerator generator(Config("angle", 29.0)("triangulate", false)("ghost_at_end", false));

    Mesh mesh = generator.generate(grid);

    // storage for edge-to-node-connectivity shape=(nb_edges,2)
    std::vector<idx_t> edge_nodes_data;

    // storage for edge-to-cell-connectivity shape=(nb_edges,2)
    std::vector<idx_t> edge_to_cell_data;

    idx_t nb_edges;
    idx_t nb_inner_edges;
    idx_t missing_value;

    // Accumulate facets of cells ( edges in 2D )
    mesh::detail::accumulate_facets(mesh.cells(), mesh.nodes(), edge_nodes_data, edge_to_cell_data, nb_edges,
                                    nb_inner_edges, missing_value);

    idx_t edge_nodes_check[] = {
        0,  21, 21, 22, 22, 1,  1,  0,  22, 23, 23, 2,  2,  1,  3,  25, 25, 26, 26, 4,  4,  3,  26, 27, 27, 5,  5,
        4,  27, 28, 28, 6,  6,  5,  28, 29, 29, 7,  7,  6,  8,  31, 31, 32, 32, 9,  9,  8,  32, 33, 33, 10, 10, 9,
        33, 34, 34, 11, 11, 10, 34, 35, 35, 12, 12, 11, 13, 37, 37, 38, 38, 14, 14, 13, 38, 39, 39, 15, 15, 14, 39,
        40, 40, 16, 16, 15, 40, 41, 41, 17, 17, 16, 18, 43, 43, 44, 44, 19, 19, 18, 44, 45, 45, 20, 20, 19, 21, 46,
        46, 47, 47, 22, 47, 48, 48, 23, 48, 49, 49, 24, 24, 23, 49, 50, 50, 25, 25, 24, 50, 51, 51, 26, 51, 52, 52,
        27, 52, 53, 53, 28, 53, 54, 54, 29, 54, 55, 55, 30, 30, 29, 55, 56, 56, 31, 31, 30, 56, 57, 57, 32, 57, 58,
        58, 33, 58, 59, 59, 34, 59, 60, 60, 35, 60, 61, 61, 36, 36, 35, 61, 62, 62, 37, 37, 36, 62, 63, 63, 38, 63,
        64, 64, 39, 64, 65, 65, 40, 65, 66, 66, 41, 66, 67, 67, 42, 42, 41, 67, 68, 68, 43, 43, 42, 68, 69, 69, 44,
        69, 70, 70, 45, 46, 71, 71, 72, 72, 47, 72, 73, 73, 48, 50, 74, 74, 75, 75, 51, 75, 76, 76, 52, 76, 77, 77,
        53, 77, 78, 78, 54, 56, 79, 79, 80, 80, 57, 80, 81, 81, 58, 81, 82, 82, 59, 82, 83, 83, 60, 62, 84, 84, 85,
        85, 63, 85, 86, 86, 64, 86, 87, 87, 65, 87, 88, 88, 66, 68, 89, 89, 90, 90, 69, 90, 91, 91, 70, 24, 2,  24,
        3,  3,  2,  30, 7,  30, 8,  8,  7,  36, 12, 36, 13, 13, 12, 42, 17, 42, 18, 18, 17, 73, 49, 73, 74, 74, 49,
        78, 55, 78, 79, 79, 55, 83, 61, 83, 84, 84, 61, 88, 67, 88, 89, 89, 67};
    EXPECT(edge_nodes_data == eckit::testing::make_view(edge_nodes_check, edge_nodes_check + 2 * nb_edges));

    idx_t edge_to_cell_check[] = {0,  missing_value,
                                  0,  16,
                                  0,  1,
                                  0,  missing_value,
                                  1,  17,
                                  1,  56,
                                  1,  missing_value,
                                  2,  58,
                                  2,  20,
                                  2,  3,
                                  2,  missing_value,
                                  3,  21,
                                  3,  4,
                                  3,  missing_value,
                                  4,  22,
                                  4,  5,
                                  4,  missing_value,
                                  5,  23,
                                  5,  59,
                                  5,  missing_value,
                                  6,  61,
                                  6,  26,
                                  6,  7,
                                  6,  missing_value,
                                  7,  27,
                                  7,  8,
                                  7,  missing_value,
                                  8,  28,
                                  8,  9,
                                  8,  missing_value,
                                  9,  29,
                                  9,  62,
                                  9,  missing_value,
                                  10, 64,
                                  10, 32,
                                  10, 11,
                                  10, missing_value,
                                  11, 33,
                                  11, 12,
                                  11, missing_value,
                                  12, 34,
                                  12, 13,
                                  12, missing_value,
                                  13, 35,
                                  13, 65,
                                  13, missing_value,
                                  14, 67,
                                  14, 38,
                                  14, 15,
                                  14, missing_value,
                                  15, 39,
                                  15, missing_value,
                                  15, missing_value,
                                  16, missing_value,
                                  16, 40,
                                  16, 17,
                                  17, 41,
                                  17, 18,
                                  18, 68,
                                  18, 19,
                                  18, 56,
                                  19, 70,
                                  19, 20,
                                  19, 58,
                                  20, 42,
                                  20, 21,
                                  21, 43,
                                  21, 22,
                                  22, 44,
                                  22, 23,
                                  23, 45,
                                  23, 24,
                                  24, 71,
                                  24, 25,
                                  24, 59,
                                  25, 73,
                                  25, 26,
                                  25, 61,
                                  26, 46,
                                  26, 27,
                                  27, 47,
                                  27, 28,
                                  28, 48,
                                  28, 29,
                                  29, 49,
                                  29, 30,
                                  30, 74,
                                  30, 31,
                                  30, 62,
                                  31, 76,
                                  31, 32,
                                  31, 64,
                                  32, 50,
                                  32, 33,
                                  33, 51,
                                  33, 34,
                                  34, 52,
                                  34, 35,
                                  35, 53,
                                  35, 36,
                                  36, 77,
                                  36, 37,
                                  36, 65,
                                  37, 79,
                                  37, 38,
                                  37, 67,
                                  38, 54,
                                  38, 39,
                                  39, 55,
                                  39, missing_value,
                                  40, missing_value,
                                  40, missing_value,
                                  40, 41,
                                  41, missing_value,
                                  41, 68,
                                  42, 70,
                                  42, missing_value,
                                  42, 43,
                                  43, missing_value,
                                  43, 44,
                                  44, missing_value,
                                  44, 45,
                                  45, missing_value,
                                  45, 71,
                                  46, 73,
                                  46, missing_value,
                                  46, 47,
                                  47, missing_value,
                                  47, 48,
                                  48, missing_value,
                                  48, 49,
                                  49, missing_value,
                                  49, 74,
                                  50, 76,
                                  50, missing_value,
                                  50, 51,
                                  51, missing_value,
                                  51, 52,
                                  52, missing_value,
                                  52, 53,
                                  53, missing_value,
                                  53, 77,
                                  54, 79,
                                  54, missing_value,
                                  54, 55,
                                  55, missing_value,
                                  55, missing_value,
                                  56, 57,
                                  57, 58,
                                  57, missing_value,
                                  59, 60,
                                  60, 61,
                                  60, missing_value,
                                  62, 63,
                                  63, 64,
                                  63, missing_value,
                                  65, 66,
                                  66, 67,
                                  66, missing_value,
                                  68, 69,
                                  69, missing_value,
                                  69, 70,
                                  71, 72,
                                  72, missing_value,
                                  72, 73,
                                  74, 75,
                                  75, missing_value,
                                  75, 76,
                                  77, 78,
                                  78, missing_value,
                                  78, 79};
    EXPECT(edge_to_cell_data == eckit::testing::make_view(edge_to_cell_check, edge_to_cell_check + 2 * nb_edges));
}

CASE("test_build_edges") {
    idx_t missing_value = -1;
    Grid grid("O2");
    StructuredMeshGenerator generator(Config("angle", 29.0)("triangulate", false)("ghost_at_end", false));
    Mesh mesh = generator.generate(grid);

    // Accumulate facets of cells ( edges in 2D )
    mesh::actions::build_edges(mesh, option::pole_edges(false));

    std::vector<idx_t> edge_nodes_check{
        0,  21, 21, 22, 22, 1,  1,  0,  22, 23, 23, 2,  2,  1,  3,  25, 25, 26, 26, 4,  4,  3,  26, 27, 27, 5,  5,
        4,  27, 28, 28, 6,  6,  5,  28, 29, 29, 7,  7,  6,  8,  31, 31, 32, 32, 9,  9,  8,  32, 33, 33, 10, 10, 9,
        33, 34, 34, 11, 11, 10, 34, 35, 35, 12, 12, 11, 13, 37, 37, 38, 38, 14, 14, 13, 38, 39, 39, 15, 15, 14, 39,
        40, 40, 16, 16, 15, 40, 41, 41, 17, 17, 16, 18, 43, 43, 44, 44, 19, 19, 18, 44, 45, 45, 20, 20, 19, 21, 46,
        46, 47, 47, 22, 47, 48, 48, 23, 48, 49, 49, 24, 24, 23, 49, 50, 50, 25, 25, 24, 50, 51, 51, 26, 51, 52, 52,
        27, 52, 53, 53, 28, 53, 54, 54, 29, 54, 55, 55, 30, 30, 29, 55, 56, 56, 31, 31, 30, 56, 57, 57, 32, 57, 58,
        58, 33, 58, 59, 59, 34, 59, 60, 60, 35, 60, 61, 61, 36, 36, 35, 61, 62, 62, 37, 37, 36, 62, 63, 63, 38, 63,
        64, 64, 39, 64, 65, 65, 40, 65, 66, 66, 41, 66, 67, 67, 42, 42, 41, 67, 68, 68, 43, 43, 42, 68, 69, 69, 44,
        69, 70, 70, 45, 46, 71, 71, 72, 72, 47, 72, 73, 73, 48, 50, 74, 74, 75, 75, 51, 75, 76, 76, 52, 76, 77, 77,
        53, 77, 78, 78, 54, 56, 79, 79, 80, 80, 57, 80, 81, 81, 58, 81, 82, 82, 59, 82, 83, 83, 60, 62, 84, 84, 85,
        85, 63, 85, 86, 86, 64, 86, 87, 87, 65, 87, 88, 88, 66, 68, 89, 89, 90, 90, 69, 90, 91, 91, 70, 24, 2,  24,
        3,  3,  2,  30, 7,  30, 8,  8,  7,  36, 12, 36, 13, 13, 12, 42, 17, 42, 18, 18, 17, 73, 49, 73, 74, 74, 49,
        78, 55, 78, 79, 79, 55, 83, 61, 83, 84, 84, 61, 88, 67, 88, 89, 89, 67};

    {
        const mesh::HybridElements::Connectivity& edge_node_connectivity = mesh.edges().node_connectivity();
        EXPECT(mesh.projection().units() == "degrees");
        const util::UniqueLonLat compute_uid(mesh);
        EXPECT_EQ(mesh.edges().size(), edge_nodes_check.size() / 2);
        for (idx_t jedge = 0; jedge < mesh.edges().size(); ++jedge) {
            if (compute_uid(edge_nodes_check[2 * jedge + 0]) < compute_uid(edge_nodes_check[2 * jedge + 1])) {
                EXPECT_EQ(edge_nodes_check[2 * jedge + 0], edge_node_connectivity(jedge, 0));
                EXPECT_EQ(edge_nodes_check[2 * jedge + 1], edge_node_connectivity(jedge, 1));
            }
            else {
                EXPECT_EQ(edge_nodes_check[2 * jedge + 0], edge_node_connectivity(jedge, 1));
                EXPECT_EQ(edge_nodes_check[2 * jedge + 1], edge_node_connectivity(jedge, 0));
            }
        }
    }

    std::vector<idx_t> edge_to_cell_check{0,  missing_value,
                                          0,  16,
                                          0,  1,
                                          0,  missing_value,
                                          1,  17,
                                          1,  56,
                                          1,  missing_value,
                                          2,  58,
                                          2,  20,
                                          2,  3,
                                          2,  missing_value,
                                          3,  21,
                                          3,  4,
                                          3,  missing_value,
                                          4,  22,
                                          4,  5,
                                          4,  missing_value,
                                          5,  23,
                                          5,  59,
                                          5,  missing_value,
                                          6,  61,
                                          6,  26,
                                          6,  7,
                                          6,  missing_value,
                                          7,  27,
                                          7,  8,
                                          7,  missing_value,
                                          8,  28,
                                          8,  9,
                                          8,  missing_value,
                                          9,  29,
                                          9,  62,
                                          9,  missing_value,
                                          10, 64,
                                          10, 32,
                                          10, 11,
                                          10, missing_value,
                                          11, 33,
                                          11, 12,
                                          11, missing_value,
                                          12, 34,
                                          12, 13,
                                          12, missing_value,
                                          13, 35,
                                          13, 65,
                                          13, missing_value,
                                          14, 67,
                                          14, 38,
                                          14, 15,
                                          14, missing_value,
                                          15, 39,
                                          15, missing_value,
                                          15, missing_value,
                                          16, missing_value,
                                          16, 40,
                                          16, 17,
                                          17, 41,
                                          17, 18,
                                          18, 68,
                                          18, 19,
                                          18, 56,
                                          19, 70,
                                          19, 20,
                                          19, 58,
                                          20, 42,
                                          20, 21,
                                          21, 43,
                                          21, 22,
                                          22, 44,
                                          22, 23,
                                          23, 45,
                                          23, 24,
                                          24, 71,
                                          24, 25,
                                          24, 59,
                                          25, 73,
                                          25, 26,
                                          25, 61,
                                          26, 46,
                                          26, 27,
                                          27, 47,
                                          27, 28,
                                          28, 48,
                                          28, 29,
                                          29, 49,
                                          29, 30,
                                          30, 74,
                                          30, 31,
                                          30, 62,
                                          31, 76,
                                          31, 32,
                                          31, 64,
                                          32, 50,
                                          32, 33,
                                          33, 51,
                                          33, 34,
                                          34, 52,
                                          34, 35,
                                          35, 53,
                                          35, 36,
                                          36, 77,
                                          36, 37,
                                          36, 65,
                                          37, 79,
                                          37, 38,
                                          37, 67,
                                          38, 54,
                                          38, 39,
                                          39, 55,
                                          39, missing_value,
                                          40, missing_value,
                                          40, missing_value,
                                          40, 41,
                                          41, missing_value,
                                          41, 68,
                                          42, 70,
                                          42, missing_value,
                                          42, 43,
                                          43, missing_value,
                                          43, 44,
                                          44, missing_value,
                                          44, 45,
                                          45, missing_value,
                                          45, 71,
                                          46, 73,
                                          46, missing_value,
                                          46, 47,
                                          47, missing_value,
                                          47, 48,
                                          48, missing_value,
                                          48, 49,
                                          49, missing_value,
                                          49, 74,
                                          50, 76,
                                          50, missing_value,
                                          50, 51,
                                          51, missing_value,
                                          51, 52,
                                          52, missing_value,
                                          52, 53,
                                          53, missing_value,
                                          53, 77,
                                          54, 79,
                                          54, missing_value,
                                          54, 55,
                                          55, missing_value,
                                          55, missing_value,
                                          56, 57,
                                          57, 58,
                                          57, missing_value,
                                          59, 60,
                                          60, 61,
                                          60, missing_value,
                                          62, 63,
                                          63, 64,
                                          63, missing_value,
                                          65, 66,
                                          66, 67,
                                          66, missing_value,
                                          68, 69,
                                          69, missing_value,
                                          69, 70,
                                          71, 72,
                                          72, missing_value,
                                          72, 73,
                                          74, 75,
                                          75, missing_value,
                                          75, 76,
                                          77, 78,
                                          78, missing_value,
                                          78, 79};

    {
        const mesh::HybridElements::Connectivity& cell_node_connectivity = mesh.cells().node_connectivity();
        const mesh::HybridElements::Connectivity& edge_cell_connectivity = mesh.edges().cell_connectivity();
        const util::UniqueLonLat compute_uid(mesh);
        EXPECT_EQ(mesh.edges().size(), edge_to_cell_check.size() / 2);
        for (idx_t jedge = 0; jedge < mesh.edges().size(); ++jedge) {
            idx_t e1 = edge_to_cell_check[2 * jedge + 0];
            idx_t e2 = edge_to_cell_check[2 * jedge + 1];
            if (e2 == edge_cell_connectivity.missing_value() ||
                compute_uid(cell_node_connectivity.row(e1)) < compute_uid(cell_node_connectivity.row(e2))) {
                EXPECT_EQ(edge_to_cell_check[2 * jedge + 0], edge_cell_connectivity(jedge, 0));
                EXPECT_EQ(edge_to_cell_check[2 * jedge + 1], edge_cell_connectivity(jedge, 1));
            }
            else {
                std::cout << "jedge " << jedge << std::endl;
                EXPECT_EQ(edge_to_cell_check[2 * jedge + 0], edge_cell_connectivity(jedge, 1));
                EXPECT_EQ(edge_to_cell_check[2 * jedge + 1], edge_cell_connectivity(jedge, 0));
            }
        }
    }

    {
        const MultiBlockConnectivity& elem_edge_connectivity = mesh.cells().edge_connectivity();
        for (idx_t jelem = 0; jelem < mesh.cells().size(); ++jelem) {
            std::cout << std::setw(3) << jelem << " : ";
            for (idx_t jedge = 0; jedge < elem_edge_connectivity.cols(jelem); ++jedge) {
                std::cout << std::setw(3) << elem_edge_connectivity(jelem, jedge) << "  ";
            }
            std::cout << std::endl;
        }
    }
}

CASE("test_build_edges_triangles_only") {
    Grid grid("O2");
    StructuredMeshGenerator generator(Config("angle", 29.0)("triangulate", true)("ghost_at_end", false));
    Mesh mesh = generator.generate(grid);

    // Accumulate facets of cells ( edges in 2D )
    mesh::actions::build_edges(mesh, option::pole_edges(false));

    {
        const MultiBlockConnectivity& elem_edge_connectivity = mesh.cells().edge_connectivity();
        const MultiBlockConnectivity& elem_node_connectivity = mesh.cells().node_connectivity();
        for (idx_t jelem = 0; jelem < mesh.cells().size(); ++jelem) {
            std::cout << "elem" << std::setw(3) << jelem << " :   edges (  ";
            for (idx_t jedge = 0; jedge < elem_edge_connectivity.cols(jelem); ++jedge) {
                std::cout << std::setw(3) << elem_edge_connectivity(jelem, jedge) << "  ";
            }
            std::cout << ")      |      nodes ( ";
            for (idx_t jnode = 0; jnode < elem_node_connectivity.cols(jelem); ++jnode) {
                std::cout << std::setw(3) << elem_node_connectivity(jelem, jnode) << "  ";
            }
            std::cout << ")" << std::endl;
        }
    }
    std::cout << "( if you see all -1 entries, those are patch elements at the pole )" << std::endl;
}

//-----------------------------------------------------------------------------

CASE("test_pole_edge_default") {
    auto pole_edges = [](const Grid& grid) {
        auto mesh = StructuredMeshGenerator().generate(grid);
        mesh::actions::build_edges(mesh);
        return mesh.edges().metadata().getBool("pole_edges");
    };
    EXPECT(pole_edges(Grid("L10x11")) == false);
    EXPECT(pole_edges(Grid("F4")) == true);
    EXPECT(pole_edges(Grid("S4")) == true);
    EXPECT(pole_edges(Grid("Slat4")) == true);
    EXPECT(pole_edges(Grid("Slon4")) == false);
    EXPECT(pole_edges(Grid(Config("type", "regional")("nx", 35)("ny", 25)("north", -10)("south", -50)("east", 170)(
               "west", 100))) == false);
}

//-----------------------------------------------------------------------------

}  // namespace test
}  // namespace atlas

int main(int argc, char** argv) {
    return atlas::test::run(argc, argv);
}
