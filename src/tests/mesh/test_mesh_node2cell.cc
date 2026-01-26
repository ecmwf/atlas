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
#include "atlas/mesh/actions/BuildHalo.h"
#include "atlas/mesh/actions/BuildNode2CellConnectivity.h"
#include "atlas/mesh/actions/BuildParallelFields.h"
#include "atlas/mesh/actions/BuildPeriodicBoundaries.h"
#include "atlas/mesh/detail/AccumulateFacets.h"
#include "atlas/meshgenerator.h"
#include "atlas/option.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/util/Unique.h"

#include "tests/AtlasTestEnvironment.h"

using namespace atlas::mesh;
using namespace atlas::util;

namespace atlas {
namespace test {

//-----------------------------------------------------------------------------

/** print_check()
 * This routine is used to print code that can be copy-pasted to check results
 */
#if 0
void print_check( Mesh& mesh ) {
    std::ostream& out = std::cout;
    auto cell_glb_idx = array::make_view<gidx_t, 1>( mesh.cells().global_index() );

    idx_t size{0};
    for ( idx_t jnode = 0; jnode < mesh.nodes().size(); ++jnode ) {
        for ( idx_t jcell = 0; jcell < mesh.nodes().cell_connectivity().cols( jnode ); ++jcell ) {
            ++size;
        }
    }

    for ( size_t p = 0; p < mpi::comm().size(); ++p ) {
        mpi::comm().barrier();
        if ( mpi::comm().rank() == p ) {
            out << newl << "if( mpi::comm().rank() == " << mpi::comm().rank()
                << " && mpi::comm().size() == " << mpi::comm().size() << ") {" << newl;
            out << "auto node_cell_connectivity_check = std::array<gidx_t," << size << "> { ";
            for ( idx_t jnode = 0; jnode < mesh.nodes().size(); ++jnode ) {
                for ( idx_t jcell = 0; jcell < mesh.nodes().cell_connectivity().cols( jnode ); ++jcell ) {
                    idx_t icell = mesh.nodes().cell_connectivity()( jnode, jcell );
                    out << cell_glb_idx( icell ) << ", ";
                }
            }
            out << " };" << newl;
            out << "check_node_cell_connectivity( mesh, node_cell_connectivity_check );" << newl;
            out << "}" << std::endl;
        }
    }
}
#endif

//-----------------------------------------------------------------------------

// Routine that checks if the connectivity-table matches a given array "check"
template <typename T>
void check_node_cell_connectivity(Mesh& mesh, const T& check) {
    auto cell_glb_idx = array::make_view<gidx_t, 1>(mesh.cells().global_index());

    size_t c{0};
    for (idx_t jnode = 0; jnode < mesh.nodes().size(); ++jnode) {
        for (idx_t jcell = 0; jcell < mesh.nodes().cell_connectivity().cols(jnode); ++jcell) {
            idx_t icell = mesh.nodes().cell_connectivity()(jnode, jcell);
            EXPECT(c < check.size());
            EXPECT(cell_glb_idx(icell) == check[c++]);
        }
    }
    EXPECT(c == check.size());
}

//-----------------------------------------------------------------------------

struct PrettyPrintNodeCellConnectivity {
    const Mesh& mesh;
    PrettyPrintNodeCellConnectivity(const Mesh& _mesh): mesh{_mesh} {}

    void print(std::ostream& out) const {
        auto node_glb_idx = array::make_view<gidx_t, 1>(mesh.nodes().global_index());
        auto cell_glb_idx = array::make_view<gidx_t, 1>(mesh.cells().global_index());

        for (idx_t jnode = 0; jnode < mesh.nodes().size(); ++jnode) {
            out << node_glb_idx(jnode) << " : ";
            for (idx_t jcell = 0; jcell < mesh.nodes().cell_connectivity().cols(jnode); ++jcell) {
                idx_t icell = mesh.nodes().cell_connectivity()(jnode, jcell);
                out << cell_glb_idx(icell) << " ";
            }
            out << '\n';
        }
    }

    friend std::ostream& operator<<(std::ostream& out, const PrettyPrintNodeCellConnectivity& p) {
        p.print(out);
        return out;
    }
};


//-----------------------------------------------------------------------------

CASE("test_node2cell") {
    // Create mesh
    Grid grid("O2");
    StructuredMeshGenerator generator;
    Mesh mesh = generator.generate(grid);

    // Create node to cell connectivity within mesh
    mesh::actions::build_node_to_cell_connectivity(mesh);

    Log::info() << "Connectivity without halo :\n" << PrettyPrintNodeCellConnectivity{mesh} << std::endl;

    // Following command generates (prints) code between braces following
    //print_check( mesh );

    // mpi-size 1
    {
        if (mpi::comm().rank() == 0 && mpi::comm().size() == 1) {
            auto node_cell_connectivity_check = std::array<gidx_t, 468>{
                113, 26,  25,  113, 115, 114, 26,  28,  27,  115, 117, 116, 28,  30,  29,  117, 119, 118, 30,  32,  31,
                119, 121, 120, 32,  34,  33,  121, 123, 122, 34,  37,  35,  36,  123, 125, 124, 37,  39,  38,  125, 127,
                126, 39,  41,  40,  127, 129, 128, 41,  43,  42,  129, 130, 43,  45,  44,  130, 45,  48,  46,  47,  129,
                130, 128, 48,  50,  49,  127, 128, 126, 50,  52,  51,  125, 126, 124, 52,  54,  53,  123, 124, 122, 54,
                56,  55,  121, 122, 120, 56,  59,  57,  58,  119, 120, 118, 59,  61,  60,  117, 118, 116, 61,  63,  62,
                115, 116, 114, 63,  65,  64,  113, 114, 65,  67,  66,  25,  1,   26,  25,  27,  1,   2,   28,  27,  29,
                2,   3,   30,  29,  31,  3,   4,   32,  31,  33,  4,   5,   34,  33,  35,  5,   6,   35,  36,  6,   7,
                37,  36,  38,  7,   8,   39,  38,  40,  8,   9,   41,  40,  42,  9,   10,  43,  42,  44,  10,  11,  45,
                44,  46,  11,  12,  46,  47,  12,  13,  48,  47,  49,  13,  14,  50,  49,  51,  14,  15,  52,  51,  53,
                15,  16,  54,  53,  55,  16,  17,  56,  55,  57,  17,  18,  57,  58,  18,  19,  59,  58,  60,  19,  20,
                61,  60,  62,  20,  21,  63,  62,  64,  21,  22,  65,  64,  66,  22,  23,  67,  66,  68,  23,  24,  1,
                69,  1,   2,   69,  71,  70,  2,   3,   71,  73,  72,  3,   4,   73,  75,  74,  4,   5,   75,  77,  76,
                5,   6,   77,  79,  78,  6,   7,   79,  80,  7,   8,   80,  82,  81,  8,   9,   82,  84,  83,  9,   10,
                84,  86,  85,  10,  11,  86,  88,  87,  11,  12,  88,  90,  89,  12,  13,  90,  91,  13,  14,  91,  93,
                92,  14,  15,  93,  95,  94,  15,  16,  95,  97,  96,  16,  17,  97,  99,  98,  17,  18,  99,  101, 100,
                18,  19,  101, 102, 19,  20,  102, 104, 103, 20,  21,  104, 106, 105, 21,  22,  106, 108, 107, 22,  23,
                108, 110, 109, 23,  24,  110, 112, 111, 69,  70,  131, 71,  70,  72,  131, 133, 132, 73,  72,  74,  133,
                135, 134, 75,  74,  76,  135, 137, 136, 77,  76,  78,  137, 139, 138, 79,  80,  78,  81,  139, 141, 140,
                82,  81,  83,  141, 143, 142, 84,  83,  85,  143, 145, 144, 86,  85,  87,  145, 147, 146, 88,  87,  89,
                147, 148, 90,  91,  89,  92,  148, 93,  92,  94,  147, 148, 146, 95,  94,  96,  145, 146, 144, 97,  96,
                98,  143, 144, 142, 99,  98,  100, 141, 142, 140, 101, 102, 100, 103, 139, 140, 138, 104, 103, 105, 137,
                138, 136, 106, 105, 107, 135, 136, 134, 108, 107, 109, 133, 134, 132, 110, 109, 111, 131, 132, 67,  68,
                68,  24,  24,  112, 112, 111,
            };
            check_node_cell_connectivity(mesh, node_cell_connectivity_check);
        }
    }
    // mpi-size 4
    {
        if (mpi::comm().rank() == 0 && mpi::comm().size() == 4) {
            auto node_cell_connectivity_check = std::array<gidx_t, 187>{
                45, 3,  2,  45, 47, 46, 3,  5,  4,  47, 49, 48, 5,  7,  6,  49, 51, 50, 7,  9,  8,  51, 53, 52,
                9,  11, 10, 53, 55, 54, 11, 14, 12, 13, 55, 57, 56, 14, 16, 15, 57, 59, 58, 16, 18, 17, 59, 61,
                60, 18, 20, 19, 61, 62, 20, 22, 21, 62, 22, 25, 23, 24, 61, 62, 60, 25, 27, 26, 59, 60, 58, 27,
                29, 28, 57, 58, 56, 29, 31, 30, 55, 56, 54, 31, 33, 32, 53, 54, 52, 33, 36, 34, 35, 51, 52, 50,
                36, 38, 37, 49, 50, 48, 38, 40, 39, 47, 48, 46, 40, 42, 41, 45, 46, 42, 44, 43, 2,  1,  3,  2,
                4,  1,  44, 5,  4,  6,  7,  6,  8,  9,  8,  10, 11, 10, 12, 12, 13, 14, 13, 15, 16, 15, 17, 18,
                17, 19, 20, 19, 21, 22, 21, 23, 23, 24, 25, 24, 26, 27, 26, 28, 29, 28, 30, 31, 30, 32, 33, 32,
                34, 34, 35, 36, 35, 37, 38, 37, 39, 40, 39, 41, 42, 41, 43, 44, 43, 1,  1,
            };
            check_node_cell_connectivity(mesh, node_cell_connectivity_check);
        }

        if (mpi::comm().rank() == 1 && mpi::comm().size() == 4) {
            auto node_cell_connectivity_check = std::array<gidx_t, 44>{
                63, 64, 64, 65, 65, 66, 66, 67, 67, 68, 68, 69, 69, 70, 70, 71, 71, 72, 72, 73, 63, 63,
                64, 64, 65, 65, 66, 66, 67, 67, 68, 68, 69, 69, 70, 70, 71, 71, 72, 72, 73, 63, 73, 73,
            };
            check_node_cell_connectivity(mesh, node_cell_connectivity_check);
        }

        if (mpi::comm().rank() == 2 && mpi::comm().size() == 4) {
            auto node_cell_connectivity_check = std::array<gidx_t, 51>{
                74, 74, 75, 75, 76, 76, 77, 77, 78, 78, 79, 79, 80, 80, 81, 81, 82, 82, 83, 83, 84, 86, 84, 85, 74, 74,
                75, 75, 76, 76, 77, 77, 78, 78, 79, 79, 80, 80, 81, 81, 82, 82, 83, 86, 86, 85, 83, 84, 84, 85, 85,
            };
            check_node_cell_connectivity(mesh, node_cell_connectivity_check);
        }

        if (mpi::comm().rank() == 3 && mpi::comm().size() == 4) {
            auto node_cell_connectivity_check = std::array<gidx_t, 186>{
                126, 128, 127, 128, 130, 129, 87,  88,  131, 89,  88,  90,  131, 133, 132, 91,  90,  92,  133, 135, 134,
                93,  92,  94,  135, 137, 136, 95,  94,  96,  137, 139, 138, 97,  98,  96,  99,  139, 141, 140, 100, 99,
                101, 141, 143, 142, 102, 101, 103, 143, 145, 144, 104, 103, 105, 145, 147, 146, 106, 105, 107, 147, 148,
                108, 109, 107, 110, 148, 111, 110, 112, 147, 148, 146, 113, 112, 114, 145, 146, 144, 115, 114, 116, 143,
                144, 142, 117, 116, 118, 141, 142, 140, 119, 120, 118, 121, 139, 140, 138, 122, 121, 123, 137, 138, 136,
                124, 123, 125, 135, 136, 134, 126, 125, 127, 133, 134, 132, 128, 127, 129, 131, 132, 87,  87,  89,  88,
                89,  91,  90,  91,  93,  92,  93,  95,  94,  95,  97,  96,  97,  98,  98,  100, 99,  100, 102, 101, 102,
                104, 103, 104, 106, 105, 106, 108, 107, 108, 109, 109, 111, 110, 111, 113, 112, 113, 115, 114, 115, 117,
                116, 117, 119, 118, 119, 120, 120, 122, 121, 122, 124, 123, 124, 126, 125, 130, 130, 129,
            };
            check_node_cell_connectivity(mesh, node_cell_connectivity_check);
        }
    }
}

//-----------------------------------------------------------------------------

CASE("test_node2cell_with_halo") {
    Grid grid("O2");
    StructuredMeshGenerator generator;
    Mesh mesh = generator.generate(grid);

    mesh::actions::build_parallel_fields(mesh);
    mesh::actions::build_periodic_boundaries(mesh);
    mesh::actions::build_halo(mesh, 1);

    mesh::actions::build_node_to_cell_connectivity(mesh);

    Log::info() << "Connectivity with halo :\n" << PrettyPrintNodeCellConnectivity{mesh} << std::endl;


    // Following command generates (prints) code between braces following
    // print_check( mesh );

    // mpi-size 1
    {
        if (mpi::comm().rank() == 0 && mpi::comm().size() == 1) {
            auto node_cell_connectivity_check = std::array<gidx_t, 500>{
                113, 149, 26,  151, 25,  113, 115, 114, 26,  28,  27,  115, 117, 116, 28,  30,  29,  117, 119, 118, 30,
                32,  31,  119, 121, 120, 32,  34,  33,  121, 123, 122, 34,  37,  35,  36,  123, 125, 124, 37,  39,  38,
                125, 127, 126, 39,  41,  40,  127, 129, 128, 41,  43,  42,  129, 130, 43,  45,  44,  130, 45,  48,  46,
                47,  129, 130, 128, 48,  50,  49,  127, 128, 126, 50,  52,  51,  125, 126, 124, 52,  54,  53,  123, 124,
                122, 54,  56,  55,  121, 122, 120, 56,  59,  57,  58,  119, 120, 118, 59,  61,  60,  117, 118, 116, 61,
                63,  62,  115, 116, 114, 63,  65,  64,  113, 114, 65,  67,  66,  151, 25,  153, 1,   26,  25,  27,  1,
                2,   28,  27,  29,  2,   3,   30,  29,  31,  3,   4,   32,  31,  33,  4,   5,   34,  33,  35,  5,   6,
                35,  36,  6,   7,   37,  36,  38,  7,   8,   39,  38,  40,  8,   9,   41,  40,  42,  9,   10,  43,  42,
                44,  10,  11,  45,  44,  46,  11,  12,  46,  47,  12,  13,  48,  47,  49,  13,  14,  50,  49,  51,  14,
                15,  52,  51,  53,  15,  16,  54,  53,  55,  16,  17,  56,  55,  57,  17,  18,  57,  58,  18,  19,  59,
                58,  60,  19,  20,  61,  60,  62,  20,  21,  63,  62,  64,  21,  22,  65,  64,  66,  22,  23,  67,  66,
                68,  23,  24,  153, 1,   155, 69,  1,   2,   69,  71,  70,  2,   3,   71,  73,  72,  3,   4,   73,  75,
                74,  4,   5,   75,  77,  76,  5,   6,   77,  79,  78,  6,   7,   79,  80,  7,   8,   80,  82,  81,  8,
                9,   82,  84,  83,  9,   10,  84,  86,  85,  10,  11,  86,  88,  87,  11,  12,  88,  90,  89,  12,  13,
                90,  91,  13,  14,  91,  93,  92,  14,  15,  93,  95,  94,  15,  16,  95,  97,  96,  16,  17,  97,  99,
                98,  17,  18,  99,  101, 100, 18,  19,  101, 102, 19,  20,  102, 104, 103, 20,  21,  104, 106, 105, 21,
                22,  106, 108, 107, 22,  23,  108, 110, 109, 23,  24,  110, 112, 111, 155, 69,  157, 70,  131, 71,  70,
                72,  131, 133, 132, 73,  72,  74,  133, 135, 134, 75,  74,  76,  135, 137, 136, 77,  76,  78,  137, 139,
                138, 79,  80,  78,  81,  139, 141, 140, 82,  81,  83,  141, 143, 142, 84,  83,  85,  143, 145, 144, 86,
                85,  87,  145, 147, 146, 88,  87,  89,  147, 148, 90,  91,  89,  92,  148, 93,  92,  94,  147, 148, 146,
                95,  94,  96,  145, 146, 144, 97,  96,  98,  143, 144, 142, 99,  98,  100, 141, 142, 140, 101, 102, 100,
                103, 139, 140, 138, 104, 103, 105, 137, 138, 136, 106, 105, 107, 135, 136, 134, 108, 107, 109, 133, 134,
                132, 110, 109, 111, 131, 132, 67,  150, 68,  152, 68,  152, 24,  154, 24,  154, 112, 156, 112, 156, 111,
                158, 149, 149, 151, 153, 153, 155, 157, 157, 150, 150, 152, 154, 154, 156, 158, 158,
            };
            check_node_cell_connectivity(mesh, node_cell_connectivity_check);
        }
    }
    // mpi-size 4
    {
        if (mpi::comm().rank() == 0 && mpi::comm().size() == 4) {
            auto node_cell_connectivity_check = std::array<gidx_t, 310>{
                45, 149, 3,   151, 2,   45, 47, 46,  3,   5,   4,   47,  49,  48,  5,   7,   6,  49, 51, 50,  7,
                9,  8,   51,  53,  52,  9,  11, 10,  53,  55,  54,  11,  14,  12,  13,  55,  57, 56, 14, 16,  15,
                57, 59,  58,  16,  18,  17, 59, 61,  60,  18,  20,  19,  61,  62,  20,  22,  21, 62, 22, 25,  23,
                24, 61,  62,  60,  25,  27, 26, 59,  60,  58,  27,  29,  28,  57,  58,  56,  29, 31, 30, 55,  56,
                54, 31,  33,  32,  53,  54, 52, 33,  36,  34,  35,  51,  52,  50,  36,  38,  37, 49, 50, 48,  38,
                40, 39,  47,  48,  46,  40, 42, 41,  45,  46,  42,  44,  43,  151, 2,   153, 1,  3,  2,  4,   1,
                63, 44,  150, 86,  152, 5,  4,  6,   63,  64,  7,   6,   8,   64,  65,  9,   8,  10, 65, 66,  11,
                10, 12,  66,  67,  12,  13, 67, 68,  14,  13,  15,  68,  69,  16,  15,  17,  69, 70, 18, 17,  19,
                70, 71,  20,  19,  21,  71, 72, 22,  21,  23,  72,  73,  23,  24,  73,  74,  25, 24, 26, 74,  75,
                27, 26,  28,  75,  76,  29, 28, 30,  76,  77,  31,  30,  32,  77,  78,  33,  32, 34, 78, 79,  34,
                35, 79,  80,  36,  35,  37, 80, 81,  38,  37,  39,  81,  82,  40,  39,  41,  82, 83, 42, 41,  43,
                83, 84,  44,  43,  86,  84, 85, 153, 1,   155, 87,  1,   63,  87,  89,  88,  63, 64, 89, 64,  65,
                65, 66,  66,  67,  67,  68, 68, 69,  69,  70,  70,  71,  71,  72,  72,  73,  73, 74, 86, 152, 85,
                74, 75,  75,  76,  76,  77, 77, 78,  78,  79,  79,  80,  80,  81,  81,  82,  82, 83, 83, 84,  84,
                85, 85,  155, 87,  88,  89, 88, 149, 149, 151, 153, 153, 155, 150, 150, 152,
            };
            check_node_cell_connectivity(mesh, node_cell_connectivity_check);
        }

        if (mpi::comm().rank() == 1 && mpi::comm().size() == 4) {
            auto node_cell_connectivity_check = std::array<gidx_t, 190>{
                5,  4,  6,   63,  64,  7,  6,   8,   64,  65,  9,   8,   10,  65,  66,  11,  10,  12,  66,
                67, 12, 13,  67,  68,  14, 13,  15,  68,  69,  16,  15,  17,  69,  70,  18,  17,  19,  70,
                71, 20, 19,  21,  71,  72, 22,  21,  23,  72,  73,  1,   87,  1,   63,  87,  89,  88,  63,
                64, 89, 91,  90,  64,  65, 91,  93,  92,  65,  66,  93,  95,  94,  66,  67,  95,  97,  96,
                67, 68, 97,  98,  68,  69, 98,  100, 99,  69,  70,  100, 102, 101, 70,  71,  102, 104, 103,
                71, 72, 104, 106, 105, 72, 73,  106, 108, 107, 3,   2,   4,   1,   63,  23,  24,  73,  74,
                73, 74, 108, 109, 3,   2,  3,   5,   4,   5,   7,   6,   7,   9,   8,   9,   11,  10,  11,
                14, 12, 13,  14,  16,  15, 16,  18,  17,  18,  20,  19,  20,  22,  21,  22,  23,  24,  2,
                1,  24, 74,  74,  109, 87, 88,  89,  88,  90,  91,  90,  92,  93,  92,  94,  95,  94,  96,
                97, 98, 96,  99,  100, 99, 101, 102, 101, 103, 104, 103, 105, 106, 105, 107, 108, 109, 107,
            };
            check_node_cell_connectivity(mesh, node_cell_connectivity_check);
        }

        if (mpi::comm().rank() == 2 && mpi::comm().size() == 4) {
            auto node_cell_connectivity_check = std::array<gidx_t, 203>{
                23,  24,  73,  74,  25,  24,  26,  74,  75,  27,  26,  28,  75,  76,  29,  28,  30,  76,  77,  31,  30,
                32,  77,  78,  33,  32,  34,  78,  79,  34,  35,  79,  80,  36,  35,  37,  80,  81,  38,  37,  39,  81,
                82,  40,  39,  41,  82,  83,  42,  41,  43,  83,  84,  44,  43,  86,  84,  85,  73,  74,  108, 109, 74,
                75,  109, 111, 110, 75,  76,  111, 113, 112, 76,  77,  113, 115, 114, 77,  78,  115, 117, 116, 78,  79,
                117, 119, 118, 79,  80,  119, 120, 80,  81,  120, 122, 121, 81,  82,  122, 124, 123, 82,  83,  124, 126,
                125, 44,  150, 86,  152, 86,  152, 85,  154, 83,  84,  126, 128, 127, 84,  85,  128, 130, 129, 85,  154,
                130, 156, 25,  23,  24,  25,  27,  26,  27,  29,  28,  29,  31,  30,  31,  33,  32,  33,  36,  34,  35,
                36,  38,  37,  38,  40,  39,  40,  42,  41,  42,  44,  43,  23,  73,  73,  108, 108, 109, 110, 111, 110,
                112, 113, 112, 114, 115, 114, 116, 117, 116, 118, 119, 120, 118, 121, 122, 121, 123, 124, 123, 125, 126,
                125, 127, 128, 127, 129, 130, 156, 129, 150, 150, 152, 154, 154, 156,
            };
            check_node_cell_connectivity(mesh, node_cell_connectivity_check);
        }

        if (mpi::comm().rank() == 3 && mpi::comm().size() == 4) {
            auto node_cell_connectivity_check = std::array<gidx_t, 302>{
                83,  84,  126, 128, 127, 84,  85,  128, 130, 129, 155, 87,  157, 88,  131, 89,  88,  90,  131, 133, 132,
                91,  90,  92,  133, 135, 134, 93,  92,  94,  135, 137, 136, 95,  94,  96,  137, 139, 138, 97,  98,  96,
                99,  139, 141, 140, 100, 99,  101, 141, 143, 142, 102, 101, 103, 143, 145, 144, 104, 103, 105, 145, 147,
                146, 106, 105, 107, 147, 148, 108, 109, 107, 110, 148, 111, 110, 112, 147, 148, 146, 113, 112, 114, 145,
                146, 144, 115, 114, 116, 143, 144, 142, 117, 116, 118, 141, 142, 140, 119, 120, 118, 121, 139, 140, 138,
                122, 121, 123, 137, 138, 136, 124, 123, 125, 135, 136, 134, 126, 125, 127, 133, 134, 132, 128, 127, 129,
                131, 132, 153, 1,   155, 87,  1,   63,  87,  89,  88,  63,  64,  89,  91,  90,  64,  65,  91,  93,  92,
                65,  66,  93,  95,  94,  66,  67,  95,  97,  96,  67,  68,  97,  98,  68,  69,  98,  100, 99,  69,  70,
                100, 102, 101, 70,  71,  102, 104, 103, 71,  72,  104, 106, 105, 72,  73,  106, 108, 107, 73,  74,  108,
                109, 74,  75,  109, 111, 110, 75,  76,  111, 113, 112, 76,  77,  113, 115, 114, 77,  78,  115, 117, 116,
                78,  79,  117, 119, 118, 79,  80,  119, 120, 80,  81,  120, 122, 121, 81,  82,  122, 124, 123, 82,  83,
                124, 126, 125, 85,  154, 130, 156, 130, 156, 129, 158, 153, 1,   1,   63,  63,  64,  64,  65,  65,  66,
                66,  67,  67,  68,  68,  69,  69,  70,  70,  71,  71,  72,  72,  73,  73,  74,  74,  75,  75,  76,  76,
                77,  77,  78,  78,  79,  79,  80,  80,  81,  81,  82,  82,  83,  83,  84,  84,  85,  85,  154, 153, 153,
                155, 157, 157, 154, 154, 156, 158, 158,
            };
            check_node_cell_connectivity(mesh, node_cell_connectivity_check);
        }
    }
}

CASE("test_node2cell_AtlasGrids") {
    std::vector<std::string> gridnames{"N16", "O2", "F2", "L2", "H2"};
    for (auto& gridname : gridnames) {
        SECTION(gridname) {
            auto grid = Grid(gridname);
            Mesh mesh = MeshGenerator(grid.meshgenerator()).generate(grid);
            mesh::actions::build_parallel_fields(mesh);
            mesh::actions::build_periodic_boundaries(mesh);
            mesh::actions::build_halo(mesh, 1);
            mesh::actions::build_node_to_cell_connectivity(mesh);
            for (idx_t jnode = 0; jnode < mesh.nodes().size(); ++jnode) {
                EXPECT(mesh.nodes().cell_connectivity().cols(jnode) > 0);
            }
            for (idx_t jcell = 0; jcell < mesh.cells().size(); ++jcell) {
                EXPECT(mesh.cells().node_connectivity().cols(jcell) > 0);
            }
        }
    }
}

//-----------------------------------------------------------------------------

}  // namespace test
}  // namespace atlas

int main(int argc, char** argv) {
    return atlas::test::run(argc, argv);
}
