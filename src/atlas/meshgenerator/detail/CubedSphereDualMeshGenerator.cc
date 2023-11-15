/*
 * (C) Crown Copyright 2021 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <string>
#include <vector>

#include "eckit/utils/Hash.h"

#include "atlas/functionspace/CubedSphereColumns.h"
#include "atlas/grid/CubedSphereGrid.h"
#include "atlas/grid/Distribution.h"
#include "atlas/grid/Partitioner.h"
#include "atlas/library/config.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/mesh/Nodes.h"
#include "atlas/meshgenerator/MeshGenerator.h"
#include "atlas/meshgenerator/detail/CubedSphereDualMeshGenerator.h"
#include "atlas/meshgenerator/detail/MeshGeneratorFactory.h"
#include "atlas/meshgenerator/detail/cubedsphere/CubedSphereUtility.h"
#include "atlas/parallel/mpi/mpi.h"

#define DEBUG_OUTPUT 0
#define DEBUG_OUTPUT_DETAIL 0

namespace atlas {
namespace meshgenerator {

// -----------------------------------------------------------------------------

CubedSphereDualMeshGenerator::CubedSphereDualMeshGenerator(const eckit::Parametrisation& p) {
    // Get mpi_comm
    std::string mpi_comm = mpi::comm().name();
    p.get("mpi_comm", mpi_comm);
    options.set("mpi_comm", mpi_comm);

    configure_defaults();

    // Get number of partitions.
    size_t nb_parts;
    if (p.get("nb_parts", nb_parts)) {
        options.set("nb_parts", nb_parts);
    }

    // Get this partition.
    int part;
    if (p.get("part", part)) {
        options.set("part", part);
    }

    // Get halo size.
    int halo;
    if (p.get("halo", halo)) {
        options.set("halo", halo);
    }

    // Get partitioner.
    std::string partitioner;
    if (p.get("partitioner", partitioner) && partitioner.size()) {
        options.set("partitioner", partitioner);
    }
}

// -----------------------------------------------------------------------------


void CubedSphereDualMeshGenerator::configure_defaults() {
    auto& comm = mpi::comm(options.getString("mpi_comm"));

    // This option sets number of partitions.
    options.set("nb_parts", comm.size());

    // This option sets the part that will be generated.
    options.set("part", comm.rank());

    // This options sets the number of halo elements around each partition.
    options.set("halo", 0);

    // This options sets the default partitioner.
    options.set<std::string>("partitioner", "cubedsphere");
}

// -----------------------------------------------------------------------------

void CubedSphereDualMeshGenerator::generate(const Grid& grid, Mesh& mesh) const {
    // Get partitioner type and number of partitions from config.
    const idx_t nParts         = static_cast<idx_t>(options.get<size_t>("nb_parts"));
    const std::string partType = options.get<std::string>("partitioner");

    auto partConfig = util::Config{};
    partConfig.set("type", partType);
    partConfig.set("partitions", nParts);
    partConfig.set("mpi_comm",options.getString("mpi_comm"));

    // Use lonlat instead of xy for non cubedsphere partitioner.
    if (partType != "cubedsphere") {
        partConfig.set("coordinates", "lonlat");
    }

    // Set distribution.
    mpi::Scope mpi_scope(options.getString("mpi_comm"));
    const auto partitioner  = grid::Partitioner(partConfig);
    const auto distribution = grid::Distribution(grid, partitioner);

    generate(grid, distribution, mesh);
}

// -----------------------------------------------------------------------------

void CubedSphereDualMeshGenerator::generate(const Grid& grid, const grid::Distribution& distribution,
                                            Mesh& mesh) const {
    // Check for correct grid and need for mesh
    ATLAS_ASSERT(!mesh.generated());

    // Cast grid to cubed sphere grid.
    const auto csGrid = CubedSphereGrid(grid);

    // Check for successful cast.
    if (!csGrid) {
        throw_Exception("CubedSphereDualMeshGenerator can only work with a cubedsphere grid.", Here());
    }

    // Check for correct grid stagger.
    if (csGrid.stagger() != "C") {
        throw_Exception("CubedSphereDualMeshGenerator will only work with a cell-centroid grid.", Here());
    }

    // Clone some grid properties.
    setGrid(mesh, csGrid, distribution);

    mpi::Scope mpi_scope(options.getString("mpi_comm"));
    generate_mesh(csGrid, distribution, mesh);
}

// -----------------------------------------------------------------------------

namespace {

// (i, j) pair
class IJ {
public:
    IJ(idx_t i, idx_t j): i_(i), j_(j) {}
    idx_t i() const { return i_; }
    idx_t j() const { return j_; }
    IJ operator+(const IJ& ij) const { return IJ{i() + ij.i(), j() + ij.j()}; }
    IJ operator-(const IJ& ij) const { return IJ{i() - ij.i(), j() - ij.j()}; }

private:
    idx_t i_{};
    idx_t j_{};
};

// Helper function to copy fields.
template <typename Value, int Rank>
void copyField(const Field& sourceField, Field& targetField) {
    // Assign source field values to target field.
    array::make_view<Value, Rank>(targetField).assign(array::make_view<Value, Rank>(sourceField));
}

// Get the surrounding node (i, j) pairs from a cell (i, j) pair.
std::vector<IJ> getIjNodes(const IJ& ijCell, idx_t N) {
    // Rotate ij 90 degrees anitclockwise about ijPivot.
    auto rotateAnticlockwise = [&](const IJ& ij, const IJ& ijPivot) {
        const auto ijTemp = ij - ijPivot;
        return IJ{-ijTemp.j(), ijTemp.i()} + ijPivot;
    };

    // Rotate ij 90 degrees clockwise about ijPivot.
    auto rotateClockwise = [&](const IJ& ij, const IJ& ijPivot) {
        const auto ijTemp = ij - ijPivot;
        return IJ{ijTemp.j(), -ijTemp.i()} + ijPivot;
    };

    // Set standard surrounding nodes.
    auto ijNodes = std::vector<IJ>{{ijCell.i() - 1, ijCell.j() - 1},
                                   {ijCell.i(), ijCell.j() - 1},
                                   {ijCell.i(), ijCell.j()},
                                   {ijCell.i() - 1, ijCell.j()}};

    // Modify nodes that lie in invalid corners of ij space.
    // Either remove a node to make cell triangular, or rotate two of the
    // nodes out of the forbidden space.

    // Bottom-left corner.
    if (ijCell.i() <= 0 && ijCell.j() <= 0) {
        // Triangle.
        if (ijCell.i() == 0 && ijCell.j() == 0) {
            ijNodes.erase(ijNodes.begin());
        }
        // Quad (i)
        else if (ijCell.i() == 0) {
            ijNodes[0] = rotateClockwise(ijNodes[1], IJ{0, 0});
            ijNodes[3] = rotateClockwise(ijNodes[2], IJ{0, 0});
        }
        else {
            // Quad (ii)
            ijNodes[0] = rotateAnticlockwise(ijNodes[3], IJ{0, 0});
            ijNodes[1] = rotateAnticlockwise(ijNodes[2], IJ{0, 0});
        }
    }
    // Bottom-right corner.
    else if (ijCell.i() >= N && ijCell.j() <= 0) {
        // Triangle.
        if (ijCell.i() == N && ijCell.j() == 0) {
            ijNodes.erase(ijNodes.begin() + 1);
        }
        // Quad (i)
        else if (ijCell.j() == 0) {
            ijNodes[0] = rotateClockwise(ijNodes[3], IJ{N - 1, 0});
            ijNodes[1] = rotateClockwise(ijNodes[2], IJ{N - 1, 0});
        }
        // Quad (ii)
        else {
            ijNodes[1] = rotateAnticlockwise(ijNodes[0], IJ{N - 1, 0});
            ijNodes[2] = rotateAnticlockwise(ijNodes[3], IJ{N - 1, 0});
        }
    }
    // Top-right corner.
    else if (ijCell.i() >= N && ijCell.j() >= N) {
        // Triangle.
        if (ijCell.i() == N && ijCell.j() == N) {
            ijNodes.erase(ijNodes.begin() + 2);
        }
        // Quad (i)
        else if (ijCell.i() == N) {
            ijNodes[1] = rotateClockwise(ijNodes[0], IJ{N - 1, N - 1});
            ijNodes[2] = rotateClockwise(ijNodes[3], IJ{N - 1, N - 1});
        }
        // Quad (ii)
        else {
            ijNodes[2] = rotateAnticlockwise(ijNodes[1], IJ{N - 1, N - 1});
            ijNodes[3] = rotateAnticlockwise(ijNodes[0], IJ{N - 1, N - 1});
        }
    }
    // Top-left corner.
    else if (ijCell.i() <= 0 && ijCell.j() >= N) {
        // Triangle.
        if (ijCell.i() == 0 && ijCell.j() == N) {
            ijNodes.erase(ijNodes.begin() + 3);
        }
        // Quad (i)
        else if (ijCell.j() == N) {
            ijNodes[2] = rotateClockwise(ijNodes[1], IJ{0, N - 1});
            ijNodes[3] = rotateClockwise(ijNodes[0], IJ{0, N - 1});
        }
        // Quad (ii)
        else {
            ijNodes[0] = rotateAnticlockwise(ijNodes[1], IJ{0, N - 1});
            ijNodes[3] = rotateAnticlockwise(ijNodes[2], IJ{0, N - 1});
        }
    }

    return ijNodes;
}

}  // namespace

void CubedSphereDualMeshGenerator::generate_mesh(const CubedSphereGrid& csGrid, const grid::Distribution& distribution,
                                                 Mesh& mesh) const {
    ATLAS_TRACE("CubedSphereDualMeshGenerator::generate");

    using Topology = atlas::mesh::Nodes::Topology;
    using namespace detail::cubedsphere;

    const idx_t N     = csGrid.N();
    const idx_t nHalo = options.get<int>("halo");

    //--------------------------------------------------------------------------
    // Create a cubed-sphere mesh.
    //--------------------------------------------------------------------------

    // Generate cubed sphere primal mesh.
    auto primalOptions = options;
    primalOptions.set("halo", nHalo + 1);
    const auto primalMesh = MeshGenerator("cubedsphere", primalOptions).generate(csGrid, distribution);

    // Generate fucntionspaces for cubed sphere primal mesh.
    const auto primalCellsFunctionSpace = functionspace::CubedSphereCellColumns(primalMesh);
    const auto primalNodesFunctionSpace = functionspace::CubedSphereNodeColumns(primalMesh);
    const auto& primalCells             = primalCellsFunctionSpace.cells();
    const auto& primalNodes             = primalNodesFunctionSpace.nodes();

    //--------------------------------------------------------------------------
    // Set dual mesh nodes (easy part).
    //--------------------------------------------------------------------------

    auto& nodes = mesh.nodes();
    nodes.resize(primalCellsFunctionSpace.size());

    nodes.add(Field("tij", array::make_datatype<idx_t>(), array::make_shape(nodes.size(), 3)));

    // Copy mesh fields to dual mesh.
    copyField<gidx_t, 1>(primalCells.global_index(), nodes.global_index());
    copyField<idx_t, 1>(primalCells.remote_index(), nodes.remote_index());
    copyField<int, 1>(primalCells.partition(), nodes.partition());
    copyField<int, 1>(primalCells.halo(), nodes.halo());
    copyField<int, 1>(primalCells.flags(), nodes.flags());
    copyField<int, 1>(primalCells.field("ghost"), nodes.ghost());
    copyField<double, 2>(primalCells.field("xy"), nodes.xy());
    copyField<double, 2>(primalCells.field("lonlat"), nodes.lonlat());
    copyField<idx_t, 2>(primalCells.field("tij"), nodes.field("tij"));

    // Need to decrement halo by one.
    auto nodesHalo = array::make_view<int, 1>(nodes.halo());

    for (idx_t idx = 0; idx < nodes.size(); ++idx) {
        nodesHalo(idx) = std::max(0, nodesHalo(idx) - 1);
    }

    //--------------------------------------------------------------------------
    // Set dual mesh cells (not so easy part).
    //--------------------------------------------------------------------------

    // Make views to cubed sphere nodes.
    const auto primalNodesHalo = array::make_view<int, 1>(primalNodes.halo());
    const auto primalNodesTij  = array::make_view<idx_t, 2>(primalNodes.field("tij"));

    // Loop over all nodes of primal mesh, excluding outermost halo.
    // Find dual mesh nodes around primal mesh nodes.
    // Note: some halo cells near corners may be incomplete.

    // Set of nodes around a cell.
    struct NodeList {
        std::vector<idx_t> nodes{};  // Node indices.
        bool incomplete{};           // True if nodes are missing.
    };

    auto nodeLists = std::vector<NodeList>{};

    for (idx_t idx = 0; idx < primalNodes.size(); ++idx) {
        // Exclude outer ring of cubed sphere mesh halo.
        if (primalNodesHalo(idx) == nHalo + 1) {
            break;
        }

        nodeLists.emplace_back();
        auto& nodeList = nodeLists.back();

        // Get tij of cell.
        const auto tCell  = primalNodesTij(idx, Coordinates::T);
        const auto ijCell = IJ{primalNodesTij(idx, Coordinates::I), primalNodesTij(idx, Coordinates::J)};

        // Get ij of surrounding nodes.
        auto ijNodes = getIjNodes(ijCell, N);

        // Add indices to nodes vector.
        for (const auto& ijNode : ijNodes) {
            if (primalCellsFunctionSpace.is_valid_index(tCell, ijNode.i(), ijNode.j())) {
                nodeList.nodes.push_back(primalCellsFunctionSpace.index(tCell, ijNode.i(), ijNode.j()));
            }
            else {
                nodeList.incomplete = true;
            }
        }
    }

    // Figure out element types.
    // Set first element type. ( first = type, second = count )
    enum struct ElemType : size_t
    {
        LINE          = 2,
        TRIANGLE      = 3,
        QUADRILATERAL = 4
    };
    auto typeCounts =
        std::vector<std::pair<ElemType, idx_t>>{std::make_pair(static_cast<ElemType>(nodeLists[0].nodes.size()), 1)};

    // Count the number of consecutive lines, triangles or quadtrilaterals in dual mesh.
    // This is an attempt to keep dual mesh cells in the same order as mesh nodes.
    // Otherwise, the halo exchange bookkeeping is invalidated.
    for (size_t idx = 1; idx < nodeLists.size(); ++idx) {
        // Get the element type.
        const auto elemType = static_cast<ElemType>(nodeLists[idx].nodes.size());

        // Increment counter if this elemType is the same as last one
        if (elemType == typeCounts.back().first) {
            ++typeCounts.back().second;
        }
        // Otherwise add a new counter.
        else {
            typeCounts.emplace_back(elemType, 1);
        }
    }

    // Add cells to mesh.
    auto& cells = mesh.cells();
    idx_t nCells{};
    // Loop through type counters.
    for (const auto& typeCount : typeCounts) {
        // Select element type.
        switch (typeCount.first) {
            // Add a block of lines.
            case ElemType::LINE: {
                cells.add(mesh::ElementType::create("Line"), typeCount.second);
                break;
            }
            // Add a block of triangles.
            case ElemType::TRIANGLE: {
                cells.add(mesh::ElementType::create("Triangle"), typeCount.second);
                break;
            }
                // Add a block of quadrilaterals.
            case ElemType::QUADRILATERAL: {
                cells.add(mesh::ElementType::create("Quadrilateral"), typeCount.second);
                break;
            }
            default: {
                ATLAS_THROW_EXCEPTION("Unknown element type with " +
                                      std::to_string(static_cast<size_t>(typeCount.first)) + " nodes.");
                break;
            }
        }
        // Increment the total number of cells.
        nCells += typeCount.second;
    }

    // Add extra fields to cells.
    cells.add(Field("tij", array::make_datatype<idx_t>(), array::make_shape(cells.size(), 3)));
    cells.add(Field("xy", array::make_datatype<double>(), array::make_shape(cells.size(), 2)));
    cells.add(Field("lonlat", array::make_datatype<double>(), array::make_shape(cells.size(), 2)));
    cells.add(Field("ghost", array::make_datatype<int>(), array::make_shape(cells.size())));

    // Copy dual cells fields from nodes.
    copyField<gidx_t, 1>(primalNodes.global_index(), cells.global_index());
    copyField<idx_t, 1>(primalNodes.remote_index(), cells.remote_index());
    copyField<int, 1>(primalNodes.partition(), cells.partition());
    copyField<int, 1>(primalNodes.ghost(), cells.field("ghost"));
    copyField<int, 1>(primalNodes.halo(), cells.halo());
    copyField<int, 1>(primalNodes.flags(), cells.flags());
    copyField<idx_t, 2>(primalNodes.field("tij"), cells.field("tij"));
    copyField<double, 2>(primalNodes.xy(), cells.field("xy"));
    copyField<double, 2>(primalNodes.lonlat(), cells.field("lonlat"));

    // Get view of flags field.
    auto dualCellsFlags = array::make_view<int, 1>(cells.flags());

    // Get node connectivity.
    auto& nodeConnectivity = cells.node_connectivity();

    // Loop over dual mesh cells and set connectivity.
    for (idx_t idx = 0; idx < nCells; ++idx) {
        // Set connectivity.
        nodeConnectivity.set(idx, nodeLists[idx].nodes.data());

        // Set invalid flag if cell is incomplete.
        if (nodeLists[idx].incomplete) {
            Topology::set(dualCellsFlags(idx), Topology::INVALID);
        }
    }

    // Set metadata.
    set_metadata(mesh);

    return;
}

// -----------------------------------------------------------------------------

void CubedSphereDualMeshGenerator::set_metadata(Mesh& mesh) const {
    const auto nHalo = options.get<int>("halo");

    // Set basic halo metadata.
    mesh.metadata().set("halo", nHalo);
    mesh.metadata().set("halo_locked", true);
    mesh.nodes().metadata().set("parallel", true);
    mesh.cells().metadata().set("parallel", true);

    mesh.metadata().set("mpi_comm",options.getString("mpi_comm"));


    // Loop over nodes and count number of halo elements.
    auto nNodes         = std::vector<idx_t>(nHalo + 2, 0);
    const auto nodeHalo = array::make_view<int, 1>(mesh.nodes().halo());
    for (idx_t i = 0; i < mesh.nodes().size(); ++i) {
        ++nNodes[static_cast<size_t>(nodeHalo(i))];
    }
    std::partial_sum(nNodes.begin(), nNodes.end(), nNodes.begin());

    // Set node halo metadata.
    for (size_t i = 0; i < nNodes.size(); ++i) {
        const auto str = "nb_nodes_including_halo[" + std::to_string(i) + "]";
        mesh.metadata().set(str, nNodes[i]);
    }


    // Loop over cells and count number of halo elements.
    auto nCells         = std::vector<std::vector<idx_t>>(mesh.cells().nb_types(), std::vector<idx_t>(nHalo + 1, 0));
    const auto cellHalo = array::make_view<int, 1>(mesh.cells().halo());

    for (idx_t i = 0; i < mesh.cells().nb_types(); ++i) {
        const auto& elems = mesh.cells().elements(i);
        for (idx_t j = elems.begin(); j < elems.end(); ++j) {
            ++nCells[static_cast<size_t>(i)][static_cast<size_t>(cellHalo(j))];
        }
        std::partial_sum(nCells[static_cast<size_t>(i)].begin(), nCells[static_cast<size_t>(i)].end(),
                         nCells[static_cast<size_t>(i)].begin());
    }

    // Set cell halo metadata.
    for (size_t i = 0; i < nCells.size(); ++i) {
        for (size_t j = 0; j < nCells[i].size(); ++j) {
            const auto str = "nb_cells_including_halo[" + std::to_string(i) + "][" + std::to_string(j) + "]";
            mesh.metadata().set(str, nCells[i][j]);
        }
    }
}

// -----------------------------------------------------------------------------

void CubedSphereDualMeshGenerator::hash(eckit::Hash& h) const {
    h.add("CubedSphereDualMeshGenerator");
    options.hash(h);
}

// -----------------------------------------------------------------------------

namespace {
static MeshGeneratorBuilder<CubedSphereDualMeshGenerator> CubedSphereDualMeshGenerator(
    CubedSphereDualMeshGenerator::static_type());
}

// -----------------------------------------------------------------------------

}  // namespace meshgenerator
}  // namespace atlas
