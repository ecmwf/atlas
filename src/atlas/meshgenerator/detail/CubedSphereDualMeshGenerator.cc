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
    if (p.get("partitioner", partitioner)) {
        options.set("partitioner", partitioner);
    }
}

// -----------------------------------------------------------------------------


void CubedSphereDualMeshGenerator::configure_defaults() {
    // This option sets number of partitions.
    options.set("nb_parts", mpi::size());

    // This option sets the part that will be generated.
    options.set("part", mpi::rank());

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

    // Use lonlat instead of xy for non cubedsphere partitioner.
    if (partType != "cubedsphere") {
        partConfig.set("coordinates", "lonlat");
    }

    // Set distribution.
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

    // Enforce compatible halo size.
    if (options.get<idx_t>("halo") != 0) {
        throw_Exception(
            "Halo size CubedSphereDualMeshGenerator is currently "
            "limited to 0.",
            Here());
    }

    // Clone some grid properties.
    setGrid(mesh, csGrid, distribution);

    generate_mesh(csGrid, distribution, mesh);
}

// -----------------------------------------------------------------------------

namespace {

// Helper function to copy fields.
template <typename Value, int Rank>
void copyField(const Field& sourceField, Field& targetField) {
    // Assign source field values to target field.
    array::make_view<Value, Rank>(targetField).assign(array::make_view<Value, Rank>(sourceField));
}

}  // namespace

void CubedSphereDualMeshGenerator::generate_mesh(const CubedSphereGrid& csGrid, const grid::Distribution& distribution,
                                                 Mesh& mesh) const {
    ATLAS_TRACE("CubedSphereDualMeshGenerator::generate");

    using namespace detail::cubedsphere;

    const idx_t N     = csGrid.N();
    const idx_t nHalo = options.get<int>("halo");

    //--------------------------------------------------------------------------
    // Create a cubed-sphere mesh.
    //--------------------------------------------------------------------------

    // Generate cubed sphere mesh.
    auto csOptions = options;
    csOptions.set("halo", nHalo + 1);
    const auto csMesh = MeshGenerator("cubedsphere", csOptions).generate(csGrid, distribution);

    // Generate fucntionspaces cubed sphere mesh.
    const auto csCellsFunctionSpace = functionspace::CubedSphereCellColumns(csMesh);
    const auto csNodesFunctionSpace = functionspace::CubedSphereNodeColumns(csMesh);
    const auto& csCells             = csCellsFunctionSpace.cells();
    const auto& csNodes             = csNodesFunctionSpace.nodes();

    //--------------------------------------------------------------------------
    // Set dual mesh nodes (easy part).
    //--------------------------------------------------------------------------

    auto& nodes = mesh.nodes();
    nodes.resize(csCellsFunctionSpace.size());

    nodes.add(Field("tij", array::make_datatype<idx_t>(), array::make_shape(nodes.size(), 3)));

    // Copy mesh fields to dual mesh.
    copyField<gidx_t, 1>(csCells.global_index(), nodes.global_index());
    copyField<idx_t, 1>(csCells.remote_index(), nodes.remote_index());
    copyField<int, 1>(csCells.partition(), nodes.partition());
    copyField<int, 1>(csCells.halo(), nodes.halo());
    copyField<int, 1>(csCells.flags(), nodes.flags());
    copyField<int, 1>(csCells.field("ghost"), nodes.ghost());
    copyField<double, 2>(csCells.field("xy"), nodes.xy());
    copyField<double, 2>(csCells.field("lonlat"), nodes.lonlat());
    copyField<idx_t, 2>(csCells.field("tij"), nodes.field("tij"));

    // Need to decrement halo by one.
    auto nodesHalo = array::make_view<int, 1>(nodes.halo());

    for (idx_t idx = 0; idx < nodes.size(); ++idx) {
        nodesHalo(idx) = std::max(0, nodesHalo(idx) - 1);
    }

    //--------------------------------------------------------------------------
    // Set dual mesh cells (not so easy part).
    //--------------------------------------------------------------------------

    // Make views to cubed sphere nodes.
    const auto csNodesHalo = array::make_view<int, 1>(csNodes.halo());
    const auto csNodesTij  = array::make_view<idx_t, 2>(csNodes.field("tij"));

    // Figure out element types.
    enum struct ElemType
    {
        tri,
        quad
    };
    const auto getType = [&](idx_t idx) -> ElemType {
        // Nodes of csMesh are cells of dual mesh.
        const idx_t i   = csNodesTij(idx, Coordinates::I);
        const idx_t j   = csNodesTij(idx, Coordinates::J);
        auto cornerCell = (i == 0 && j == 0) || (i == N && j == 0) || (i == N && j == N) || (i == 0 && j == N);
        return cornerCell ? ElemType::tri : ElemType::quad;
    };

    // Set first element type. ( first = type, second = count )
    auto typeCounts = std::vector<std::pair<ElemType, idx_t> >{std::make_pair(getType(0), 1)};

    // Count the number of consecutive triangles and quadtrilaterals in dual mesh.
    // This is an attempt to keep dual mesh cells in the same order as mesh nodes.
    // Otherwise, the halo exchange bookkeeping is invalidated.
    for (idx_t idx = 1; idx < csNodes.size(); ++idx) {
        // Exclude outer ring of cubed sphere mesh halo.
        if (csNodesHalo(idx) == nHalo + 1) {
            break;
        }

        // Get the element type.
        const auto elemType = getType(idx);

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
            // Add a batch of triangles.
            case ElemType::tri: {
                cells.add(new mesh::temporary::Triangle(), typeCount.second);
                break;
            }
            // Add a batch quadrilaterals.
            case ElemType::quad: {
                cells.add(new mesh::temporary::Quadrilateral(), typeCount.second);
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
    copyField<gidx_t, 1>(csNodes.global_index(), cells.global_index());
    copyField<idx_t, 1>(csNodes.remote_index(), cells.remote_index());
    copyField<int, 1>(csNodes.partition(), cells.partition());
    copyField<int, 1>(csNodes.ghost(), cells.field("ghost"));
    copyField<int, 1>(csNodes.halo(), cells.halo());
    copyField<int, 1>(csNodes.flags(), cells.flags());
    copyField<idx_t, 2>(csNodes.field("tij"), cells.field("tij"));
    copyField<double, 2>(csNodes.xy(), cells.field("xy"));
    copyField<double, 2>(csNodes.lonlat(), cells.field("lonlat"));

    // Get node connectivity.
    auto& nodeConnectivity = cells.node_connectivity();

    // Loop over dual mesh cells and set connectivity.
    for (idx_t idx = 0; idx < nCells; ++idx) {
        const idx_t t = csNodesTij(idx, Coordinates::T);
        const idx_t i = csNodesTij(idx, Coordinates::I);
        const idx_t j = csNodesTij(idx, Coordinates::J);


        // Declare vector of surrounding nodes.
        auto cellNodes = std::vector<idx_t>{};

        // Get number of nodes per element.
        const auto nNodes = cells.node_connectivity().row(idx).size();
        // Element is a triangle.
        if (nNodes == 3) {
            // Bottom-left corner.
            if (i == 0 && j == 0) {
                cellNodes = {csCellsFunctionSpace.index(t, 0, -1), csCellsFunctionSpace.index(t, 0, 0),
                             csCellsFunctionSpace.index(t, -1, 0)};
            }
            // Bottom-right corner.
            else if (i == N && j == 0) {
                cellNodes = {csCellsFunctionSpace.index(t, N - 1, -1), csCellsFunctionSpace.index(t, N, 0),
                             csCellsFunctionSpace.index(t, N - 1, 0)};
            }
            // Top-right corner.
            else if (i == N && j == N) {
                cellNodes = {csCellsFunctionSpace.index(t, N - 1, N - 1), csCellsFunctionSpace.index(t, N, N - 1),
                             csCellsFunctionSpace.index(t, N - 1, N)};
            }
            // Top-left corner.
            else {
                cellNodes = {csCellsFunctionSpace.index(t, -1, N - 1), csCellsFunctionSpace.index(t, 0, N - 1),
                             csCellsFunctionSpace.index(t, 0, N)};
            }
        }
        // Element is a quadtrilateral.
        else if (nNodes == 4) {
            cellNodes = {csCellsFunctionSpace.index(t, i - 1, j - 1), csCellsFunctionSpace.index(t, i, j - 1),
                         csCellsFunctionSpace.index(t, i, j), csCellsFunctionSpace.index(t, i - 1, j)};
        }
        // Couldn't determine element type.
        else {
            ATLAS_THROW_EXCEPTION("Could not determine element type for cell " + std::to_string(idx) + ".";);
        }

        nodeConnectivity.set(idx, cellNodes.data());
    }

    // Set metadata.
    mesh.metadata().set("halo", nHalo);
    mesh.metadata().set("halo_locked", true);
    mesh.nodes().metadata().set("parallel", true);
    mesh.cells().metadata().set("parallel", true);

    return;
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
