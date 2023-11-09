/*
 * (C) Crown Copyright 2021 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <cmath>

#include "atlas/field/FieldSet.h"
#include "atlas/functionspace.h"
#include "atlas/functionspace/CubedSphereColumns.h"
#include "atlas/grid.h"
#include "atlas/grid/Partitioner.h"
#include "atlas/meshgenerator.h"
#include "atlas/output/Gmsh.h"
#include "atlas/redistribution/Redistribution.h"
#include "atlas/util/Config.h"
#include "atlas/util/Constants.h"
#include "atlas/util/Unique.h"

#include "atlas/mesh/actions/BuildHalo.h"

#include "tests/AtlasTestEnvironment.h"

namespace atlas {
namespace test {


int mpi_color() {
    static int c = mpi::comm("world").rank()%2;
    return c;
}

struct Fixture {
    Fixture() {
        mpi::comm().split(mpi_color(),"split");
    }
    ~Fixture() {
        if (eckit::mpi::hasComm("split")) {
            eckit::mpi::deleteComm("split");
        }
    }
};


// Set floating point tolerance.
template <typename Value>
Value tolerance() {
    return std::numeric_limits<Value>::epsilon() * 64;
}

// Set field config for different ranks.
template <int Rank>
util::Config fieldConfig();
// Rank 1 config.
template <>
util::Config fieldConfig<1>() {
    return util::Config();
}
// Rank 2 config.
template <>
util::Config fieldConfig<2>() {
    auto config = util::Config();
    config.set("levels", 10);
    return config;
}
// Rank 3 config.
template <>
util::Config fieldConfig<3>() {
    auto config = util::Config();
    config.set("levels", 10);
    config.set("variables", 2);
    return config;
}

// Helper types to distinguish floating point and integral types.
template <typename Value>
using IsIntegral = typename std::enable_if<std::is_integral<Value>::value>::type*;

template <typename Value>
using IsFloatingPoint = typename std::enable_if<std::is_floating_point<Value>::value>::type*;

// Aplitude for test function/
constexpr double testAmplitude = 5.;

// Convert value to integral.
template <typename Value, IsIntegral<Value> = nullptr>
Value castValue(const double& inVal) {
    return static_cast<Value>(std::round(inVal));
}

// Cast value to different float.
template <typename Value, IsFloatingPoint<Value> = nullptr>
Value castValue(const double& inVal) {
    return static_cast<Value>(inVal);
}

// Check integral types are equal.
template <typename Value, IsIntegral<Value> = nullptr>
bool checkValue(const Value& valA, const Value& valB) {
    return valA == valB;
}

// Check floating point types are almost equal.
template <typename Value, IsFloatingPoint<Value> = nullptr>
bool checkValue(const Value& valA, const Value& valB) {
    const auto tol = testAmplitude * tolerance<Value>();
    return std::abs(valA - valB) <= tol;
}

// Define test pattern for nodes.
template <typename Value>
Value testPattern(double lambda, double phi, idx_t level) {
    return castValue<Value>(testAmplitude * std::cos(lambda * (1 + level) * M_PI / 180.) *
                            std::cos(phi * (1 + level) * M_PI / 180.));
}
// Define test pattern for HybridElements.
template <typename Value>
Value testPattern(const mesh::Connectivity::Row& elem, const array::ArrayView<double, 2>& lonLatView) {
    Value f{};
    for (idx_t i = 0; i < elem.size(); ++i) {
        f += testPattern<Value>(lonLatView(elem(i), LON), lonLatView(elem(i), LAT), 2);
    }
    return f / elem.size();
}

// Try and get a mesh from functionspace.
Mesh getMesh(const FunctionSpace& functionSpace) {
    // Try and create a pointer to one of the mesh based functionspaces.
    const auto cellColumns       = functionspace::CellColumns(functionSpace);
    const auto edgeColumns       = functionspace::EdgeColumns(functionSpace);
    const auto nodeColumns       = functionspace::NodeColumns(functionSpace);
    const auto structuredColumns = functionspace::StructuredColumns(functionSpace);

    auto mesh = Mesh();

    if (cellColumns) {
        mesh = cellColumns.mesh();
    }
    else if (edgeColumns) {
        mesh = edgeColumns.mesh();
    }
    else if (nodeColumns) {
        mesh = nodeColumns.mesh();
    }
    else if (structuredColumns) {
        auto mpi_comm = util::Config("mpi_comm",functionSpace.mpi_comm());
        auto grid = structuredColumns.grid();
        auto partitioner = grid::Partitioner(functionSpace.distribution(),mpi_comm);
        mesh = MeshGenerator("structured",mpi_comm).generate(grid, partitioner);
    }
    return mesh;
}

// Try and get cells or edges node-connectivity from functionspace.
const mesh::HybridElements::Connectivity* getConnectivity(const FunctionSpace& functionSpace) {
    // Try and create a pointer to CellColumns or EdgeColumns.
    const auto cellColumns = functionspace::CellColumns(functionSpace);
    const auto edgeColumns = functionspace::EdgeColumns(functionSpace);

    if (cellColumns) {
        return &(cellColumns.cells().node_connectivity());
    }
    else if (edgeColumns) {
        return &(edgeColumns.edges().node_connectivity());
    }
    else {
        return nullptr;
    }
}

// Class to test functionspace redistribution.
template <typename Value, int Rank>
struct TestRedistribution {
public:
    TestRedistribution(const FunctionSpace& sourceFunctionSpace, const FunctionSpace& targetFunctionSpace):
        sourceFunctionSpace_(sourceFunctionSpace),
        targetFunctionSpace_(targetFunctionSpace),
        redist_(sourceFunctionSpace, targetFunctionSpace),
        sourceFieldSet_(sourceFunctionSpace_.createField<Value>(fieldConfig<Rank>())),
        targetFieldSet_(targetFunctionSpace_.createField<Value>(fieldConfig<Rank>())),
        sourceView_(array::make_view<Value, Rank>(sourceFieldSet_[0])),
        targetView_(array::make_view<Value, Rank>(targetFieldSet_[0])) {}

    void outputFields(const std::string& fileStr) const {
        // Try and create a pointer to one of the mesh based functionspaces.
        const auto sourceMesh = getMesh(sourceFunctionSpace_);
        const auto targetMesh = getMesh(targetFunctionSpace_);

        // Output gmsh if both pointers are not null.
        if (sourceMesh && targetMesh) {
            // Set gmsh config.
            const auto gmshConfigXy = util::Config("coordinates", "xy");

            // Set source gmsh object.
            auto sourceGmshXy = output::Gmsh(fileStr + "_source_xy.msh", gmshConfigXy);

            // Write out source mesh and field.
            sourceGmshXy.write(sourceMesh);
            sourceGmshXy.write(sourceFieldSet_, sourceFunctionSpace_);


            // Set target gmsh object.
            auto targetGmshXy = output::Gmsh(fileStr + "_target_xy.msh", gmshConfigXy);

            // Write out source mesh and field.
            targetGmshXy.write(targetMesh);
            targetGmshXy.write(targetFieldSet_, targetFunctionSpace_);
        }
    }

protected:
    FunctionSpace sourceFunctionSpace_;
    FunctionSpace targetFunctionSpace_;

    Redistribution redist_;

    FieldSet sourceFieldSet_;
    FieldSet targetFieldSet_;

    array::ArrayView<Value, Rank> sourceView_;
    array::ArrayView<Value, Rank> targetView_;
};

// Test rank 1 fields with lonlat method.
template <typename Value>
struct TestRedistributionPoints1 : public TestRedistribution<Value, 1> {
    using TestRedistribution<Value, 1>::TestRedistribution;
    void execute() {
        auto sourceLonlatView = array::make_view<double, 2>(this->sourceFunctionSpace_.lonlat());
        auto targetLonlatView = array::make_view<double, 2>(this->targetFunctionSpace_.lonlat());

        // Set source field.
        for (idx_t i = 0; i < this->sourceView_.shape(0); ++i) {
            this->sourceView_(i) = testPattern<Value>(sourceLonlatView(i, 0), sourceLonlatView(i, 1), 0);
        }

        // Perform redistribution.
        this->redist_.execute(this->sourceFieldSet_, this->targetFieldSet_);

        // Perform halo exchange;
        this->targetFunctionSpace_.haloExchange(this->targetFieldSet_);

        // Check target field.
        int nCheck{};
        for (idx_t i = 0; i < this->targetView_.shape(0); ++i) {
            EXPECT(checkValue(this->targetView_(i),
                              testPattern<Value>(targetLonlatView(i, LON), targetLonlatView(i, LAT), 0)));
            ++nCheck;
        }
        const auto& comm = mpi::comm(this->sourceFunctionSpace_.mpi_comm());
        comm.allReduceInPlace(nCheck, eckit::mpi::Operation::SUM);
        Log::debug() << "Checked " << nCheck << " elements." << std::endl;
    }
};

// Test rank 2 fields.
template <typename Value>
struct TestRedistributionPoints2 : public TestRedistribution<Value, 2> {
    using TestRedistribution<Value, 2>::TestRedistribution;
    void execute() {
        auto sourceLonlatView = array::make_view<double, 2>(this->sourceFunctionSpace_.lonlat());
        auto targetLonlatView = array::make_view<double, 2>(this->targetFunctionSpace_.lonlat());

        // Set source field.
        for (idx_t i = 0; i < this->sourceView_.shape(0); ++i) {
            for (idx_t j = 0; j < this->sourceView_.shape(1); ++j) {
                this->sourceView_(i, j) = testPattern<Value>(sourceLonlatView(i, LON), sourceLonlatView(i, LAT), j);
            }
        }

        // Perform redistribution.
        this->redist_.execute(this->sourceFieldSet_, this->targetFieldSet_);

        // Perform halo exchange;
        this->targetFunctionSpace_.haloExchange(this->targetFieldSet_);

        // Check target field.
        int nCheck{};
        for (idx_t i = 0; i < this->targetView_.shape(0); ++i) {
            for (idx_t j = 0; j < this->targetView_.shape(1); ++j) {
                EXPECT(checkValue(this->targetView_(i, j),
                                  testPattern<Value>(targetLonlatView(i, LON), targetLonlatView(i, LAT), j)));
                ++nCheck;
            }
        }
        const auto& comm = mpi::comm(this->sourceFunctionSpace_.mpi_comm());
        comm.allReduceInPlace(nCheck, eckit::mpi::Operation::SUM);
        Log::debug() << "Checked " << nCheck << " elements." << std::endl;
    }
};

// Test rank 3 fields .
template <typename Value>
struct TestRedistributionPoints3 : public TestRedistribution<Value, 3> {
    using TestRedistribution<Value, 3>::TestRedistribution;
    void execute() {
        auto sourceLonlatView = array::make_view<double, 2>(this->sourceFunctionSpace_.lonlat());
        auto targetLonlatView = array::make_view<double, 2>(this->targetFunctionSpace_.lonlat());

        // Set source field.
        for (idx_t i = 0; i < this->sourceView_.shape(0); ++i) {
            for (idx_t j = 0; j < this->sourceView_.shape(1); ++j) {
                this->sourceView_(i, j, 0) = testPattern<Value>(sourceLonlatView(i, LON), sourceLonlatView(i, LAT), j);
                this->sourceView_(i, j, 1) = -testPattern<Value>(sourceLonlatView(i, LON), sourceLonlatView(i, LAT), j);
            }
        }

        // Perform redistribution.
        this->redist_.execute(this->sourceFieldSet_, this->targetFieldSet_);

        // Perform halo exchange;
        this->targetFunctionSpace_.haloExchange(this->targetFieldSet_);


        // Check target field.
        int nCheck{};
        for (idx_t i = 0; i < this->targetView_.shape(0); ++i) {
            for (idx_t j = 0; j < this->targetView_.shape(1); ++j) {
                EXPECT(checkValue(this->targetView_(i, j, 0),
                                  testPattern<Value>(targetLonlatView(i, LON), targetLonlatView(i, LAT), j)));
                ++nCheck;
                EXPECT(checkValue(this->targetView_(i, j, 1),
                                  -testPattern<Value>(targetLonlatView(i, LON), targetLonlatView(i, LAT), j)));
                ++nCheck;
            }
        }
        const auto& comm = mpi::comm(this->sourceFunctionSpace_.mpi_comm());
        comm.allReduceInPlace(nCheck, eckit::mpi::Operation::SUM);
        Log::debug() << "Checked " << nCheck << " elements." << std::endl;
    }
};

// Test CellColumns or EdgeColumns fields.
template <typename Value>
struct TestRedistributionElems : public TestRedistribution<Value, 1> {
    using TestRedistribution<Value, 1>::TestRedistribution;
    void execute() {
        // Get source node connectivity and lonlat view.
        const auto* sourceConnectivity = getConnectivity(this->sourceFunctionSpace_);
        const auto sourceLonLatView = array::make_view<double, 2>(getMesh(this->sourceFunctionSpace_).nodes().lonlat());

        // Get target node connectivity and lonlat view.
        const auto* targetConnectivity = getConnectivity(this->targetFunctionSpace_);
        const auto targetLonLatView = array::make_view<double, 2>(getMesh(this->targetFunctionSpace_).nodes().lonlat());

        // Set source field.
        for (idx_t i = 0; i < this->sourceView_.shape(0); ++i) {
            this->sourceView_(i) = testPattern<Value>(sourceConnectivity->row(i), sourceLonLatView);
        }

        // Perform redistribution.
        this->redist_.execute(this->sourceFieldSet_, this->targetFieldSet_);

        // Perform halo exchange;
        this->targetFunctionSpace_.haloExchange(this->targetFieldSet_);

        // Check target field.
        int nCheck{};
        for (idx_t i = 0; i < this->targetView_.shape(0); ++i) {
            EXPECT(checkValue(this->targetView_(i), testPattern<Value>(targetConnectivity->row(i), targetLonLatView)));
            ++nCheck;
        }
        const auto& comm = mpi::comm(this->sourceFunctionSpace_.mpi_comm());
        comm.allReduceInPlace(nCheck, eckit::mpi::Operation::SUM);
        Log::debug() << "Checked " << nCheck << " elements." << std::endl;
    }
};

CASE("Structured grid") {
    auto grid = atlas::Grid("L24x19");

    // Set mesh config.
    const auto sourceMeshConfig = util::Config("partitioner", "equal_regions");
    const auto targetMeshConfig = util::Config("partitioner", "equal_bands");

    auto sourceMesh = MeshGenerator("structured", sourceMeshConfig).generate(grid);
    auto targetMesh = MeshGenerator("structured", targetMeshConfig).generate(grid);

    SECTION("NodeColumns") {
        const auto sourceFunctionSpace = functionspace::NodeColumns(sourceMesh, util::Config("halo", 2));
        const auto targetFunctionSpace = functionspace::NodeColumns(targetMesh, util::Config("halo", 2));

        // Test double for different ranks.
        auto test1 = TestRedistributionPoints1<double>(sourceFunctionSpace, targetFunctionSpace);
        auto test2 = TestRedistributionPoints2<double>(sourceFunctionSpace, targetFunctionSpace);
        auto test3 = TestRedistributionPoints3<double>(sourceFunctionSpace, targetFunctionSpace);

        // Test float.
        auto test4 = TestRedistributionPoints1<float>(sourceFunctionSpace, targetFunctionSpace);

        // Test int.
        auto test5 = TestRedistributionPoints1<int>(sourceFunctionSpace, targetFunctionSpace);

        // Test long.
        auto test6 = TestRedistributionPoints1<long>(sourceFunctionSpace, targetFunctionSpace);

        test1.execute();
        test2.execute();
        test3.execute();
        test4.execute();
        test5.execute();
        test6.execute();

        test2.outputFields("StructuredGrid_NodeColumns");
    }

    SECTION("CellColumns") {
        // No build_cells_global_idx method implemented in mesh/actions/BuildParallelFields.cc.

        Log::debug() << "Structured Grid Cell Columns currently unsupported." << std::endl;

        //const auto sourceFunctionSpace = functionspace::CellColumns( sourceMesh, util::Config( "halo", 0 ) );
        //const auto targetFunctionSpace = functionspace::CellColumns( targetMesh, util::Config( "halo", 0 ) );

        //auto test = TestRedistributionElems<double>( sourceFunctionSpace, targetFunctionSpace );

        //test.execute();

        //test.outputFields( "StructuredGird_CellColumns");
    }

    SECTION("EdgeColumns") {
        // Note StructuredGrid EdegColumns redistribution only works for halo = 0;

        const auto sourceFunctionSpace = functionspace::EdgeColumns(sourceMesh, util::Config("halo", 0));
        const auto targetFunctionSpace = functionspace::EdgeColumns(targetMesh, util::Config("halo", 0));

        // Test long int.
        auto test = TestRedistributionElems<double>(sourceFunctionSpace, targetFunctionSpace);

        test.execute();

        // EdgeColumns not currently supported by gmsh IO.
    }

    SECTION("Structured Columns") {
        const auto sourceFunctionSpace = functionspace::StructuredColumns(
            grid, grid::Partitioner("equal_regions"), util::Config("halo", 2) | util::Config("periodic_points", true));
        const auto targetFunctionSpace = functionspace::StructuredColumns(
            grid, grid::Partitioner("regular_bands"), util::Config("halo", 2) | util::Config("periodic_points", true));

        auto test = TestRedistributionPoints2<double>(sourceFunctionSpace, targetFunctionSpace);

        test.execute();

        test.outputFields("StructuredGrid_StructuredColumns");
    }

    SECTION("Point Cloud") {
        // Make a point cloud from NodeColumns functionspace.
        const auto sourceFunctionSpace = functionspace::NodeColumns(sourceMesh, util::Config("halo", 0));
        const auto targetFunctionSpace = functionspace::NodeColumns(targetMesh, util::Config("halo", 0));

        // Make a vector of lonlats.
        auto sourceLonLat = std::vector<PointXY>{};
        auto targetLonLat = std::vector<PointXY>{};

        const auto sourceGhostView  = array::make_view<int, 1>(sourceFunctionSpace.ghost());
        const auto sourceLonLatView = array::make_view<double, 2>(sourceFunctionSpace.lonlat());
        const auto targetGhostView  = array::make_view<int, 1>(targetFunctionSpace.ghost());
        const auto targetLonLatView = array::make_view<double, 2>(targetFunctionSpace.lonlat());

        // Add non-ghost lonlats to vector.
        sourceLonLat.reserve(sourceFunctionSpace.size());
        for (idx_t i = 0; i < sourceFunctionSpace.size(); ++i) {
            if (!sourceGhostView(i)) {
                sourceLonLat.emplace_back(sourceLonLatView(i, LON), sourceLonLatView(i, LAT));
            }
        }
        targetLonLat.reserve(targetFunctionSpace.size());
        for (idx_t i = 0; i < targetFunctionSpace.size(); ++i) {
            if (!targetGhostView(i)) {
                targetLonLat.emplace_back(targetLonLatView(i, LON), targetLonLatView(i, LAT));
            }
        }

        // Make point cloud functionspaces.
        const auto sourcePointCloud = functionspace::PointCloud(sourceLonLat);
        const auto targetPointCloud = functionspace::PointCloud(targetLonLat);

        auto test = TestRedistributionPoints2<double>(sourcePointCloud, targetPointCloud);

        test.execute();
    }
}

CASE("Cubed sphere grid") {
    auto grid = atlas::Grid("CS-LFR-C-8");

    // Set mesh config.
    const auto sourceMeshConfig = util::Config("partitioner", "equal_regions") | util::Config("halo", "2");
    const auto targetMeshConfig = util::Config("partitioner", "cubedsphere") | util::Config("halo", "2");

    auto sourceMesh = MeshGenerator("cubedsphere", sourceMeshConfig).generate(grid);
    auto targetMesh = MeshGenerator("cubedsphere", targetMeshConfig).generate(grid);

    SECTION("CubedSphereNodeColumns") {
        const auto sourceFunctionSpace = functionspace::CubedSphereNodeColumns(sourceMesh);
        const auto targetFunctionSpace = functionspace::CubedSphereNodeColumns(targetMesh);

        auto test = TestRedistributionPoints2<double>(sourceFunctionSpace, targetFunctionSpace);

        test.execute();

        test.outputFields("CubedSphere_NodeColumns");
    }

    SECTION("CubedSphereCellColumns") {
        const auto sourceFunctionSpace = functionspace::CubedSphereCellColumns(sourceMesh);
        const auto targetFunctionSpace = functionspace::CubedSphereCellColumns(targetMesh);

        auto test = TestRedistributionElems<double>(sourceFunctionSpace, targetFunctionSpace);

        test.execute();

        test.outputFields("CubedSphere_CellColumns");
    }
}

CASE("Cubed sphere dual grid") {
    auto grid = atlas::Grid("CS-LFR-C-8");

    // Set mesh config.
    const auto sourceMeshConfig = util::Config("partitioner", "equal_regions") | util::Config("halo", "0");
    const auto targetMeshConfig = util::Config("partitioner", "cubedsphere") | util::Config("halo", "0");

    auto sourceMesh = MeshGenerator("cubedsphere_dual", sourceMeshConfig).generate(grid);
    auto targetMesh = MeshGenerator("cubedsphere_dual", targetMeshConfig).generate(grid);

    SECTION("CubedSphereDualNodeColumns") {
        const auto sourceFunctionSpace = functionspace::CubedSphereNodeColumns(sourceMesh);
        const auto targetFunctionSpace = functionspace::CubedSphereNodeColumns(targetMesh);

        auto test = TestRedistributionPoints2<double>(sourceFunctionSpace, targetFunctionSpace);

        test.execute();

        test.outputFields("CubedSphereDual_NodeColumns");
    }

    SECTION("CubedSphereDualCellColumns") {
        const auto sourceFunctionSpace = functionspace::CubedSphereCellColumns(sourceMesh);
        const auto targetFunctionSpace = functionspace::CubedSphereCellColumns(targetMesh);

        auto test = TestRedistributionElems<double>(sourceFunctionSpace, targetFunctionSpace);

        test.execute();

        test.outputFields("CubedSphereDual_CellColumns");
    }
}

CASE("Structured grid with split comms") {
    Fixture fixture;

    auto grid = mpi_color() == 0 ? atlas::Grid("L24x13") : atlas::Grid("O16");
    auto mpi_comm = util::Config("mpi_comm","split");

    // auto grid = atlas::Grid("O48");
    // auto mpi_comm = util::Config("mpi_comm","world");

    // Set mesh config.
    const auto sourceMeshConfig = util::Config("partitioner", "equal_regions") | mpi_comm;
    const auto targetMeshConfig = util::Config("partitioner", "equal_bands") | mpi_comm;

    auto sourceMesh = MeshGenerator("structured", sourceMeshConfig).generate(grid);
    auto targetMesh = MeshGenerator("structured", targetMeshConfig).generate(grid);

    SECTION("NodeColumns") {
        const auto sourceFunctionSpace = functionspace::NodeColumns(sourceMesh, util::Config("halo", 2));
        const auto targetFunctionSpace = functionspace::NodeColumns(targetMesh, util::Config("halo", 2));

        // Test double for different ranks.
        auto test1 = TestRedistributionPoints1<double>(sourceFunctionSpace, targetFunctionSpace);
        auto test2 = TestRedistributionPoints2<double>(sourceFunctionSpace, targetFunctionSpace);
        auto test3 = TestRedistributionPoints3<double>(sourceFunctionSpace, targetFunctionSpace);

        // Test float.
        auto test4 = TestRedistributionPoints1<float>(sourceFunctionSpace, targetFunctionSpace);

        // Test int.
        auto test5 = TestRedistributionPoints1<int>(sourceFunctionSpace, targetFunctionSpace);

        // Test long.
        auto test6 = TestRedistributionPoints1<long>(sourceFunctionSpace, targetFunctionSpace);

        test1.execute();
        test2.execute();
        test3.execute();
        test4.execute();
        test5.execute();
        test6.execute();

        test2.outputFields("StructuredGrid_NodeColumns_"+std::to_string(mpi_color()));
    }
    SECTION("CellColumns") {
        // No build_cells_global_idx method implemented in mesh/actions/BuildParallelFields.cc.

        Log::debug() << "Structured Grid Cell Columns currently unsupported." << std::endl;

        //const auto sourceFunctionSpace = functionspace::CellColumns( sourceMesh, util::Config( "halo", 0 ) );
        //const auto targetFunctionSpace = functionspace::CellColumns( targetMesh, util::Config( "halo", 0 ) );

        //auto test = TestRedistributionElems<double>( sourceFunctionSpace, targetFunctionSpace );

        //test.execute();

        //test.outputFields( "StructuredGird_CellColumns_"+std::to_string(mpi_color()));
    }
    SECTION("EdgeColumns") {
        // Note StructuredGrid EdegColumns redistribution only works for halo = 0;

        const auto sourceFunctionSpace = functionspace::EdgeColumns(sourceMesh, util::Config("halo", 0));
        const auto targetFunctionSpace = functionspace::EdgeColumns(targetMesh, util::Config("halo", 0));

        // Test long int.
        auto test = TestRedistributionElems<double>(sourceFunctionSpace, targetFunctionSpace);

        test.execute();

        // EdgeColumns not currently supported by gmsh IO.
    }
    SECTION("Structured Columns") {
        const auto sourceFunctionSpace = functionspace::StructuredColumns(
            grid, grid::Partitioner("equal_regions",mpi_comm), util::Config("halo", 2) | util::Config("periodic_points", true) | mpi_comm);
        const auto targetFunctionSpace = functionspace::StructuredColumns(
            grid, grid::Partitioner("regular_bands",mpi_comm), util::Config("halo", 2) | util::Config("periodic_points", true) | mpi_comm);

        auto test = TestRedistributionPoints2<double>(sourceFunctionSpace, targetFunctionSpace);

        test.execute();

        test.outputFields("StructuredGrid_StructuredColumns_"+std::to_string(mpi_color()));
    }
    SECTION("Point Cloud") {
        // Make a point cloud from NodeColumns functionspace.
        const auto sourceFunctionSpace = functionspace::NodeColumns(sourceMesh, util::Config("halo", 0));
        const auto targetFunctionSpace = functionspace::NodeColumns(targetMesh, util::Config("halo", 0));

        // Make a vector of lonlats.
        auto sourceLonLat = std::vector<PointXY>{};
        auto targetLonLat = std::vector<PointXY>{};

        const auto sourceGhostView  = array::make_view<int, 1>(sourceFunctionSpace.ghost());
        const auto sourceLonLatView = array::make_view<double, 2>(sourceFunctionSpace.lonlat());
        const auto targetGhostView  = array::make_view<int, 1>(targetFunctionSpace.ghost());
        const auto targetLonLatView = array::make_view<double, 2>(targetFunctionSpace.lonlat());

        // Add non-ghost lonlats to vector.
        sourceLonLat.reserve(sourceFunctionSpace.size());
        for (idx_t i = 0; i < sourceFunctionSpace.size(); ++i) {
            if (!sourceGhostView(i)) {
                sourceLonLat.emplace_back(sourceLonLatView(i, LON), sourceLonLatView(i, LAT));
            }
        }
        targetLonLat.reserve(targetFunctionSpace.size());
        for (idx_t i = 0; i < targetFunctionSpace.size(); ++i) {
            if (!targetGhostView(i)) {
                targetLonLat.emplace_back(targetLonLatView(i, LON), targetLonLatView(i, LAT));
            }
        }

        // Make point cloud functionspaces.
        const auto sourcePointCloud = functionspace::PointCloud(sourceLonLat);
        const auto targetPointCloud = functionspace::PointCloud(targetLonLat);

        auto test = TestRedistributionPoints2<double>(sourcePointCloud, targetPointCloud);

        test.execute();
    }
}

CASE("Cubed sphere grid with split comms") {
    Fixture fixture;

    auto grid = mpi_color() == 0 ? atlas::Grid("CS-LFR-C-8") : atlas::Grid("CS-LFR-C-16");
    auto mpi_comm = util::Config("mpi_comm","split");

    // Set mesh config.
    const auto sourceMeshConfig = util::Config("partitioner", "equal_regions") | util::Config("halo", "2") | mpi_comm;
    const auto targetMeshConfig = util::Config("partitioner", "cubedsphere") | util::Config("halo", "2") | mpi_comm;

    auto sourceMesh = MeshGenerator("cubedsphere", sourceMeshConfig).generate(grid);
    auto targetMesh = MeshGenerator("cubedsphere", targetMeshConfig).generate(grid);

    SECTION("CubedSphereNodeColumns") {
        const auto sourceFunctionSpace = functionspace::CubedSphereNodeColumns(sourceMesh);
        const auto targetFunctionSpace = functionspace::CubedSphereNodeColumns(targetMesh);

        EXPECT_EQ( sourceFunctionSpace.mpi_comm(), "split" );

        auto test = TestRedistributionPoints2<double>(sourceFunctionSpace, targetFunctionSpace);

        test.execute();

        test.outputFields("CubedSphere_NodeColumns_"+std::to_string(mpi_color()));
    }

    SECTION("CubedSphereCellColumns") {
        const auto sourceFunctionSpace = functionspace::CubedSphereCellColumns(sourceMesh);
        const auto targetFunctionSpace = functionspace::CubedSphereCellColumns(targetMesh);

        EXPECT_EQ( sourceFunctionSpace.mpi_comm(), "split" );

        auto test = TestRedistributionElems<double>(sourceFunctionSpace, targetFunctionSpace);

        test.execute();

        test.outputFields("CubedSphere_CellColumns_"+std::to_string(mpi_color()));
    }
}

CASE("Cubed sphere dual grid with split comms") {
    Fixture fixture;

    auto grid = mpi_color() == 0 ? atlas::Grid("CS-LFR-C-8") : atlas::Grid("CS-LFR-C-16");
    auto mpi_comm = util::Config("mpi_comm","split");

    // Set mesh config.
    const auto sourceMeshConfig = util::Config("partitioner", "equal_regions") | util::Config("halo", "0") | mpi_comm;
    const auto targetMeshConfig = util::Config("partitioner", "cubedsphere") | util::Config("halo", "0") | mpi_comm;

    auto sourceMesh = MeshGenerator("cubedsphere_dual", sourceMeshConfig).generate(grid);
    auto targetMesh = MeshGenerator("cubedsphere_dual", targetMeshConfig).generate(grid);

    EXPECT_EQ( sourceMesh.mpi_comm(), "split" );

    SECTION("CubedSphereDualNodeColumns") {
        const auto sourceFunctionSpace = functionspace::CubedSphereNodeColumns(sourceMesh);
        const auto targetFunctionSpace = functionspace::CubedSphereNodeColumns(targetMesh);

        EXPECT_EQ( sourceFunctionSpace.mpi_comm(), "split" );


        auto test = TestRedistributionPoints2<double>(sourceFunctionSpace, targetFunctionSpace);

        test.execute();

        test.outputFields("CubedSphereDual_NodeColumns_"+std::to_string(mpi_color()));
    }

    SECTION("CubedSphereDualCellColumns") {
        const auto sourceFunctionSpace = functionspace::CubedSphereCellColumns(sourceMesh);
        const auto targetFunctionSpace = functionspace::CubedSphereCellColumns(targetMesh);

        EXPECT_EQ( sourceFunctionSpace.mpi_comm(), "split" );

        auto test = TestRedistributionElems<double>(sourceFunctionSpace, targetFunctionSpace);

        test.execute();

        test.outputFields("CubedSphereDual_CellColumns_"+std::to_string(mpi_color()));
    }
}

}  // namespace test
}  // namespace atlas


int main(int argc, char** argv) {
    return atlas::test::run(argc, argv);
}
