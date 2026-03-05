/*
 * (C) Crown Copyright 2026 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */


#include "atlas/functionspace/NodeColumns.h"
#include "atlas/functionspace/StructuredColumns.h"
#include "atlas/grid/Distribution.h"
#include "atlas/grid/Grid.h"
#include "atlas/grid/Partitioner.h"
#include "atlas/interpolation/Interpolation.h"
#include "atlas/meshgenerator/MeshGenerator.h"
#include "atlas/redistribution/Redistribution.h"
#include "atlas/util/Config.h"
#include "atlas/util/CoordinateEnums.h"
#include "tests/AtlasTestEnvironment.h"

namespace atlas {
namespace test {

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

// Amplitude for test function/
constexpr double testAmplitude = 5.;

// Convert value to integral.
template <typename Value, IsIntegral<Value> = nullptr>
Value castValue(const double& inVal) {
    constexpr double eps     = 1.e-12;
    Value rounded_minus_half = std::round(inVal - 0.5);
    if (std::abs(rounded_minus_half - (inVal - 0.5)) < eps) {
        return static_cast<Value>(rounded_minus_half);
    }
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

// Class to test functionspace redistribution.
template <typename Value, int Rank>
struct TestRedistribution {
public:
    TestRedistribution(const FunctionSpace& sourceFunctionSpace, const FunctionSpace& targetFunctionSpace, const Redistribution& redistribution):
        sourceFunctionSpace_(sourceFunctionSpace),
        targetFunctionSpace_(targetFunctionSpace),
        redist_(redistribution),
        sourceFieldSet_(sourceFunctionSpace_.createField<Value>(fieldConfig<Rank>())),
        targetFieldSet_(targetFunctionSpace_.createField<Value>(fieldConfig<Rank>())),
        sourceView_(array::make_view<Value, Rank>(sourceFieldSet_[0])),
        targetView_(array::make_view<Value, Rank>(targetFieldSet_[0])) {}

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

CASE("Cubesphere to Gauss") {
    // Target function space (Default Gauss)
    const auto gauss_grid = Grid("O32");
    const auto gauss_functionspace = functionspace::StructuredColumns(gauss_grid, util::Config("halo",2));

    // Construct source function space (Cubespherey Gauss)
    const auto cubedsphere_grid = Grid("CS-LFR-14");
    const auto cubedsphere_functionspace = functionspace::NodeColumns(cubedsphere_grid, util::Config("halo",2));
    atlas::StructuredMeshGenerator mesh_generator;
    const auto cubedsphere_partitioner =
        atlas::grid::MatchingPartitioner(cubedsphere_functionspace.mesh(), atlas::option::type("cubedsphere"));
    
    const auto cubedsphere_to_gauss_distribution = cubedsphere_partitioner.partition(gauss_grid);
    const auto cubedsphere_to_gauss_mesh = mesh_generator.generate(gauss_grid, cubedsphere_to_gauss_distribution);
    const auto cubedsphere_to_gauss_functionspace = atlas::functionspace::NodeColumns(cubedsphere_to_gauss_mesh, util::Config("halo",2));

    // Fail to create redistribution from source to target
    EXPECT_THROWS(atlas::Redistribution(cubedsphere_to_gauss_functionspace, gauss_functionspace));

    // Succeed in creating redistribution from source to target
    const util::Config succeedsRedistributionConfig = util::Config("heterogeneous_redistribution", "true");
    const auto succeedsRedistribution =
        Redistribution(cubedsphere_to_gauss_functionspace, gauss_functionspace, succeedsRedistributionConfig);
    
    // Test redistribution on doubles (rank 1, 2, 3), floats, ints, longs
    auto test1 = TestRedistributionPoints1<double>(cubedsphere_to_gauss_functionspace, gauss_functionspace, succeedsRedistribution);
    auto test2 = TestRedistributionPoints2<double>(cubedsphere_to_gauss_functionspace, gauss_functionspace, succeedsRedistribution);
    auto test3 = TestRedistributionPoints3<double>(cubedsphere_to_gauss_functionspace, gauss_functionspace, succeedsRedistribution);
    auto test4 = TestRedistributionPoints1<float>(cubedsphere_to_gauss_functionspace, gauss_functionspace, succeedsRedistribution);
    auto test5 = TestRedistributionPoints1<int>(cubedsphere_to_gauss_functionspace, gauss_functionspace, succeedsRedistribution);
    auto test6 = TestRedistributionPoints1<long>(cubedsphere_to_gauss_functionspace, gauss_functionspace, succeedsRedistribution);

    test1.execute();
    test2.execute();
    test3.execute();
    test4.execute();
    test5.execute();
    test6.execute();
}

}  // namespace test
}  // namespace atlas

int main(int argc, char** argv) {
    return atlas::test::run(argc, argv);
}

