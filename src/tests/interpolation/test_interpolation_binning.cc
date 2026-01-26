/*
 * (C) Crown Copyright 2024 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */


#include <string>
#include <utility>

#include "atlas/field/FieldSet.h"
#include "atlas/functionspace/CubedSphereColumns.h"
#include "atlas/functionspace/FunctionSpace.h"
#include "atlas/functionspace/NodeColumns.h"
#include "atlas/functionspace/StructuredColumns.h"
#include "atlas/interpolation.h"
#include "atlas/meshgenerator/MeshGenerator.h"
#include "atlas/option/Options.h"
#include "atlas/output/Gmsh.h"
#include "atlas/runtime/Exception.h"
#include "atlas/runtime/Log.h"
#include "atlas/util/Config.h"
#include "atlas/util/CoordinateEnums.h"
#include "atlas/util/function/VortexRollup.h"

#include "tests/AtlasTestEnvironment.h"


namespace atlas {
namespace test {


/// function to generate a synthetic field (vortex field)
Field getVortexField(const FunctionSpace& functionSpace) {
    const auto lonlat = array::make_view<double, 2>(functionSpace.lonlat());

    auto vField     = functionSpace.createField<double>(option::name{"vortex field"});
    auto vFieldView = array::make_view<double, 1>(vField);
    for (idx_t idx = 0; idx < vField.shape(0); ++idx) {
        vFieldView(idx) = util::function::vortex_rollup(lonlat(idx, atlas::LON), lonlat(idx, atlas::LAT), 1.);
    }

    return vField;
}

/// Generate and the difference between a target field and the analytic solution
Field getErrorField(const Field& targetField) {
    auto errorField               = targetField.functionspace().createField<double>(option::name("error field"));
    const auto referenceField     = getVortexField(targetField.functionspace());
    auto errorFieldView           = array::make_view<double, 1>(errorField);
    const auto referenceFieldView = array::make_view<double, 1>(referenceField);
    const auto targetFieldView    = array::make_view<double, 1>(targetField);
    for (idx_t idx = 0; idx < errorField.shape(0); ++idx) {
        errorFieldView(idx) = targetFieldView(idx) - referenceFieldView(idx);
    }
    return errorField;
}

/// function to write a field set in a Gmsh file
void makeGmshOutput(const std::string& fileNamePrefix, const FieldSet& fields) {
    const auto& fs = fields[0].functionspace();

    const auto mesh = [&]() {
        if (const auto nc = functionspace::NodeColumns(fs)) {
            return nc.mesh();
        }
        if (const auto sc = functionspace::StructuredColumns(fs)) {
            return MeshGenerator{"structured", option::halo(sc.halo())}.generate(sc.grid());
        }
        throw_Exception("Unsupported function space type for Gmsh output: " + fs.type());
    }();

    const auto gmshConfig =
        util::Config("coordinates", "xyz") | util::Config("ghost", true) | util::Config("info", true);
    const auto gmsh = output::Gmsh(fileNamePrefix, gmshConfig);
    gmsh.write(mesh);
    gmsh.write(fields, fs);
}


/// function to carry out a dot product
double dotProd(const Field& fieldA, const Field& fieldB) {
    double dprod{};

    const auto field01_view = array::make_view<double, 1>(fieldA);
    const auto field02_view = array::make_view<double, 1>(fieldB);

    for (idx_t i = 0; i < field01_view.shape(0); ++i) {
        dprod += field01_view(i) * field02_view(i);
    }
    eckit::mpi::comm().allReduceInPlace(dprod, eckit::mpi::Operation::SUM);

    return dprod;
}


void regriddingTest(const util::Config& config) {
    const auto sourceGridName = config.getString("source_grid");
    const auto targetGridName = config.getString("target_grid");
    const auto sourceGrid     = Grid{sourceGridName};
    const auto targetGrid     = Grid{targetGridName};

    const auto [sourceFunctionSpace, targetFunctionSpace] = [&]() -> std::pair<FunctionSpace, FunctionSpace> {
        const auto fSpaceName = config.getString("functionspace");
        const size_t haloSize = config.getInt("halo");
        if (fSpaceName == "NodeColumns") {
            const auto meshGenerator = MeshGenerator{config.getString("mesh_generator"), option::halo{haloSize}};

            const auto sourceMesh = meshGenerator.generate(sourceGrid);
            const auto targetMesh = meshGenerator.generate(targetGrid);

            return {functionspace::NodeColumns{sourceMesh, option::halo{haloSize}},
                    functionspace::NodeColumns{targetMesh}};
        }
        if (fSpaceName == "StructuredColumns") {
            return {functionspace::StructuredColumns{sourceGrid, option::halo{haloSize}},
                    functionspace::StructuredColumns{targetGrid, option::halo{3}}};
        }
        throw_Exception("Unsupported functionspace type: " + fSpaceName);
    }();

    const auto sourceField = getVortexField(sourceFunctionSpace);
    auto targetField       = targetFunctionSpace.createField<double>(option::name("field_target"));

    const auto binningScheme = option::type{"binning"} | util::Config{"scheme", config.getSubConfiguration("scheme")};
    auto binning             = Interpolation{binningScheme, sourceFunctionSpace, targetFunctionSpace};

    sourceField.haloExchange();
    binning.execute(sourceField, targetField);
    targetField.haloExchange();

    auto targetFieldSet = FieldSet{};
    targetFieldSet.add(targetField);
    targetFieldSet.add(getErrorField(targetField));

    auto sourceFieldSet = FieldSet{};
    sourceFieldSet.add(sourceField);

    const auto binningType = binningScheme.getString("scheme.type");
    makeGmshOutput(sourceGridName + "_to_" + targetGridName + "_" + binningType + "_source.msh", sourceFieldSet);
    makeGmshOutput(sourceGridName + "_to_" + targetGridName + "_" + binningType + "_target.msh", targetFieldSet);
}

CASE("Regridding from high to low resolution: cubed sphere, bilinear") {
    const auto config = util::Config{"source_grid", "CS-LFR-112"} | util::Config{"target_grid", "CS-LFR-28"} |
                        util::Config{"functionspace", "NodeColumns"} | util::Config{"halo", 1} |
                        util::Config{"mesh_generator", "cubedsphere_dual"} |
                        util::Config{"scheme", option::type{"cubedsphere-bilinear"}};

    regriddingTest(config);
}

CASE("Regridding from high to low resolution: cubed sphere, nearest neighbour") {
    const auto config = util::Config{"source_grid", "CS-LFR-112"} | util::Config{"target_grid", "CS-LFR-28"} |
                        util::Config{"functionspace", "NodeColumns"} | util::Config{"halo", 0} |
                        util::Config{"mesh_generator", "cubedsphere_dual"} |
                        util::Config{"scheme", option::type{"nearest-neighbour"}};

    regriddingTest(config);
}

CASE("Regridding from high to low resolution: cubed sphere, finite element") {
    const auto config = util::Config{"source_grid", "CS-LFR-112"} | util::Config{"target_grid", "CS-LFR-28"} |
                        util::Config{"functionspace", "NodeColumns"} | util::Config{"halo", 1} |
                        util::Config{"mesh_generator", "cubedsphere_dual"} |
                        util::Config{"scheme", option::type{"finite-element"}};

    regriddingTest(config);
}

CASE("Regridding from high to low resolution: gaussian, structured bilinear") {
    const auto config = util::Config{"source_grid", "O48"} | util::Config{"target_grid", "O24"} |
                        util::Config{"functionspace", "StructuredColumns"} | util::Config{"halo", 3} |
                        util::Config{"scheme", option::type{"structured-bilinear"}};

    regriddingTest(config);
}

CASE("plot binning kernel") {
    const auto sourceGrid = Grid{"CS-LFR-20"};
    const auto targetGrid = Grid{"CS-LFR-5"};

    const auto haloOption = option::halo(1);

    const auto sourceMesh = MeshGenerator("cubedsphere_dual", haloOption).generate(sourceGrid);
    const auto targetMesh = MeshGenerator("cubedsphere_dual").generate(targetGrid);

    const auto sourceFunctionSpace = functionspace::NodeColumns{sourceMesh, haloOption};
    const auto targetFunctionSpace = functionspace::NodeColumns{targetMesh};

    const auto binningScheme = option::type{"binning"} | util::Config{"scheme", option::type{"cubedsphere-bilinear"}};
    auto binning             = Interpolation{binningScheme, sourceFunctionSpace, targetFunctionSpace};

    const auto binningMatrixStorage = interpolation::MatrixCache(binning).matrix();
    const auto binningMatrix        = linalg::make_non_owning_eckit_sparse_matrix(binningMatrixStorage);

    const auto targetGhostView = array::make_view<int, 1>(targetFunctionSpace.ghost());
    const auto targetRidxView  = array::make_indexview<idx_t, 1>(targetFunctionSpace.remote_index());
    const auto targetPartView  = array::make_view<int, 1>(targetFunctionSpace.partition());

    const auto targetRemoteIndices = std::vector{0, 2, 12};
    const auto targetPartition     = 0;

    auto kernelField     = sourceFunctionSpace.createField<double>(option::name("kernel_field"));
    auto kernelFieldView = array::make_view<double, 1>(kernelField);
    kernelFieldView.assign(0.);

    for (auto matrixItr = binningMatrix.begin(); matrixItr != binningMatrix.end(); ++matrixItr) {
        const auto rowIdx = matrixItr.row();
        const auto colIdx = matrixItr.col();

        if (targetGhostView(rowIdx)) {
            continue;
        }

        for (const auto targetRemoteIndex : targetRemoteIndices) {
            if (targetRidxView(rowIdx) == targetRemoteIndex && targetPartView(rowIdx) == targetPartition) {
                kernelFieldView(colIdx) += *matrixItr;
            }
        }
    }
    // Force halo weights into owned region of field.
    sourceFunctionSpace.adjointHaloExchange(kernelField);
    sourceFunctionSpace.haloExchange(kernelField);

    makeGmshOutput("binning_kernel.msh", kernelField);
}

/// test to carry out the 'dot-product' test for the rigridding from
/// 'high' to 'low' resolution, for a given type of grid (CS-LFR)
///
CASE("dot-product test for the rigridding from high to low resolution; grid type: Cubed Sphere") {
    // source grid (high res.)
    const auto sourceGrid          = Grid("CS-LFR-100");
    const auto sourceMesh          = MeshGenerator("cubedsphere_dual").generate(sourceGrid);
    const auto sourceFunctionSpace = functionspace::NodeColumns(sourceMesh);

    // target grid (low res.)
    const auto targetGrid          = Grid("CS-LFR-50");
    const auto targetMesh          = MeshGenerator("cubedsphere_dual").generate(targetGrid);
    const auto targetFunctionSpace = functionspace::NodeColumns(targetMesh);

    // source field
    const auto sourceField = getVortexField(sourceFunctionSpace);

    auto sourceFieldSet = FieldSet{};
    sourceFieldSet.add(sourceField);

    // target field
    auto targetField = targetFunctionSpace.createField<double>(option::name("field_01_t"));

    auto targetFieldSet = FieldSet{};
    targetFieldSet.add(targetField);

    const auto scheme = util::Config("type", "binning") | util::Config("scheme", option::type("cubedsphere-bilinear")) |
                        util::Config("adjoint", true);

    Interpolation binning(scheme, sourceFunctionSpace, targetFunctionSpace);

    // performing the regridding from high to low resolution
    binning.execute(sourceFieldSet, targetFieldSet);


    targetFieldSet["field_01_t"].haloExchange();

    // target field (adjoint)
    auto targetFieldStar = targetFunctionSpace.createField<double>(option::name("field_01_ad_t"));
    array::make_view<double, 1>(targetFieldStar).assign(array::make_view<double, 1>(targetField));
    targetFieldStar.adjointHaloExchange();

    auto targetFieldSetStar = FieldSet{};
    targetFieldSetStar.add(targetFieldStar);

    // source field (adjoint)
    auto sourceFieldStar = sourceFunctionSpace.createField<double>(option::name("field_01_ad_s"));
    array::make_view<double, 1>(sourceFieldStar).assign(0.);

    auto sourceFieldSetStar = FieldSet{};
    sourceFieldSetStar.add(sourceFieldStar);

    // performing adjoint operation
    binning.execute_adjoint(sourceFieldSetStar, targetFieldSetStar);


    const auto targetDotTarget     = dotProd(targetField, targetField);
    const auto sourceDotSourceStar = dotProd(sourceField, sourceFieldStar);

    double scaled_diff = std::abs(targetDotTarget - sourceDotSourceStar) / std::abs(targetDotTarget);

    // carrrying out a dot-product test ...
    Log::info() << "\n- dot-product test:\n"
                << "(Ax) . (Ax) = " << targetDotTarget << "; "
                << "x . (A^t A x) = " << sourceDotSourceStar << "; "
                << "scaled difference = " << scaled_diff << "\n"
                << std::endl;

    EXPECT(scaled_diff < 1e-12);
}

}  // namespace test
}  // namespace atlas


//--

int main(int argc, char** argv) {
    return atlas::test::run(argc, argv);
}
