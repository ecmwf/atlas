/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "eckit/types/Types.h"

#include "atlas/array/ArrayView.h"
#include "atlas/array/MakeView.h"
#include "atlas/field/Field.h"
#include "atlas/field/FieldSet.h"
#include "atlas/functionspace/BlockStructuredColumns.h"
#include "atlas/functionspace/CellColumns.h"
#include "atlas/functionspace/EdgeColumns.h"
#include "atlas/functionspace/NodeColumns.h"
#include "atlas/functionspace/PointCloud.h"
#include "atlas/functionspace/Spectral.h"
#include "atlas/functionspace/StructuredColumns.h"
#include "atlas/grid/Grid.h"
#include "atlas/grid/StructuredGrid.h"
#include "atlas/grid/UnstructuredGrid.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/mesh/MeshBuilder.h"
#include "atlas/mesh/Nodes.h"
#include "atlas/meshgenerator.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/trans/Trans.h"

#include "tests/AtlasTestEnvironment.h"

using namespace eckit;
using namespace atlas::functionspace;
using namespace atlas::util;

namespace atlas {
namespace test {

//-----------------------------------------------------------------------------
CASE("test_functionspace_NodeColumns_no_halo") {
    Grid grid("O8");
    Mesh mesh = StructuredMeshGenerator().generate(grid);
    functionspace::NodeColumns nodes_fs(mesh);
    Field field(nodes_fs.createField<int>());
    array::ArrayView<int, 1> value = array::make_view<int, 1>(field);
    array::ArrayView<int, 1> ghost = array::make_view<int, 1>(mesh.nodes().ghost());
    const size_t nb_nodes          = mesh.nodes().size();
    for (size_t j = 0; j < nb_nodes; ++j) {
        if (ghost(j)) {
            value(j) = -1;
        }
        else {
            value(j) = 1;
        }
    }
    nodes_fs.haloExchange(field);
    for (size_t j = 0; j < nb_nodes; ++j) {
        EXPECT(value(j) == 1);
    }
}

CASE("test_functionspace_NodeColumns") {
    ReducedGaussianGrid grid({4, 8, 8, 4});

    StructuredMeshGenerator generator;
    // generator.options.set("3d",true);
    Mesh mesh = generator.generate(grid);

    // grid.reset();

    idx_t nb_levels = 10;

    functionspace::NodeColumns nodes_fs(mesh, option::halo(1) | option::levels(nb_levels));
    // NodesColumnFunctionSpace columns_fs("columns",mesh,nb_levels,Halo(1));

    // EXPECT( nodes_fs.nb_nodes() == columns_fs.nb_nodes() );
    // EXPECT( columns_fs.nb_levels() == 10 );

    Field surface_scalar_field = nodes_fs.createField<double>(option::name("scalar") | option::levels(false));
    Field surface_vector_field =
        nodes_fs.createField<double>(option::name("vector") | option::levels(false) | option::variables(2));
    Field surface_tensor_field =
        nodes_fs.createField<double>(option::name("tensor") | option::levels(false) | option::variables(2 * 2));

    EXPECT(surface_scalar_field.name() == std::string("scalar"));
    EXPECT(surface_vector_field.name() == std::string("vector"));
    EXPECT(surface_tensor_field.name() == std::string("tensor"));

    EXPECT(surface_scalar_field.size() == nodes_fs.nb_nodes());
    EXPECT(surface_vector_field.size() == nodes_fs.nb_nodes() * 2);
    EXPECT(surface_tensor_field.size() == nodes_fs.nb_nodes() * 2 * 2);

    EXPECT(surface_scalar_field.rank() == 1);
    EXPECT(surface_vector_field.rank() == 2);
    EXPECT(surface_tensor_field.rank() == 2);

    EXPECT(surface_scalar_field.levels() == 0);
    EXPECT(surface_vector_field.levels() == 0);
    EXPECT(surface_tensor_field.levels() == 0);

    array::ArrayView<double, 1> surface_scalar = array::make_view<double, 1>(surface_scalar_field);
    array::ArrayView<double, 2> surface_vector = array::make_view<double, 2>(surface_vector_field);
    array::ArrayView<double, 2> surface_tensor = array::make_view<double, 2>(surface_tensor_field);

    EXPECT(surface_scalar.shape(0) == nodes_fs.nb_nodes());
    EXPECT(surface_vector.shape(0) == nodes_fs.nb_nodes());
    EXPECT(surface_tensor.shape(0) == nodes_fs.nb_nodes());
    EXPECT(surface_vector.shape(1) == 2);
    EXPECT(surface_tensor.shape(1) == 2 * 2);

    Field columns_scalar_field = nodes_fs.createField<double>(option::name("scalar"));
    Field columns_vector_field = nodes_fs.createField<double>(option::name("vector") | option::variables(2));
    Field columns_tensor_field = nodes_fs.createField<double>(option::name("tensor") | option::variables(2 * 2));

    EXPECT(columns_scalar_field.name() == std::string("scalar"));
    EXPECT(columns_vector_field.name() == std::string("vector"));
    EXPECT(columns_tensor_field.name() == std::string("tensor"));

    EXPECT(columns_scalar_field.size() == nodes_fs.nb_nodes() * nb_levels);
    EXPECT(columns_vector_field.size() == nodes_fs.nb_nodes() * nb_levels * 2);
    EXPECT(columns_tensor_field.size() == nodes_fs.nb_nodes() * nb_levels * 2 * 2);

    EXPECT(columns_scalar_field.rank() == 2);
    EXPECT(columns_vector_field.rank() == 3);
    EXPECT(columns_tensor_field.rank() == 3);

    EXPECT(columns_scalar_field.levels() == nb_levels);
    EXPECT(columns_vector_field.levels() == nb_levels);
    EXPECT(columns_tensor_field.levels() == nb_levels);

    array::ArrayView<double, 2> columns_scalar = array::make_view<double, 2>(columns_scalar_field);
    array::ArrayView<double, 3> columns_vector = array::make_view<double, 3>(columns_vector_field);
    array::ArrayView<double, 3> columns_tensor = array::make_view<double, 3>(columns_tensor_field);

    EXPECT(columns_scalar.shape(0) == nodes_fs.nb_nodes());
    EXPECT(columns_vector.shape(0) == nodes_fs.nb_nodes());
    EXPECT(columns_tensor.shape(0) == nodes_fs.nb_nodes());
    EXPECT(columns_scalar.shape(1) == nb_levels);
    EXPECT(columns_vector.shape(1) == nb_levels);
    EXPECT(columns_tensor.shape(1) == nb_levels);
    EXPECT(columns_vector.shape(2) == 2);
    EXPECT(columns_tensor.shape(2) == 2 * 2);

    Field field                  = nodes_fs.createField<int>(option::name("partition"));
    array::ArrayView<int, 2> arr = array::make_view<int, 2>(field);
    arr.assign(int(mpi::comm().rank()));
    // field->dump( Log::info() );
    nodes_fs.haloExchange(field);
    // field->dump( Log::info() );

    Field field2 = nodes_fs.createField<int>(option::name("partition2") | option::variables(2));
    Log::info() << "field2.rank() = " << field2.rank() << std::endl;
    array::ArrayView<int, 3> arr2 = array::make_view<int, 3>(field2);
    arr2.assign(int(mpi::comm().rank()));

    // field2->dump( Log::info() );
    nodes_fs.haloExchange(field2);
    // field2->dump( Log::info() );

    Log::info() << nodes_fs.checksum(field) << std::endl;

    size_t root     = mpi::comm().size() - 1;
    Field glb_field = nodes_fs.createField(option::name("partition") | option::datatype(field.datatype()) |
                                           option::levels(field.levels()) | option::variables(field.variables()) |
                                           option::global(root));
    nodes_fs.gather(field, glb_field);

    EXPECT(glb_field.rank() == field.rank());
    EXPECT(glb_field.levels() == nb_levels);
    EXPECT(field.levels() == nb_levels);

    Log::info() << "field = " << field << std::endl;
    Log::info() << "global_field = " << glb_field << std::endl;

    Log::info() << "local points = " << nodes_fs.nb_nodes() << std::endl;
    Log::info() << "grid points = " << grid.size() << std::endl;
    Log::info() << "glb_field.shape(0) = " << glb_field.shape(0) << std::endl;

    EXPECT(glb_field.metadata().get<bool>("global") == true);
    EXPECT(glb_field.metadata().get<int>("owner") == int(root));

    // glb_field->dump( Log::info() );

    if (mpi::comm().rank() == root) {
        glb_field.metadata().set("test_broadcast", 123);
    }

    arr.assign(-1);
    nodes_fs.scatter(glb_field, field);
    EXPECT(field.metadata().get<int>("test_broadcast") == 123);
    nodes_fs.haloExchange(field);
    // field->dump( Log::info() );

    Log::info() << nodes_fs.checksum(field) << std::endl;

    FieldSet fields;
    fields.add(field);
    fields.add(field2);
    Log::info() << nodes_fs.checksum(fields) << std::endl;

    Log::info() << "Testing collectives for nodes scalar field" << std::endl;
    {
        Field field                         = surface_scalar_field;
        const functionspace::NodeColumns fs = nodes_fs;

        double max;
        double min;
        double sum;
        double mean;
        double stddev;
        idx_t N;
        gidx_t gidx_max;
        gidx_t gidx_min;

        auto sfc_arr = array::make_view<double, 1>(field);
        sfc_arr.assign(mpi::comm().rank() + 1);
        fs.maximum(surface_scalar_field, max);
        EXPECT(max == double(mpi::comm().size()));

        fs.minimum(surface_scalar_field, min);
        EXPECT(min == 1);

        fs.maximumAndLocation(field, max, gidx_max);
        EXPECT(max == double(mpi::comm().size()));
        Log::info() << "global index for maximum: " << gidx_max << std::endl;

        fs.minimumAndLocation(field, min, gidx_min);
        EXPECT(min == 1);
        Log::info() << "global index for minimum: " << gidx_min << std::endl;

        fs.orderIndependentSum(field, sum, N);
        Log::info() << "oisum: " << sum << std::endl;
        Log::info() << "oiN: " << N << std::endl;

        fs.sum(field, sum, N);
        Log::info() << "sum: " << sum << std::endl;
        Log::info() << "N: " << N << std::endl;

        fs.mean(field, mean, N);
        Log::info() << "mean: " << mean << std::endl;
        Log::info() << "N: " << N << std::endl;

        fs.meanAndStandardDeviation(field, mean, stddev, N);
        Log::info() << "mean: " << mean << std::endl;
        Log::info() << "standard deviation: " << stddev << std::endl;
        Log::info() << "N: " << N << std::endl;

        int sumint;
        fs.orderIndependentSum(field, sumint, N);
        Log::info() << "sumint: " << sumint << std::endl;

        fs.sum(field, sumint, N);
        Log::info() << "sumint: " << sumint << std::endl;
    }

    Log::info() << "Testing collectives for nodes vector field" << std::endl;
    {
        Field& field                        = surface_vector_field;
        const functionspace::NodeColumns fs = nodes_fs;

        std::vector<double> max;
        std::vector<double> min;
        std::vector<double> sum;
        std::vector<double> mean;
        std::vector<double> stddev;
        idx_t N;
        std::vector<gidx_t> gidx_max;
        std::vector<gidx_t> gidx_min;

        auto vec_arr = array::make_view<double, 2>(field);
        vec_arr.assign(mpi::comm().rank() + 1);
        fs.maximum(field, max);
        std::vector<double> check_max(field.variables(), mpi::comm().size());
        EXPECT(max == check_max);

        fs.minimum(field, min);
        std::vector<double> check_min(field.variables(), 1);
        EXPECT(min == check_min);

        fs.maximumAndLocation(field, max, gidx_max);
        EXPECT(max == check_max);
        Log::info() << "global index for maximum: " << gidx_max << std::endl;

        fs.minimumAndLocation(field, min, gidx_min);
        EXPECT(min == check_min);
        Log::info() << "global index for minimum: " << gidx_min << std::endl;

        fs.orderIndependentSum(field, sum, N);
        Log::info() << "oisum: " << sum << std::endl;
        Log::info() << "oiN: " << N << std::endl;

        fs.mean(field, mean, N);
        Log::info() << "mean: " << mean << std::endl;
        Log::info() << "N: " << N << std::endl;

        fs.meanAndStandardDeviation(field, mean, stddev, N);
        Log::info() << "mean: " << mean << std::endl;
        Log::info() << "standard deviation: " << stddev << std::endl;
        Log::info() << "N: " << N << std::endl;

        std::vector<int> sumint;
        fs.orderIndependentSum(field, sumint, N);
        Log::info() << "sumint: " << sumint << std::endl;

        fs.sum(field, sumint, N);
        Log::info() << "sumint: " << sumint << std::endl;
    }

    Log::info() << "Testing collectives for columns scalar field" << std::endl;
    if (1) {
        Field& field                        = columns_scalar_field;
        const functionspace::NodeColumns fs = nodes_fs;
        double max;
        double min;
        double sum;
        double mean;
        double stddev;
        idx_t N;
        gidx_t gidx_max;
        gidx_t gidx_min;
        idx_t level;

        EXPECT(field.levels() == nb_levels);

        auto arr = array::make_view<double, 2>(field);
        arr.assign(mpi::comm().rank() + 1);
        fs.maximum(field, max);
        EXPECT(max == double(mpi::comm().size()));

        fs.minimum(field, min);
        EXPECT(min == 1);

        fs.maximumAndLocation(field, max, gidx_max, level);
        EXPECT(max == double(mpi::comm().size()));
        Log::info() << "global index for maximum: " << gidx_max << std::endl;
        Log::info() << "level for maximum: " << level << std::endl;

        fs.minimumAndLocation(field, min, gidx_min, level);
        EXPECT(min == 1);
        Log::info() << "global index for minimum: " << gidx_min << std::endl;
        Log::info() << "level for minimum: " << level << std::endl;

        fs.orderIndependentSum(field, sum, N);
        Log::info() << "order independent sum: " << sum << std::endl;
        Log::info() << "N: " << N << std::endl;

        fs.sum(field, sum, N);
        Log::info() << "sum: " << sum << std::endl;
        Log::info() << "N: " << N << std::endl;

        fs.mean(field, mean, N);
        Log::info() << "mean: " << mean << std::endl;
        Log::info() << "N: " << N << std::endl;

        fs.meanAndStandardDeviation(field, mean, stddev, N);
        Log::info() << "mean: " << mean << std::endl;
        Log::info() << "standard deviation: " << stddev << std::endl;
        Log::info() << "N: " << N << std::endl;

        int sumint;
        fs.orderIndependentSum(field, sumint, N);
        Log::info() << "order independent sum in int: " << sumint << std::endl;

        fs.sum(field, sumint, N);
        Log::info() << "sum in int: " << sumint << std::endl;

        Field max_per_level("max", array::make_datatype<double>(), array::make_shape(nb_levels));
        Field min_per_level("min", array::make_datatype<double>(), array::make_shape(nb_levels));
        Field sum_per_level("sum", array::make_datatype<double>(), array::make_shape(nb_levels));
        Field mean_per_level("mean", array::make_datatype<double>(), array::make_shape(nb_levels));
        Field stddev_per_level("stddev", array::make_datatype<double>(), array::make_shape(nb_levels));
        Field gidx_per_level("gidx", array::make_datatype<gidx_t>(), array::make_shape(nb_levels));

        fs.maximumPerLevel(field, max_per_level);
        // max_per_level.dump(Log::info());
        fs.minimumPerLevel(field, min_per_level);
        // min_per_level.dump(Log::info());
        fs.sumPerLevel(field, sum_per_level, N);
        // sum_per_level.dump(Log::info());
        fs.meanPerLevel(field, mean_per_level, N);
        // mean_per_level.dump(Log::info());
        fs.meanAndStandardDeviationPerLevel(field, mean_per_level, stddev_per_level, N);
        // mean_per_level.dump(Log::info());
        // stddev_per_level.dump(Log::info());
        fs.orderIndependentSumPerLevel(field, sum_per_level, N);
        // sum_per_level.dump(Log::info());
    }

    Log::info() << "Testing collectives for columns vector field" << std::endl;
    if (1) {
        Field& field                        = columns_vector_field;
        const functionspace::NodeColumns fs = nodes_fs;
        idx_t nvar                          = field.variables();
        std::vector<double> max;
        std::vector<double> min;
        std::vector<double> sum;
        std::vector<double> mean;
        std::vector<double> stddev;
        idx_t N;
        std::vector<gidx_t> gidx_max;
        std::vector<gidx_t> gidx_min;
        std::vector<idx_t> levels;

        auto vec_arr = array::make_view<double, 3>(field);
        vec_arr.assign(mpi::comm().rank() + 1);
        fs.maximum(field, max);
        std::vector<double> check_max(nvar, mpi::comm().size());
        EXPECT(max == check_max);

        fs.minimum(field, min);
        std::vector<double> check_min(nvar, 1);
        EXPECT(min == check_min);

        fs.maximumAndLocation(field, max, gidx_max, levels);
        EXPECT(max == check_max);
        Log::info() << "global index for maximum: " << gidx_max << std::endl;

        fs.minimumAndLocation(field, min, gidx_min, levels);
        EXPECT(min == check_min);
        Log::info() << "global index for minimum: " << gidx_min << std::endl;

        fs.orderIndependentSum(field, sum, N);
        Log::info() << "oisum: " << sum << std::endl;
        Log::info() << "N: " << N << std::endl;

        fs.mean(field, mean, N);
        Log::info() << "mean: " << mean << std::endl;
        Log::info() << "N: " << N << std::endl;

        fs.meanAndStandardDeviation(field, mean, stddev, N);
        Log::info() << "mean: " << mean << std::endl;
        Log::info() << "standard deviation: " << stddev << std::endl;
        Log::info() << "N: " << N << std::endl;

        std::vector<int> sumint;
        fs.orderIndependentSum(field, sumint, N);
        Log::info() << "sumint: " << sumint << std::endl;

        fs.sum(field, sumint, N);
        Log::info() << "sumint: " << sumint << std::endl;

        Field max_per_level("max", array::make_datatype<double>(), array::make_shape(nb_levels, nvar));
        Field min_per_level("min", array::make_datatype<double>(), array::make_shape(nb_levels, nvar));
        Field sum_per_level("sum", array::make_datatype<double>(), array::make_shape(nb_levels, nvar));
        Field mean_per_level("mean", array::make_datatype<double>(), array::make_shape(nb_levels, nvar));
        Field stddev_per_level("stddev", array::make_datatype<double>(), array::make_shape(nb_levels, nvar));
        Field gidx_per_level("gidx", array::make_datatype<gidx_t>(), array::make_shape(nb_levels, nvar));

        fs.maximumPerLevel(field, max_per_level);
        // max_per_level.dump(Log::info());

        fs.minimumPerLevel(field, min_per_level);
        // min_per_level.dump(Log::info());

        fs.sumPerLevel(field, sum_per_level, N);
        // sum_per_level.dump(Log::info());

        fs.meanPerLevel(field, mean_per_level, N);
        // mean_per_level.dump(Log::info());

        fs.meanAndStandardDeviationPerLevel(field, mean_per_level, stddev_per_level, N);
        // mean_per_level.dump(Log::info());
        // stddev_per_level.dump(Log::info());

        fs.orderIndependentSumPerLevel(field, sum_per_level, N);
        // sum_per_level.dump(Log::info());
    }

    Field tmp = nodes_fs.createField(option::datatypeT<double>() | option::global(0) | option::levels(10) |
                                     option::name("tmp"));
}

CASE("test_SpectralFunctionSpace") {
    idx_t truncation = 159;
    idx_t nb_levels  = 10;
    idx_t nspec2g    = (truncation + 1) * (truncation + 2);

    Spectral spectral_fs(truncation);

    Field surface_scalar_field = spectral_fs.createField<double>(option::name("scalar"));

    EXPECT(surface_scalar_field.name() == std::string("scalar"));

    EXPECT(surface_scalar_field.size() == nspec2g);

    EXPECT(surface_scalar_field.rank() == 1);

    auto surface_scalar = array::make_view<double, 1>(surface_scalar_field);

    EXPECT(surface_scalar.shape(0) == nspec2g);

    Field columns_scalar_field = spectral_fs.createField<double>(option::name("scalar") | option::levels(nb_levels));

    EXPECT(columns_scalar_field.name() == std::string("scalar"));

    EXPECT(columns_scalar_field.size() == nspec2g * nb_levels);

    EXPECT(columns_scalar_field.rank() == 2);

    auto columns_scalar = array::make_view<double, 2>(columns_scalar_field);

    EXPECT(columns_scalar.shape(0) == nspec2g);
    EXPECT(columns_scalar.shape(1) == nb_levels);
}

#if ATLAS_HAVE_TRANS

CASE("test_SpectralFunctionSpace_trans_dist") {
    trans::Trans trans(Grid("F80"), 159);
    idx_t nb_levels(10);

    Spectral spectral_fs(trans);
    idx_t nspec2 = spectral_fs.nb_spectral_coefficients();

    Field surface_scalar_field = spectral_fs.createField<double>(option::name("scalar"));

    EXPECT(surface_scalar_field.name() == std::string("scalar"));

    EXPECT(surface_scalar_field.size() == nspec2);

    EXPECT(surface_scalar_field.rank() == 1);

    auto surface_scalar = array::make_view<double, 1>(surface_scalar_field);

    EXPECT(surface_scalar.shape(0) == nspec2);
    // size_t surface_scalar_shape[] = { nspec2 };
    // EXPECT( eckit::testing::make_view(
    // surface_scalar.shape(),surface_scalar.shape()+1) ==
    // eckit::testing::make_view(surface_scalar_shape,surface_scalar_shape+1) );

    Field columns_scalar_field = spectral_fs.createField<double>(option::name("scalar") | option::levels(nb_levels));

    EXPECT(columns_scalar_field.name() == std::string("scalar"));

    EXPECT(columns_scalar_field.size() == nspec2 * nb_levels);

    EXPECT(columns_scalar_field.rank() == 2);

    auto columns_scalar = array::make_view<double, 2>(columns_scalar_field);

    EXPECT(columns_scalar.shape(0) == nspec2);
    EXPECT(columns_scalar.shape(1) == nb_levels);
    // size_t columns_scalar_shape[] = { nspec2, nb_levels };
    // EXPECT(eckit::testing::make_view(columns_scalar.shape(),columns_scalar.shape()+2)
    // == eckit::testing::make_view(columns_scalar_shape,columns_scalar_shape+2));
}
CASE("test_SpectralFunctionSpace_trans_global") {
    idx_t nb_levels(10);
    idx_t truncation = 159;

    Spectral spectral_fs(truncation, option::levels(nb_levels));
    idx_t nspec2g = spectral_fs.nb_spectral_coefficients_global();

    Field surface_scalar_field =
        spectral_fs.createField<double>(option::name("scalar") | option::levels(false) | option::global());

    EXPECT(surface_scalar_field.name() == std::string("scalar"));

    if (eckit::mpi::comm().rank() == 0) {
        EXPECT(surface_scalar_field.size() == nspec2g);
    }

    EXPECT(surface_scalar_field.rank() == 1);

    EXPECT(surface_scalar_field.metadata().get<bool>("global"));

    EXPECT(surface_scalar_field.metadata().get<size_t>("owner") == 0);

    auto surface_scalar = array::make_view<double, 1>(surface_scalar_field);

    if (eckit::mpi::comm().rank() == 0) {
        EXPECT(surface_scalar.shape(0) == nspec2g);
    }
    Field columns_scalar_field = spectral_fs.createField<double>(option::name("scalar") | option::global());

    EXPECT(columns_scalar_field.name() == std::string("scalar"));

    if (eckit::mpi::comm().rank() == 0) {
        EXPECT(columns_scalar_field.size() == nspec2g * nb_levels);
    }
    else {
        EXPECT(columns_scalar_field.size() == 0);
    }

    EXPECT(columns_scalar_field.rank() == 2);

    auto columns_scalar = array::make_view<double, 2>(columns_scalar_field);

    if (eckit::mpi::comm().rank() == 0) {
        EXPECT(columns_scalar.shape(0) == nspec2g);
        EXPECT(columns_scalar.shape(1) == nb_levels);
    }
}
CASE("test_SpectralFunctionSpace_norm") {
    trans::Trans trans(Grid("F80"), 159);
    size_t nb_levels(10);

    Spectral spectral_fs(trans);

    Field twoD_field   = spectral_fs.createField<double>(option::name("2d"));
    Field threeD_field = spectral_fs.createField<double>(option::name("3d") | option::levels(nb_levels));

    // Set first wave number
    {
        auto twoD = array::make_view<double, 1>(twoD_field);
        twoD.assign(0.);
        if (mpi::comm().rank() == 0) {
            twoD(0) = 1.;
        }

        auto threeD = array::make_view<double, 2>(threeD_field);
        threeD.assign(0.);
        for (size_t jlev = 0; jlev < nb_levels; ++jlev) {
            if (mpi::comm().rank() == 0) {
                threeD((size_t)0, jlev) = jlev;
            }
        }
    }

    double twoD_norm(0.);
    std::vector<double> threeD_norms(threeD_field.levels(), 0.);

    spectral_fs.norm(twoD_field, twoD_norm);
    if (not threeD_field.contiguous()) {
        EXPECT_THROWS_AS(spectral_fs.norm(threeD_field, threeD_norms), eckit::Exception);
        return;
    }
    else {
        EXPECT_NO_THROW(spectral_fs.norm(threeD_field, threeD_norms));
    }

    if (eckit::mpi::comm().rank() == 0) {
        EXPECT(eckit::types::is_approximately_equal(twoD_norm,
                                                    1.0));  // before is_approximately_equal tolerance was 1.e-10
        for (size_t jlev = 0; jlev < nb_levels; ++jlev) {
            EXPECT(eckit::types::is_approximately_equal(
                threeD_norms[jlev],
                double(jlev)));  // before is_approximately_equal tolerance was 1.e-10
        }
    }
}
CASE("test_functionspace_grid") {
    // Create list of points and construct Grid from them
    std::vector<PointXY> points = {
        {0.0, 5.0}, {0.0, 0.0}, {10.0, 0.0}, {15.0, 0.0}, {5.0, 5.0}, {15.0, 5.0}
    };
    UnstructuredGrid grid(points);

    // Create Mesh from Grid
    Mesh mesh_from_grid(grid);

    // Create Mesh from list of points using MeshBuilder
    std::vector<double> lons;
    std::vector<double> lats;
    for (const auto& point : points) {
        lons.push_back(point.x());
        lats.push_back(point.y());
    }
    std::vector<int> ghosts(6, 0);
    std::vector<gidx_t> global_indices(6);
    std::iota(global_indices.begin(), global_indices.end(), 1);
    const idx_t remote_index_base = 0;
    std::vector<idx_t> remote_indices(6);
    std::iota(remote_indices.begin(), remote_indices.end(), remote_index_base);
    std::vector<int> partitions(6, 0);
    std::vector<std::array<gidx_t, 3>> tri_boundary_nodes = {{{3, 6, 5}}, {{3, 4, 6}}};
    std::vector<gidx_t> tri_global_indices                = {1, 2};
    std::vector<std::array<gidx_t, 4>> quad_boundary_nodes = {{{1, 2, 3, 5}}};
    std::vector<gidx_t> quad_global_indices                = {3};
    const mesh::MeshBuilder mesh_builder{};
    const Mesh mesh_from_meshbuilder = mesh_builder(
        lons, lats, ghosts, global_indices, remote_indices, remote_index_base, partitions,
        tri_boundary_nodes, tri_global_indices, quad_boundary_nodes, quad_global_indices);

    // Create Cell/Edge/NodeColumns and PointCloud FunctionSpaces that will save the Grid on construction
    functionspace::CellColumns cells_from_grid(mesh_from_grid);
    functionspace::EdgeColumns edges_from_grid(mesh_from_grid);
    functionspace::NodeColumns nodes_from_grid(mesh_from_grid);
    functionspace::PointCloud pointcloud_from_grid(grid);

    // Create Cell/Edge/NodeColumns and PointCloud FunctionSpaces that will generate the Grid when requested
    functionspace::CellColumns cells_from_meshbuilder(mesh_from_meshbuilder);
    functionspace::EdgeColumns edges_from_meshbuilder(mesh_from_meshbuilder);
    functionspace::NodeColumns nodes_from_meshbuilder(mesh_from_meshbuilder);
    functionspace::PointCloud pointcloud_from_points(points);

    // All Grids should match original
    EXPECT(cells_from_grid.grid().uid() == grid.uid());
    EXPECT(edges_from_grid.grid().uid() == grid.uid());
    EXPECT(nodes_from_grid.grid().uid() == grid.uid());
    EXPECT(pointcloud_from_grid.grid().uid() == grid.uid());
    EXPECT(cells_from_meshbuilder.grid().uid() == grid.uid());
    EXPECT(edges_from_meshbuilder.grid().uid() == grid.uid());
    EXPECT(nodes_from_meshbuilder.grid().uid() == grid.uid());
    EXPECT(pointcloud_from_points.grid().uid() == grid.uid());

    // Repeat with a StructuredGrid to check StructuredColumns and BlockStructuredColumns
    Grid structured_grid("O8");
    Mesh structured_mesh = StructuredMeshGenerator().generate(structured_grid);
    functionspace::StructuredColumns sc(structured_grid);
    functionspace::BlockStructuredColumns bsc(structured_grid);
    EXPECT(sc.grid() == structured_grid);
    EXPECT(bsc.grid() == structured_grid);

    // Check that Spectral throws exception
    functionspace::Spectral spectral(159);
    EXPECT_THROWS(spectral.grid());
}

#endif

//-----------------------------------------------------------------------------

}  // namespace test
}  // namespace atlas

int main(int argc, char** argv) {
    return atlas::test::run(argc, argv);
}
