/*
 * (C) Crown Copyright 2026 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */
#include <functional>
#include <string>

#include "atlas/array.h"
#include "atlas/field/Field.h"
#include "atlas/functionspace/FunctionSpace.h"
#include "atlas/functionspace/Spectral.h"
#include "atlas/functionspace/StructuredColumns.h"
#include "atlas/grid/Grid.h"
#include "atlas/grid/StructuredGrid.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/runtime/Exception.h"
#include "atlas/runtime/Log.h"
#include "atlas/trans/Trans.h"
#include "atlas/util/CoordinateEnums.h"
#include "atlas/util/function/VortexRollup.h"

#include "tests/AtlasTestEnvironment.h"

namespace atlas {
namespace test {

// Test Fixture ----------------------------------------------------------------

struct adjointTestFixture {
    Grid gaussGrid_;
    functionspace::StructuredColumns gaussFunctionSpace_;
    Field gaussField_;
};

// This function sets up the function spaces and vortex rollup gauss field given the type (O or F) and number (int)
void prepareTestFixture(adjointTestFixture& testFixture, std::string& grid_type, int& grid_number) {
    std::string grid_name = grid_type + std::to_string(grid_number);
    testFixture.gaussGrid_ = StructuredGrid(grid_name);
    testFixture.gaussFunctionSpace_ = functionspace::StructuredColumns(testFixture.gaussGrid_);

    testFixture.gaussField_ = testFixture.gaussFunctionSpace_.createField<double>(option::name("gauss_field"));
    {
        const array::ArrayView<double, 2> lonlat = array::make_view<double, 2>(testFixture.gaussFunctionSpace_.lonlat());
        array::ArrayView<double, 1> view         = array::make_view<double, 1>(testFixture.gaussField_);
        for (idx_t i = 0; i < testFixture.gaussFunctionSpace_.size(); ++i) {
            view(i) = util::function::vortex_rollup(lonlat(i, LON), lonlat(i, LAT), 1.);
        }
    }
}


/// Helper functions -----------------------------------------------------------

// Dot product of atlas fields
// Note this works for any two Gauss grid-point fields as they are real
// but only works for taking the dot of a spectral field with itself as they are complex
double dotProd(const Field& a, const Field& b) {
    const auto ghost = [&] {
        ATLAS_ASSERT(a.functionspace().type() == b.functionspace().type());
        if (a.functionspace().type() == "Spectral") {
            return std::function<int(idx_t)>([](idx_t ){return 0;});
        }
        return std::function<int(idx_t)>(array::make_view<int, 1>(a.functionspace().ghost()));
    }();
    
    double prod{};

    const auto aView = array::make_view<double, 1>(a);
    const auto bView = array::make_view<double, 1>(b);

    for (size_t i = 0; i < a.size(); ++i) {
        if (ghost(i)) {
            continue;
        }
        prod += aView(i) * bView(i);
    }
    mpi::comm().allReduceInPlace(prod, eckit::mpi::Operation::SUM);
    return prod;
}

/// Main test function ------------------------------------------------------
// Takes the Gauss grid type (O or F) and the grid number (int)
// We compute all inner products in the real (Gauss grid) space rather than the complex (Spectral)
// Adjoint tests are <Tx,y> = <x, T*y> with Tx = y or T*y = x
void testFunction(std::string grid_type, int grid_number) {

    // Construct initial gauss field and spectral function space.
    adjointTestFixture testFixture;
    prepareTestFixture(testFixture, grid_type, grid_number);

    trans::Trans trans_(testFixture.gaussGrid_, grid_number - 1);
    functionspace::Spectral spectralFunctionspace(trans_);

    /// Test dirtrans ----------------------------
    // Construct spectral field by applying TL to Gauss. Compute dot product.
    Field spectralField = spectralFunctionspace.createField<double>(option::name("spectral_field"));
    trans_.dirtrans(testFixture.gaussField_, spectralField);
    double yDotY = dotProd(spectralField, spectralField);

    // Construct adjoint spectral field. Compute dot product (dirtrans_adj)
    Field adjointSpectralField =
        testFixture.gaussFunctionSpace_.createField<double>(option::name("adjoint_spectral_field"));
    trans_.dirtrans_adj(spectralField, adjointSpectralField);
    double xDotAdjY = dotProd(testFixture.gaussField_, adjointSpectralField);

    // Adjoint test <y,y> = <x, T*y>
    Log::error() << "dirtrans test" << std::endl;
    Log::error() << "<y,y>: " << yDotY << " and <x,T*y>: " << xDotAdjY << std::endl;
    EXPECT_APPROX_EQ(yDotY / xDotAdjY, 1., 1e-12);

    //// Test invtrans ---------------------------
    // Construct Gauss from Spectral field. Compute dot product (invtrans)
    Field secondGaussField = testFixture.gaussFunctionSpace_.createField<double>(option::name("second_gauss_field"));
    trans_.invtrans(spectralField, secondGaussField);
    double xDotX = dotProd(secondGaussField,secondGaussField);

    // Construct adjoint Gauss field. Compute dot product (invtrans_adj)
    Field adjointGaussField = spectralFunctionspace.createField<double>(option::name("adjoint_gauss_field"));
    trans_.invtrans_adj(secondGaussField, adjointGaussField);
    double TxDotY = dotProd(adjointGaussField, spectralField);

    // Adjoint test
    Log::error() << "invtrans test" << std::endl;
    Log::error() << "<x,x>: " << xDotX << " and <Tx,y>: " << TxDotY << std::endl;
    EXPECT_APPROX_EQ(xDotX / TxDotY, 1., 1e-12);
}

/// Test cases.
CASE("O12") {
    testFunction("O", 12);
}

CASE("F12") {
    testFunction("F", 12);
}

CASE("O15") {
    testFunction("O", 15);
}

CASE("F15") {
    testFunction("F", 15);
}

CASE("O8") {
    testFunction("O", 8);
}

CASE("F8") {
    testFunction("F", 8);
}


//-----------------------------------------------------------------------------

}  // namespace test
}  // namespace atlas

int main(int argc, char** argv) {
    return atlas::test::run(argc, argv);
}