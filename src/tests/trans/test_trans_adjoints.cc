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

/// Helper functions -----------------------------------------------------------

// Dot product of atlas fields
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

Field createVortexRollup(const FunctionSpace& gaussFunctionSpace) {
    Field gaussField = gaussFunctionSpace.createField<double>(option::name("x"));
    {
        const array::ArrayView<double, 2> lonlat =
            array::make_view<double, 2>(gaussFunctionSpace.lonlat());
        array::ArrayView<double, 1> view = array::make_view<double, 1>(gaussField);
        for (idx_t i = 0; i < gaussFunctionSpace.size(); ++i) {
            view(i) = util::function::vortex_rollup(lonlat(i, LON), lonlat(i, LAT), 1.);
        }
    }
    return gaussField;
}

/// Main test function ------------------------------------------------------
// Adjoint tests are <Tx,y> = <x, T*y> with Tx = y or T*y = x
void testFunction(const GaussianGrid& gaussGrid) {
    // Construct initial gauss field and spectral function space.
    functionspace::StructuredColumns gaussFunctionSpace = functionspace::StructuredColumns(gaussGrid);
    Field gaussField                                    = createVortexRollup(gaussFunctionSpace);
    functionspace::Spectral spectralFunctionSpace = functionspace::Spectral(gaussGrid.N()-1);
    trans::Trans trans_(gaussFunctionSpace, spectralFunctionSpace);

    /// Test dirtrans ----------------------------
    // Construct spectral field by applying TL to Gauss. Compute dot product.
    Field spectralField = spectralFunctionSpace.createField<double>(option::name("y"));
    trans_.dirtrans(gaussField, spectralField);
    double yDotY = dotProd(spectralField, spectralField);

    // Construct adjoint spectral field. Compute dot product (dirtrans_adj)
    Field adjointSpectralField = gaussFunctionSpace.createField<double>(option::name("T*y"));
    trans_.dirtrans_adj(spectralField, adjointSpectralField);
    double xDotAdjY = dotProd(gaussField, adjointSpectralField);

    // Adjoint test <y,y> = <x, T*y>
    Log::error() << "dirtrans test" << std::endl;
    Log::error() << "<y,y>: " << yDotY << " and <x,T*y>: " << xDotAdjY << std::endl;
    EXPECT_APPROX_EQ(yDotY / xDotAdjY, 1., 1e-12);

    //// Test invtrans ---------------------------
    // Construct Gauss from Spectral field. Compute dot product (invtrans)
    Field secondGaussField = gaussFunctionSpace.createField<double>(option::name("x"));
    trans_.invtrans(spectralField, secondGaussField);
    double xDotX = dotProd(secondGaussField,secondGaussField);

    // Construct adjoint Gauss field. Compute dot product (invtrans_adj)
    Field adjointGaussField = spectralFunctionSpace.createField<double>(option::name("T*x"));
    trans_.invtrans_adj(secondGaussField, adjointGaussField);
    double AdjXDotY = dotProd(adjointGaussField, spectralField);

    // Adjoint test
    Log::error() << "invtrans test" << std::endl;
    Log::error() << "<x,x>: " << xDotX << " and <Tx,y>: " << AdjXDotY << std::endl;
    EXPECT_APPROX_EQ(xDotX / AdjXDotY, 1., 1e-12);
}

/// Test cases.
CASE("O12") {
    testFunction(Grid("O12"));
}

CASE("F12") {
    testFunction(Grid("F12"));
}

CASE("O15") {
    testFunction(Grid("O15"));
}

CASE("F15") {
    testFunction(Grid("F15"));
}

CASE("O8") {
    testFunction(Grid("O8"));
}

CASE("F8") {
    testFunction(Grid("F8"));
}


//-----------------------------------------------------------------------------

}  // namespace test
}  // namespace atlas

int main(int argc, char** argv) {
    return atlas::test::run(argc, argv);
}