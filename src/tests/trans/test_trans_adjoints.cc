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


double dotProdGridPoint(const Field& a, const Field& b) {
    auto ghost = array::make_view<const int, 1>(a.functionspace().ghost());

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

double dotProdSpectral(const Field& a, const Field& b) {
   /* Because the input fields are real-valued in the spatial domain,
    * their spectral coefficients exhibit conjugate symmetry: F(n, -m) = (-1)^m * F(n, m)^*.
    *
    * To optimize storage, the layout omits negative azimuthal orders (m < 0).
    * To reconstruct the full global energy integral from this half-spectrum:
    *  - Zonal modes (m = 0) are counted once (weight = 1.0).
    *  - Non-zonal modes (m > 0) are doubled (weight = 2.0) to account for the missing
    *    negative m counterparts.
    */
    double prod{};
    const auto aView = array::make_view<double, 1>(a);
    const auto bView = array::make_view<double, 1>(b);

    functionspace::Spectral spectral = functionspace::Spectral{a.functionspace()};
    ATLAS_ASSERT(spectral);

    const auto zonal_wavenumbers = spectral.zonal_wavenumbers();
    const int truncation = spectral.truncation();
    idx_t index = 0;
    for( idx_t jm=0; jm<zonal_wavenumbers.size(); ++jm ) {
        const int m = zonal_wavenumbers(jm);
        for( int n=m; n<=truncation; ++n ) {
            for (int c=0; c<2; ++c) { // 2 complex components (real, imag)
                if (m == 0) {
                    prod += aView(index) * bView(index);
                }
                else {
                    prod += 2.0 * aView(index) * bView(index);
                }
                ++index;
            }
        }
    }
    ATLAS_ASSERT(index == a.size());
    mpi::comm().allReduceInPlace(prod, eckit::mpi::Operation::SUM);
    return prod;
}

// Dot product of atlas fields
double dotProd(const Field& a, const Field& b) {
    ATLAS_ASSERT(a.functionspace().type() == b.functionspace().type());
    if (a.functionspace().type() == "Spectral") {
        return dotProdSpectral(a, b);
    }
    return dotProdGridPoint(a, b);
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
    Field gaussField = createVortexRollup(gaussFunctionSpace);
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