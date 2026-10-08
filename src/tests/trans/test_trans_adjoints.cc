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

/// Dot Products -----------------------------------------------------------
// Dot product of 1D grid point fields
double dotProd1DGridPoint(const Field& a, const Field& b) {
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

// Dot product of 2D grid point fields
double dotProd2DGridPoint(const Field& a, const Field& b) {
    auto ghost = array::make_view<const int, 1>(a.functionspace().ghost());

    double prod{};
    const auto aView = array::make_view<double, 2>(a);
    const auto bView = array::make_view<double, 2>(b);

    for (idx_t i = 0; i < a.shape(0); ++i) {
        if (ghost(i)) {
            continue;
        }
        for (idx_t j = 0; j < a.shape(1); ++j) {
            prod += aView(i, j) * bView(i, j);
        }
    }
    mpi::comm().allReduceInPlace(prod, eckit::mpi::Operation::SUM);
    return prod;
}

// Dot product of spectral fields
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

// Dot product of grid point fields
double dotProdGridPoint(const Field& a, const Field& b) {
    if (a.rank() == 1) {
        return dotProd1DGridPoint(a, b);
    }
    else if (a.rank() == 2) {
        return dotProd2DGridPoint(a, b);
    }
    ATLAS_NOTIMPLEMENTED;
}

// Dot product of atlas fields
double dotProd(const Field& a, const Field& b) {
    ATLAS_ASSERT(a.functionspace().type() == b.functionspace().type());
    if (a.functionspace().type() == "Spectral") {
        return dotProdSpectral(a, b);
    }
    return dotProdGridPoint(a, b);
}

/// Vortex Rollups -------------------------------------------
// Vortex rollup for 1D scalar field
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

// Vortex rollup for 2D wind field
Field createVortexRollupWinds(const FunctionSpace& gaussFunctionSpace) {
    Field windField = gaussFunctionSpace.createField<double>(option::name("x_winds") | option::variables(2));
    {
        const array::ArrayView<double, 2> lonlat = array::make_view<double, 2>(gaussFunctionSpace.lonlat());
        array::ArrayView<double, 2> view = array::make_view<double, 2>(windField);
        
        for (idx_t i = 0; i < gaussFunctionSpace.size(); ++i) {
            view(i, 0) = util::function::vortex_rollup(lonlat(i, LON), lonlat(i, LAT), 1.);
            view(i, 1) = util::function::vortex_rollup(lonlat(i, LON), lonlat(i, LAT), 1.);
        }
    }
    return windField;
}

/// Test Fixture --------------------------------------------

struct testFixture {
    functionspace::StructuredColumns gaussFunctionSpace;
    functionspace::Spectral spectralFunctionSpace;
    trans::Trans trans_;
};

testFixture createTestFixture(const GaussianGrid& gaussGrid) {
    testFixture fixture;
    fixture.gaussFunctionSpace = functionspace::StructuredColumns(gaussGrid);
    fixture.spectralFunctionSpace = functionspace::Spectral(gaussGrid.N()-1);
    fixture.trans_ = trans::Trans(fixture.gaussFunctionSpace, fixture.spectralFunctionSpace);
    return fixture;
}

/// Main test functions ------------------------------------------------------
// Test dirtrans and invtrans adjoints
void testDirtrans(const testFixture& testFixture) {
    // x: Gauss field
    Field gaussField = createVortexRollup(testFixture.gaussFunctionSpace);

    /// Test dirtrans (T) ----------------------------
    // y = Tx: Spectral field constructed from Gauss field and <y,y> computed
    Field spectralField = testFixture.spectralFunctionSpace.createField<double>(option::name("y"));
    testFixture.trans_.dirtrans(gaussField, spectralField);
    const double yDotY = dotProd(spectralField, spectralField);

    // T*Tx: Adjoint spectral field constructed and <x, T*y> computed
    Field adjointSpectralField = testFixture.gaussFunctionSpace.createField<double>(option::name("T*y"));
    testFixture.trans_.dirtrans_adj(spectralField, adjointSpectralField);
    const double xDotAdjY = dotProd(gaussField, adjointSpectralField);

    // Adjoint test <y,y> = <x, T*y>
    Log::error() << "dirtrans test: <y,y>=" << yDotY << ", <x,T*y>=" << xDotAdjY << std::endl;
    EXPECT_APPROX_EQ(yDotY / xDotAdjY, 1., 1e-12);

    //// Test invtrans (T) ---------------------------
    // x = Ty: Gauss field constructed from spectral field and <x,x> computed
    Field secondGaussField = testFixture.gaussFunctionSpace.createField<double>(option::name("x"));
    testFixture.trans_.invtrans(spectralField, secondGaussField);
    const double xDotX = dotProd(secondGaussField,secondGaussField);

    // T*x: Adjoint Gauss field constructed and <T*x,y> computed.
    Field adjointGaussField = testFixture.spectralFunctionSpace.createField<double>(option::name("T*x"));
    testFixture.trans_.invtrans_adj(secondGaussField, adjointGaussField);
    const double AdjXDotY = dotProd(adjointGaussField, spectralField);

    // Adjoint test <x,x> = <T*x,y>
    Log::error() << "invtrans test: <x,x>=" << xDotX << ", <Tx,y>=" << AdjXDotY << std::endl;
    EXPECT_APPROX_EQ(xDotX / AdjXDotY, 1., 1e-12);
}

// Test dirtrans_wind2vordiv and invtrans_wind2vordiv adjoints
void testWindVorDiv(const testFixture& testFixture) {
    // x: Gauss wind field (u,v)
    Field gaussWindField = createVortexRollupWinds(testFixture.gaussFunctionSpace);

    /// Test dirtrans_wind2vordiv (T) ----------------------------
    // y = Tx (u,v): Spectral fields y = (vor,div) constructed from Gauss wind and <y,y> = <vor,vor> + <div,div> computed
    Field vor = testFixture.spectralFunctionSpace.createField<double>(option::name("vor"));
    Field div = testFixture.spectralFunctionSpace.createField<double>(option::name("div"));
    testFixture.trans_.dirtrans_wind2vordiv(gaussWindField, vor, div);
    const double yDotY = dotProd(vor, vor) + dotProd(div, div);

    // T*y: Adjoint spectral field constructed and <x,T*y> computed.
    Field adjointSpectralField =
        testFixture.gaussFunctionSpace.createField<double>(option::name("T*y") | option::variables(2));
    testFixture.trans_.dirtrans_wind2vordiv_adj(vor, div, adjointSpectralField);
    const double xDotAdjY = dotProd(gaussWindField, adjointSpectralField);

    // Adjoint test <y,y> = <x,T*y>
    Log::error() << "dirtrans_wind2vordiv test: <y,y>=" << yDotY << ", <x,T*y>=" << xDotAdjY << std::endl;
    EXPECT_APPROX_EQ(yDotY / xDotAdjY, 1., 1e-12);

    /// Test invtrans_wind2vordiv (T) ----------------------------
    // x = Ty: Gauss wind field (u,v) constructed from spectral field and <x,x> computed
    Field secondGaussWindField = testFixture.gaussFunctionSpace.createField<double>(option::name("x_winds") | option::variables(2));
    testFixture.trans_.invtrans_vordiv2wind(vor, div, secondGaussWindField);
    const double xDotX = dotProd(secondGaussWindField, secondGaussWindField);

    // T*x: Adjoint Gauss fields (div,vor) constructed from second Gauss wind field and <T*x,y> computed.
    Field adjointGaussFieldVor = testFixture.spectralFunctionSpace.createField<double>(option::name("T*x_vor"));
    Field adjointGaussFieldDiv = testFixture.spectralFunctionSpace.createField<double>(option::name("T*x_div"));
    testFixture.trans_.invtrans_vordiv2wind_adj(secondGaussWindField, adjointGaussFieldVor, adjointGaussFieldDiv);
    const double AdjXDotY = dotProd(adjointGaussFieldVor, vor) + dotProd(adjointGaussFieldDiv, div);

    // Adjoint test <x,x> = <T*x,y>
    Log::error() << "invtrans_wind2vordiv test: <x,x>=" << xDotX << ", <Tx,y>=" << AdjXDotY << std::endl;
    EXPECT_APPROX_EQ(xDotX / AdjXDotY, 1., 1e-12);
}

// Test invtrans_grad adjoint
void testInvtransGrad(const testFixture& testFixture) {
    // y: Spectral field from vortex rollup of Gauss field
    Field gaussField    = createVortexRollup(testFixture.gaussFunctionSpace);
    Field spectralField = testFixture.spectralFunctionSpace.createField<double>(option::name("y"));
    testFixture.trans_.dirtrans(gaussField, spectralField);

    // x = Ty: Gauss gradient field (d/dlambda, d/dphi) constructed from spectral field and <x,x> computed.
    Field gradField = testFixture.gaussFunctionSpace.createField<double>(option::name("grad") | option::variables(2));
    array::make_view<double, 2>(gradField).assign(0.);
    testFixture.trans_.invtrans_grad(spectralField, gradField);
    const double xDotX = dotProd(gradField, gradField);

    // T*x: Adjoint spectral field constructed from Gauss gradient field and <T*x,y> computed.
    Field adjointSpectralField = testFixture.spectralFunctionSpace.createField<double>(option::name("T*x"));
    array::make_view<double, 1>(adjointSpectralField).assign(0.);
    testFixture.trans_.invtrans_grad_adj(gradField, adjointSpectralField);
    const double adjXDotY = dotProd(adjointSpectralField, spectralField);

    // Adjoint test <x,x> = <T*x,y>
    Log::info() << "invtrans_grad test: <x,x>=" << xDotX << ", <T*x,y>=" << adjXDotY << std::endl;
    EXPECT_APPROX_EQ(xDotX / adjXDotY, 1., 1e-12);
}

// Create test fixture and run all adjoint tests for given grid
void testFunction(const GaussianGrid& gaussGrid) {
    testFixture fixture = createTestFixture(Grid(gaussGrid));
    testDirtrans(fixture);
    testWindVorDiv(fixture);
    testInvtransGrad(fixture);
}


/// Test cases ---------------------------------------------------------------
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

CASE("O12") {
    testFunction(Grid("O12"));
}

CASE("F12") {
    testFunction(Grid("F12"));
}


//-----------------------------------------------------------------------------

}  // namespace test
}  // namespace atlas

int main(int argc, char** argv) {
    return atlas::test::run(argc, argv);
}
