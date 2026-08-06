/*
 * (C) Copyright 2026 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <string>

#include "atlas/array/ArrayView.h"
#include "atlas/array/MakeView.h"
#include "atlas/field/Field.h"
#include "atlas/functionspace/Spectral.h"
#include "atlas/option.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/parallel/omp/omp.h"

#include "tests/AtlasTestEnvironment.h"

namespace atlas {
namespace test {

namespace {

constexpr int truncation = 79;
constexpr int levels     = 3;

enum class ComplexComponent { Real, Imag };

double spectral_value(int n, int m, int level, ComplexComponent component) {
    return component == ComplexComponent::Imag ? -1000. - 100. * level - 10. * m - n : 1000. + 100. * level + 10. * m + n;
}

Field create_spectral_field(const functionspace::Spectral& spectral, const std::string& name) {
    Field field = spectral.createField<double>(option::name(name) | option::levels(levels));
    array::make_view<double, 2>(field).assign(-1.);
    return field;
}

Field create_global_spectral_field(const functionspace::Spectral& spectral, const std::string& name) {
    Field field = spectral.createField<double>(option::name(name) | option::levels(levels) | option::global());
    array::make_view<double, 2>(field).assign(-1.);
    return field;
}

void set_coefficient(array::ArrayView<double, 2>& data, idx_t real, idx_t imag, int n, int m) {
    for (idx_t level = 0; level < levels; ++level) {
        data(real, level) = spectral_value(n, m, static_cast<int>(level), ComplexComponent::Real);
        data(imag, level) = spectral_value(n, m, static_cast<int>(level), ComplexComponent::Imag);
    }
}

void fill_sequential_storage_order(const functionspace::Spectral& spectral, Field& field) {
    auto data                     = array::make_view<double, 2>(field);
    const auto zonal_wavenumbers = spectral.zonal_wavenumbers();

    idx_t index = 0;
    for (idx_t jm = 0; jm < zonal_wavenumbers.size(); ++jm) {
        const int m = zonal_wavenumbers(jm);
        for (int n = m; n <= spectral.truncation(); ++n) {
            set_coefficient(data, index, index + 1, n, m);
            index += 2;
        }
    }
    EXPECT_EQ(index, field.shape(0));
}

void fill_global_sequential_storage_order(const functionspace::Spectral& spectral, Field& field) {
    auto data = array::make_view<double, 2>(field);

    idx_t index = 0;
    if (mpi::comm().rank() == 0) {
        for (int m = 0; m <= spectral.truncation(); ++m) {
            for (int n = m; n <= spectral.truncation(); ++n) {
                set_coefficient(data, index, index + 1, n, m);
                index += 2;
            }
        }
        EXPECT_EQ(index, spectral.nb_spectral_coefficients_global());
    }
    EXPECT_EQ(index, field.shape(0));
}

void fill_nasm0_storage_order(const functionspace::Spectral& spectral, Field& field) {
    auto data                     = array::make_view<double, 2>(field);
    const auto zonal_wavenumbers = spectral.zonal_wavenumbers();
    const auto nasm0             = spectral.offsets_by_zonal_wavenumber();

    for (idx_t jm = 0; jm < zonal_wavenumbers.size(); ++jm) {
        const int m = zonal_wavenumbers(jm);
        idx_t index = nasm0[m];
        for (int n = m; n <= spectral.truncation(); ++n, index+=2) {
            set_coefficient(data, index, index + 1, n, m);
        }
    }
}

void fill_inverted_nasm0_storage_order(const functionspace::Spectral& spectral, Field& field) {
    auto data                     = array::make_view<double, 2>(field);
    const auto zonal_wavenumbers = spectral.zonal_wavenumbers();
    const auto nasm0             = spectral.offsets_by_zonal_wavenumber();

    for (int n = 0; n <= spectral.truncation(); ++n) {
        for (idx_t jm = 0; jm < zonal_wavenumbers.size(); ++jm) {
            const int m = zonal_wavenumbers(jm);
            if (m > n) {
                continue;
            }
            const idx_t offset_for_zonal_wavenumber = nasm0[m];
            const idx_t index = offset_for_zonal_wavenumber + 2 * (n - m);
            set_coefficient(data, index, index + 1, n, m);
        }
    }
}

void fill_parallel_for_order(const functionspace::Spectral& spectral, Field& field) {
    auto data = array::make_view<double, 2>(field);

    spectral.parallel_for([&](idx_t real, idx_t imag, int n, int m) { set_coefficient(data, real, imag, n, m); });
}

void fill_global_parallel_for_order(const functionspace::Spectral& spectral, Field& field) {
    auto data = array::make_view<double, 2>(field);

    spectral.parallel_for(option::global(), [&](idx_t real, idx_t imag, int n, int m) { set_coefficient(data, real, imag, n, m); });
}

void expect_field_matches_function(const functionspace::Spectral& spectral, const Field& field) {
    const auto data = array::make_view<double, 2>(field);

    spectral.parallel_for([&](idx_t real, idx_t imag, int n, int m) {
        for (idx_t level = 0; level < levels; ++level) {
            EXPECT_EQ(data(real, level), spectral_value(n, m, static_cast<int>(level), ComplexComponent::Real));
            EXPECT_EQ(data(imag, level), spectral_value(n, m, static_cast<int>(level), ComplexComponent::Imag));
        }
    });
}

void expect_equivalent_fields(const Field& left, const Field& right) {
    EXPECT_EQ(left.shape(0), right.shape(0));
    EXPECT_EQ(left.shape(1), right.shape(1));

    const auto left_data  = array::make_view<double, 2>(left);
    const auto right_data = array::make_view<double, 2>(right);
    for (idx_t coefficient = 0; coefficient < left.shape(0); ++coefficient) {
        for (idx_t level = 0; level < left.shape(1); ++level) {
            EXPECT_EQ(left_data(coefficient, level), right_data(coefficient, level));
        }
    }
}

}  // namespace

//-----------------------------------------------------------------------------

CASE("test_parallel_for_count") {
    const idx_t expected_global_coefficients = (truncation + 1) * (truncation + 2);

    functionspace::Spectral spectral(truncation);
    EXPECT_EQ(spectral.nb_spectral_coefficients_global(), expected_global_coefficients);

    // Following must be updated atomically in parallel_for
    idx_t parallel_for_count = 0;
    idx_t m0_count           = 0;
    bool  imag_follows_real  = true;
    spectral.parallel_for([&](idx_t real, idx_t imag, int n, int m) {
        if (imag != real + 1) {
            atlas_omp_atomic_write
            imag_follows_real = false;
        }
        atlas_omp_atomic_update
        ++parallel_for_count;
        if (m == 0) {
            atlas_omp_atomic_update
            ++m0_count;
        }
    });
    EXPECT(imag_follows_real);
    mpi::comm().allReduceInPlace(parallel_for_count, eckit::mpi::sum());
    mpi::comm().allReduceInPlace(m0_count, eckit::mpi::sum());
    EXPECT_EQ(2 * parallel_for_count, expected_global_coefficients);
    EXPECT_EQ(m0_count, truncation + 1);
}

CASE("test_spectral_sequential_loop_order") {
    functionspace::Spectral spectral(truncation);
    Field field = create_spectral_field(spectral, "sequential");

    fill_sequential_storage_order(spectral, field);

    expect_field_matches_function(spectral, field);
}

CASE("test_spectral_nasm0_loop_order") {
    functionspace::Spectral spectral(truncation);
    Field sequential = create_spectral_field(spectral, "sequential");
    Field with_nasm0 = create_spectral_field(spectral, "with_nasm0");

    fill_sequential_storage_order(spectral, sequential);
    fill_nasm0_storage_order(spectral, with_nasm0);

    expect_equivalent_fields(sequential, with_nasm0);
}

CASE("test_spectral_inverted_loop_order") {
    functionspace::Spectral spectral(truncation);
    Field sequential = create_spectral_field(spectral, "sequential");
    Field inverted   = create_spectral_field(spectral, "inverted");

    fill_sequential_storage_order(spectral, sequential);
    fill_inverted_nasm0_storage_order(spectral, inverted);

    expect_equivalent_fields(sequential, inverted);
}

CASE("test_spectral_parallel_for_loop_order") {
    functionspace::Spectral spectral(truncation);
    Field sequential   = create_spectral_field(spectral, "sequential");
    Field parallel_for = create_spectral_field(spectral, "parallel_for");

    fill_sequential_storage_order(spectral, sequential);
    fill_parallel_for_order(spectral, parallel_for);

    expect_equivalent_fields(sequential, parallel_for);
}

CASE("test_spectral_global_parallel_for_loop_order") {
    functionspace::Spectral spectral(truncation);
    Field sequential   = create_global_spectral_field(spectral, "global_sequential");
    Field parallel_for = create_global_spectral_field(spectral, "global_parallel_for");

    fill_global_sequential_storage_order(spectral, sequential);
    fill_global_parallel_for_order(spectral, parallel_for);

    expect_equivalent_fields(sequential, parallel_for);
}

//-----------------------------------------------------------------------------

}  // namespace test
}  // namespace atlas

int main(int argc, char** argv) {
    return atlas::test::run(argc, argv);
}