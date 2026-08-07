/*
 * (C) Copyright 2026 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <array>
#include <complex>
#include <vector>

#include "eckit/types/FloatCompare.h"

#include "atlas/algorithms/spectral/filter_cutoff.h"
#include "atlas/algorithms/spectral/power_spectrum.h"
#include "atlas/array/MakeView.h"
#include "atlas/field/Field.h"
#include "atlas/functionspace/Spectral.h"
#include "atlas/library/config.h"
#include "atlas/mdspan.h"
#include "atlas/option.h"
#include "atlas/parallel/mpi/mpi.h"

#include "tests/AtlasTestEnvironment.h"

#if ATLAS_HAVE_TRANS
#if ATLAS_HAVE_ECTRANS
#include "ectrans/transi.h"
#else
#include "transi/trans.h"
#endif
#endif

namespace atlas {
namespace test {

//-----------------------------------------------------------------------------

struct AtlasSpectralPowerSpectrumEnvironment : public AtlasTestEnvironment {
    AtlasSpectralPowerSpectrumEnvironment(int argc, char* argv[]): AtlasTestEnvironment(argc, argv) {
        if (mpi::comm().size() == 1) {
            trans_use_mpi(false);
        }
        trans_init();
    }

    ~AtlasSpectralPowerSpectrumEnvironment() { trans_finalize(); }
};

namespace {

template <typename Value>
Value spherical_harmonics_power_metric(int n) {
    const Value earth_radius = 6371.229e3;
    const Value zlam         = 0.5;
    const Value zfact        = (n == 0 ? 1. : earth_radius * earth_radius / static_cast<Value>(n * (n + 1)));
    return zlam * zfact;
}

std::complex<double> spectral_function(int level, int n, int m) {
    const double re = 100. * static_cast<double>(level + 1) + 10. * static_cast<double>(m) + static_cast<double>(n);
    const double im = (m == 0 ? 0. : -50. * static_cast<double>(level + 1) - 5. * static_cast<double>(m) + 0.5 * static_cast<double>(n));
    return {re, im};
}

double coefficient_tolerance(double) {
    return 1.e-12;
}

float coefficient_tolerance(float) {
    return 1.e-6f;
}

template <typename Value>
Field make_spectral_test_field(const functionspace::Spectral& spectral, int levels, const char* name) {
    Field field = spectral.createField<Value>(option::name(name) | option::levels(levels));
    auto coefficients = array::make_view<Value, 2>(field);

    coefficients.assign(Value{0});
    spectral.parallel_for([&](idx_t real, idx_t imag, int n, int m) {
        for (int level = 0; level < levels; ++level) {
            const auto value = spectral_function(level, n, m);
            coefficients(real, level) = static_cast<Value>(value.real());
            coefficients(imag, level) = static_cast<Value>(value.imag());
        }
    });

    return field;
}

template <typename CoefficientValue>
std::vector<double> expected_spherical_harmonics_power_spectrum(int truncation, int levels) {
    std::vector<double> spectrum(static_cast<std::size_t>((truncation + 1) * levels), 0.);
    mdspan<double, dims<2>> spectrum_view{spectrum.data(), static_cast<size_t>(truncation + 1), static_cast<size_t>(levels)};

    for (int n = 0; n <= truncation; ++n) {
        const double metric = spherical_harmonics_power_metric<double>(n);
        for (int m = 0; m <= n; ++m) {
            const double multiplicity = (m == 0 ? 1. : 2.);
            for (int level = 0; level < levels; ++level) {
                const auto value = spectral_function(level, n, m);
                const double re  = static_cast<double>(static_cast<CoefficientValue>(value.real()));
                const double im  = static_cast<double>(static_cast<CoefficientValue>(value.imag()));
                spectrum_view(n, level) += multiplicity * (re * re + im * im);
            }
        }
        for (int level = 0; level < levels; ++level) {
            spectrum_view(n, level) *= metric;
        }
    }

    return spectrum;
}

template <typename Value>
void expect_spherical_harmonics_power_spectrum(const functionspace::Spectral& spectral, const Field& field, int truncation, int levels) {
    const std::vector<double> spectrum = spectral::power_spectrum(spectral, field);
    const std::vector<double> expected = expected_spherical_harmonics_power_spectrum<Value>(truncation, levels);

    EXPECT_EQ(spectrum.size(), expected.size());
    for (std::size_t j = 0; j < spectrum.size(); ++j) {
        EXPECT_APPROX_EQ(spectrum[j], expected[j], 1.e-12);
    }
}

template <typename Value>
void expect_spectral_cutoff(const functionspace::Spectral& spectral, const Field& field, int truncation, int levels, int cutoff) {
    auto coefficients = array::make_view<Value, 2>(field);
    const Value tolerance = coefficient_tolerance(Value{});

    EXPECT(field.dirty());
    spectral.parallel_for([&](idx_t real, idx_t imag, int n, int m) {
        for (int level = 0; level < levels; ++level) {
            const auto value = (n > cutoff ? std::complex<double>{0., 0.} : spectral_function(level, n, m));
            EXPECT_APPROX_EQ(coefficients(real, level), static_cast<Value>(value.real()), tolerance);
            EXPECT_APPROX_EQ(coefficients(imag, level), static_cast<Value>(value.imag()), tolerance);
        }
    });

    const std::vector<double> spectrum = spectral::power_spectrum(spectral, field);
    const std::vector<double> unfiltered_expected = expected_spherical_harmonics_power_spectrum<Value>(truncation, levels);
    mdspan<const double, dims<2>> spectrum_view{spectrum.data(), static_cast<size_t>(truncation + 1), static_cast<size_t>(levels)};
    mdspan<const double, dims<2>> expected_view{unfiltered_expected.data(), static_cast<size_t>(truncation + 1), static_cast<size_t>(levels)};

    for (int n = 0; n <= truncation; ++n) {
        for (int level = 0; level < levels; ++level) {
            const double expected = (n > cutoff ? 0. : expected_view(n, level));
            EXPECT_APPROX_EQ(spectrum_view(n, level), expected, 1.e-12);
        }
    }
}

}  // namespace

//-----------------------------------------------------------------------------

CASE("test_spectral_field_stores_real_and_imaginary_slots") {
    const int truncation = 9;
    const idx_t expected_global_coefficients = (truncation + 1) * (truncation + 2);

    functionspace::Spectral spectral(truncation);
    EXPECT_EQ(spectral.nb_spectral_coefficients_global(), expected_global_coefficients);

    idx_t parallel_for_count = 0;
    idx_t m0_count           = 0;
    spectral.parallel_for([&](idx_t real, idx_t imag, int, int m) {
        EXPECT_EQ(imag, real + 1);
        ++parallel_for_count;
        if (m == 0) {
            ++m0_count;
        }
    });
    mpi::comm().allReduceInPlace(parallel_for_count, eckit::mpi::sum());
    mpi::comm().allReduceInPlace(m0_count, eckit::mpi::sum());
    EXPECT_EQ(2 * parallel_for_count, expected_global_coefficients);
    EXPECT_EQ(m0_count, truncation + 1);
}

CASE("test_spectral_power_spectrum") {
    const int truncation = 9;
    const int levels     = 3;

    functionspace::Spectral spectral(truncation);

    const Field field = make_spectral_test_field<double>(spectral, levels, "spherical_harmonics");
    expect_spherical_harmonics_power_spectrum<double>(spectral, field, truncation, levels);

    const Field float_field = make_spectral_test_field<float>(spectral, levels, "spherical_harmonics_float");
    expect_spherical_harmonics_power_spectrum<float>(spectral, float_field, truncation, levels);
}

CASE("test_spectral_filter") {
    const int truncation = 9;
    const int levels     = 3;
    const int cutoff     = 4;

    functionspace::Spectral spectral(truncation);

    Field field = make_spectral_test_field<double>(spectral, levels, "spherical_harmonics");
    spectral::filter_cutoff(spectral, field, cutoff);
    expect_spectral_cutoff<double>(spectral, field, truncation, levels, cutoff);

    Field float_field = make_spectral_test_field<float>(spectral, levels, "spherical_harmonics_float");
    spectral::filter_cutoff(spectral, float_field, cutoff);
    expect_spectral_cutoff<float>(spectral, float_field, truncation, levels, cutoff);
}

//-----------------------------------------------------------------------------

}  // namespace test
}  // namespace atlas

int main(int argc, char** argv) {
    return atlas::test::run<atlas::test::AtlasSpectralPowerSpectrumEnvironment>(argc, argv);
}
