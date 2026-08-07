/*
 * (C) Copyright 2026 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "atlas/algorithms/spectral/power_spectrum.h"

#include <algorithm>
#include <array>
#include <vector>

#include "atlas/array/MakeView.h"
#include "atlas/field/detail/FieldImpl.h"
#include "atlas/functionspace/FunctionSpace.h"
#include "atlas/mdspan.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/runtime/Exception.h"

namespace atlas {
namespace spectral {
namespace {

template <typename Value>
Value spherical_harmonics_power_metric(int n) {
    const Value earth_radius = 6371.229e3;
    const Value zlam         = 0.5;
    const Value zfact        = (n == 0 ? 1. : earth_radius * earth_radius / static_cast<Value>(n * (n + 1)));
    return zlam * zfact;
}

template <typename Value>
void apply_spherical_harmonics_power_metric(int truncation, int levels, Value* spectrum) {
    mdspan<Value, dims<2>> spectrum_view{spectrum, static_cast<size_t>(truncation + 1), static_cast<size_t>(levels)};
    for (int n = 0; n <= truncation; ++n) {
        const Value metric = spherical_harmonics_power_metric<Value>(n);
        for (int level = 0; level < levels; ++level) {
            spectrum_view(n, level) *= metric;
        }
    }
}

template <typename Value>
void power_spectrum_impl(const functionspace::Spectral& spectral, const Field& field, Value spectrum_raw[],
                         std::array<size_t, 2> spectrum_extents) {
    const int levels = field.levels();

    if (spectrum_extents[0] != static_cast<size_t>(spectral.truncation() + 1)) {
        throw_Exception("Power spectrum wave-number dimension does not match spectral truncation", Here());
    }
    if (spectrum_extents[1] != static_cast<size_t>(std::max(1, levels))) {
        throw_Exception("Power spectrum level dimension does not match spectral field levels", Here());
    }

    mdspan<Value, dims<2>> spectrum{spectrum_raw, spectrum_extents};
    std::fill(spectrum.data_handle(), spectrum.data_handle() + spectrum.size(), Value{});

    auto implementation_for_datatype = [&](auto coefficient_value) {
        using CoefficientValue = decltype(coefficient_value);

        if (levels == 0) {
            auto coefficients = array::make_view<CoefficientValue, 1>(field);
            spectral.parallel_for([&](idx_t real, idx_t imag, int n, int m) {
                const CoefficientValue scale = (m == 0 ? CoefficientValue{1} : CoefficientValue{2});
                const CoefficientValue re    = coefficients(real);
                const CoefficientValue im    = (m == 0 ? CoefficientValue{0} : coefficients(imag));
                spectrum(n, 0) += static_cast<Value>(scale * (re * re + im * im));
            });
        }
        else {
            auto coefficients = array::make_view<CoefficientValue, 2>(field);
            spectral.parallel_for([&](idx_t real, idx_t imag, int n, int m) {
                const CoefficientValue scale = (m == 0 ? CoefficientValue{1} : CoefficientValue{2});
                for (int level = 0; level < levels; ++level) {
                    const CoefficientValue re = coefficients(real, level);
                    const CoefficientValue im = (m == 0 ? CoefficientValue{0} : coefficients(imag, level));
                    spectrum(n, level) += static_cast<Value>(scale * (re * re + im * im));
                }
            });
        }
    };

    switch (field.datatype().kind()) {
        case array::DataType::KIND_REAL32: implementation_for_datatype(float{}); break;
        case array::DataType::KIND_REAL64: implementation_for_datatype(double{}); break;
        default:
            throw_Exception("Power spectrum only supports real32 and real64 spectral fields", Here());
    }

    mpi::comm().allReduceInPlace(spectrum.data_handle(), spectrum.data_handle() + spectrum.size(), eckit::mpi::sum());
    apply_spherical_harmonics_power_metric(spectral.truncation(), std::max(1, levels), spectrum.data_handle());
}

}  // namespace

void power_spectrum(const functionspace::Spectral& spectral, const Field& field, double spectrum[],
                    std::array<size_t, 2> spectrum_extents) {
    power_spectrum_impl(spectral, field, spectrum, spectrum_extents);
}

void power_spectrum(const functionspace::Spectral& spectral, const Field& field, double spectrum[], size_t spectrum_size) {
    constexpr size_t one_level = 1;
    power_spectrum(spectral, field, spectrum, {spectrum_size, one_level});
}

std::vector<double> power_spectrum(const functionspace::Spectral& spectral, const Field& field) {
    const size_t levels     = static_cast<size_t>(std::max(1, field.levels()));
    const size_t truncation = static_cast<size_t>(spectral.truncation());
    std::vector<double> spectrum(levels * (truncation + 1));
    power_spectrum(spectral, field, spectrum.data(), {truncation + 1, levels});
    return spectrum;
}

}  // namespace spectral
}  // namespace atlas

extern "C" {
void atlas__spectral__power_spectrum(const atlas::functionspace::detail::Spectral* spectral,
                                     const atlas::field::FieldImpl* field, double spectrum[], int spectrum_extents[]) {
    ATLAS_ASSERT(spectral != nullptr);
    ATLAS_ASSERT(field != nullptr);
    ATLAS_ASSERT(spectrum != nullptr);
    ATLAS_ASSERT(spectrum_extents != nullptr);

    atlas::spectral::power_spectrum(atlas::functionspace::Spectral(atlas::FunctionSpace(spectral)), atlas::Field(field),
                                    spectrum,
                                    {static_cast<size_t>(spectrum_extents[0]), static_cast<size_t>(spectrum_extents[1])});
}
}
