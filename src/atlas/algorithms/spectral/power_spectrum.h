/*
 * (C) Copyright 2026 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#pragma once

#include <array>
#include <cstddef>
#include <vector>

#include "atlas/field/Field.h"
#include "atlas/functionspace/Spectral.h"

namespace atlas {
namespace field {
class FieldImpl;
}
namespace functionspace {
namespace detail {
class Spectral;
}
}  // namespace functionspace
}  // namespace atlas

namespace atlas {
namespace spectral {

void power_spectrum(const functionspace::Spectral& spectral, const Field& field, double spectrum[],
                    std::array<size_t, 2> spectrum_extents);

void power_spectrum(const functionspace::Spectral& spectral, const Field& field, double spectrum[], size_t spectrum_size);

std::vector<double> power_spectrum(const functionspace::Spectral& spectral, const Field& field);

}  // namespace spectral
}  // namespace atlas

extern "C" {
void atlas__spectral__power_spectrum(const atlas::functionspace::detail::Spectral* spectral,
                                     const atlas::field::FieldImpl* field, double spectrum[], int spectrum_extents[]);
}
