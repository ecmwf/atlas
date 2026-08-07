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

void filter_cutoff(const functionspace::Spectral& spectral, Field& field, int cutoff);

}  // namespace spectral
}  // namespace atlas

extern "C" {
void atlas__spectral__filter_cutoff(const atlas::functionspace::detail::Spectral* spectral, atlas::field::FieldImpl* field,
									int cutoff);
}
