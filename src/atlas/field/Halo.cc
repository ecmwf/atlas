/*
 * (C) Copyright 2025 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "atlas/library/config.h"
#include "atlas/util/Metadata.h"
#include "atlas/field/detail/FieldImpl.h"
#include "atlas/functionspace/FunctionSpace.h"

namespace atlas {
namespace field {

Halo::Halo(field::FieldImpl& f) :
    functionspace::HaloDescription(f.functionspace().halo_description()),
    field_(f) {
}

bool Halo::updated() const {
    return field_.metadata().getBool("halo_updated", false);
}

void Halo::updated(bool v) {
    field_.metadata().set("halo_updated", v);
}

void Halo::update() {
    field_.haloExchange();
}

void Halo::update(const eckit::Parametrisation& config) {
    std::string execution_space{"host"};
    config.get("execution_space", execution_space);
    bool on_device = (execution_space == "device");
    field_.haloExchange(on_device);
}

bool Halo::appended() const {
    auto size = field_.functionspace().size();
    return contiguous() && (end() == size);
}

// ------------------------------------------------------------------

}  // namespace field
}  // namespace atlas
