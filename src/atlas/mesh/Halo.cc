/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "atlas/mesh/Halo.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/runtime/Exception.h"
#include "atlas/util/Metadata.h"

namespace atlas {
namespace mesh {

Halo::Halo(const Mesh& mesh) {
    size_ = 0;
    mesh.metadata().get("halo", size_);
}

Halo::Halo(const detail::MeshImpl& mesh) {
    size_ = 0;
    mesh.metadata().get("halo", size_);
}

int Halo::size() const {
    ATLAS_ASSERT(size_ >= 0);
    return size_;
}

}  // namespace mesh
}  // namespace atlas
