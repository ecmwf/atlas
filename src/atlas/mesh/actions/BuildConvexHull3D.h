/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#pragma once

#include "eckit/config/Parametrisation.h"

#include "atlas/util/Config.h"

namespace atlas {

class Mesh;

namespace mesh {
namespace actions {

/// Creates a 3D convex-hull on the mesh points
class BuildConvexHull3D {
public:
    BuildConvexHull3D(const eckit::Parametrisation& = util::NoConfig());
    void operator()(Mesh&) const;
private:
    bool remove_duplicate_points_;
};

}  // namespace actions
}  // namespace mesh
}  // namespace atlas
