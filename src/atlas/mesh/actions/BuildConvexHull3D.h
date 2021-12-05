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

namespace atlas {

class Mesh;

namespace mesh {
namespace actions {

/// Creates a 3D convex-hull on the mesh points
class BuildConvexHull3D {
public:
    void operator()(Mesh&) const;
};

}  // namespace actions
}  // namespace mesh
}  // namespace atlas
