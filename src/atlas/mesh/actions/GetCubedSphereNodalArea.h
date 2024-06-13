/*
 * (C) Crown Copyright 2024 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include "atlas/field/Field.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/mesh/Nodes.h"


namespace atlas {

class Mesh;
class Field;

namespace mesh {
namespace actions {


/// Provide the area around nodes of cubed sphere mesh
class GetCubedSphereNodalArea {
public:
    Field& operator()(Mesh&);
};

}  // namespace actions
}  // namespace mesh
}  // namespace atlas
