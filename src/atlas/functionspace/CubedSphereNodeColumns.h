/*
 * (C) Crown Copyright 2021 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include "atlas/functionspace/NodeColumns.h"
#include "atlas/functionspace/detail/CubedSphereUtility.h"
#include "atlas/mesh.h"
#include "atlas/mesh/Nodes.h"

namespace atlas {
namespace functionspace {

class CubedSphereNodeColumns : public functionspace::NodeColumns,
  public detail::CubedSphereUtility {

public:

  inline CubedSphereNodeColumns() : NodeColumns(), CubedSphereUtility() {}

  inline CubedSphereNodeColumns(const FunctionSpace& functionSpace) :
    NodeColumns(functionSpace), CubedSphereUtility(
      this->mesh().nodes().field("ijt"),
      this->mesh().nodes().ghost()) {}

  inline CubedSphereNodeColumns(const Mesh& mesh,
    const eckit::Configuration& configuration) :
    NodeColumns(mesh, configuration), CubedSphereUtility(
      this->mesh().nodes().field("ijt"),
      this->mesh().nodes().ghost()) {}

  inline CubedSphereNodeColumns(const Mesh& mesh) :
    NodeColumns(mesh), CubedSphereUtility(
      this->mesh().nodes().field("ijt"),
      this->mesh().nodes().ghost()) {}

};

} // namespace functionspace
} // namespace atlas
