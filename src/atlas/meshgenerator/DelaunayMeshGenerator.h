/*
 * (C) Copyright 1996-2017 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#pragma once

#include "atlas/meshgenerator/MeshGenerator.h"

namespace atlas {
namespace mesh {
class Mesh;
}
}

namespace atlas {
namespace grid {
class Grid;
}
}

namespace atlas {
namespace meshgenerator {

//----------------------------------------------------------------------------------------------------------------------

class DelaunayMeshGenerator : public MeshGenerator {
public:

  DelaunayMeshGenerator();
  DelaunayMeshGenerator(const eckit::Parametrisation& p);

  virtual ~DelaunayMeshGenerator();

private: // methods

  virtual void hash(eckit::MD5&) const;

  virtual void generate(const grid::Grid& g, const grid::Distribution&, mesh::Mesh& mesh) const;
  virtual void generate(const grid::Grid& g, mesh::Mesh& mesh) const;
};

//----------------------------------------------------------------------------------------------------------------------

}  // namespace meshgenerator
}  // namespace atlas
