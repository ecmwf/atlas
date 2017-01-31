/*
 * (C) Copyright 1996-2017 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef atlas_meshgen_Delaunay_h
#define atlas_meshgen_Delaunay_h

#include "atlas/mesh/generators/MeshGenerator.h"

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
namespace mesh {
namespace generators {

//----------------------------------------------------------------------------------------------------------------------

class Delaunay : public MeshGenerator {
public:

  Delaunay();
  Delaunay(const eckit::Parametrisation& p);

  virtual ~Delaunay();

private: // methods

  virtual void hash(eckit::MD5&) const;

  virtual void generate(const grid::Grid& g, const grid::GridDistribution&, Mesh& mesh) const;
  virtual void generate(const grid::Grid& g, Mesh& mesh) const;
};

//----------------------------------------------------------------------------------------------------------------------

}  // namespace generators
}  // namespace mesh
}  // namespace atlas

#endif
