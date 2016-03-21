/*
 * (C) Copyright 1996-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef atlas_mesh_generators_Structured_h
#define atlas_mesh_generators_Structured_h

#include "atlas/mesh/generators/MeshGenerator.h"
#include "atlas/util/Metadata.h"

namespace eckit { class Parametrisation; }

namespace atlas {
namespace mesh {
    class Mesh;
} }

namespace atlas {
namespace grid {
namespace global {
    class Structured;
}
    class GridDistribution;
} }

namespace atlas {
namespace mesh {
namespace generators {
    struct Region;
} } }

namespace atlas {
namespace mesh {
namespace generators {

// -----------------------------------------------------------------------------

class Structured : public MeshGenerator {

public:

  Structured();
  Structured(const eckit::Parametrisation&);

  virtual void generate(const grid::Grid&, const grid::GridDistribution&, Mesh&) const;
  virtual void generate(const grid::Grid&, Mesh&) const;
  using MeshGenerator::generate;

private:

  void configure_defaults();

  void generate_region(
    const grid::global::Structured&,
    const std::vector<int>& parts,
    int mypart,
    Region& region) const;

  void generate_mesh_new(
    const grid::global::Structured&,
    const std::vector<int>& parts,
    const Region& region,
    Mesh& m) const;

  void generate_mesh(
    const grid::global::Structured&,
    const std::vector<int>& parts,
    const Region& region,
    Mesh& m ) const;

  void generate_global_element_numbering(
    Mesh& mesh ) const;

public:

  util::Metadata options;

};

// -----------------------------------------------------------------------------

} // namespace generators
} // namespace mesh
} // namespace atlas

#endif // atlas_mesh_generators_Structured_h
