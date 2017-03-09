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
#include "atlas/util/Metadata.h"
#include "atlas/util/Config.h"

namespace eckit { class Parametrisation; }

namespace atlas {
namespace mesh {
    class Mesh;
} }

namespace atlas {
namespace grid {
    class StructuredGrid;
    class Distribution;
} }


namespace atlas {
namespace meshgenerator {

struct Region;

//----------------------------------------------------------------------------------------------------------------------

class StructuredMeshGenerator : public MeshGenerator {

public:

  StructuredMeshGenerator(const eckit::Parametrisation& = util::NoConfig() );

  virtual void generate(const grid::Grid&, const grid::Distribution&, mesh::Mesh&) const;
  virtual void generate(const grid::Grid&, mesh::Mesh&) const;

  using MeshGenerator::generate;

private:

  virtual void hash(eckit::MD5&) const;

  void configure_defaults();

  void generate_region(
    const grid::StructuredGrid&,
    const std::vector<int>& parts,
    int mypart,
    Region& region) const;

  void generate_mesh_new(
    const grid::StructuredGrid&,
    const std::vector<int>& parts,
    const Region& region,
    mesh::Mesh& m) const;

  void generate_mesh(
    const grid::StructuredGrid&,
    const std::vector<int>& parts,
    const Region& region,
    mesh::Mesh& m ) const;

  void generate_global_element_numbering(
    mesh::Mesh& mesh ) const;

private:

  util::Metadata options;

};

//----------------------------------------------------------------------------------------------------------------------

} // namespace meshgenerator
} // namespace atlas
