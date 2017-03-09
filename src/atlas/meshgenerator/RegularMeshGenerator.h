
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
    //class Regular;
    class Distribution;
} }

namespace atlas {
namespace meshgenerator {

//----------------------------------------------------------------------------------------------------------------------

class RegularMeshGenerator : public MeshGenerator {

public:

    RegularMeshGenerator(const eckit::Parametrisation& = util::NoConfig() );

    virtual void generate(const atlas::grid::Grid&, const grid::Distribution&, mesh::Mesh&) const;
    virtual void generate(const atlas::grid::Grid&, mesh::Mesh&) const;

    using MeshGenerator::generate;

private:

    virtual void hash(eckit::MD5&) const;

    void configure_defaults();

    void generate_mesh(
      const atlas::grid::RegularGrid&,
      const std::vector<int>& parts,
      mesh::Mesh& m ) const;

    void generate_global_element_numbering(
      mesh::Mesh& mesh ) const;

private:

    util::Metadata options;

};

//----------------------------------------------------------------------------------------------------------------------

} // namespace meshgenerator
} // namespace atlas
