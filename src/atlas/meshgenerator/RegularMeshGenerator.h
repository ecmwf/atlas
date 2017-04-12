
#pragma once

#include "atlas/meshgenerator/MeshGenerator.h"
#include "atlas/util/Metadata.h"
#include "atlas/util/Config.h"

namespace eckit { class Parametrisation; }

namespace atlas {
    class Mesh;
}

namespace atlas {
namespace grid {
    //class Regular;
    class Distribution;
} }

namespace atlas {
namespace meshgenerator {

//----------------------------------------------------------------------------------------------------------------------

class RegularMeshGenerator : public MeshGenerator::meshgenerator_t {

public:

    RegularMeshGenerator(const eckit::Parametrisation& = util::NoConfig() );

    virtual void generate(const atlas::grid::Grid&, const grid::Distribution&, Mesh&) const;
    virtual void generate(const atlas::grid::Grid&, Mesh&) const;

    using MeshGenerator::meshgenerator_t::generate;

private:

    virtual void hash(eckit::MD5&) const;

    void configure_defaults();

    void generate_mesh(
      const atlas::grid::RegularGrid&,
      const std::vector<int>& parts,
      Mesh& m ) const;

private:

    util::Metadata options;

};

//----------------------------------------------------------------------------------------------------------------------

} // namespace meshgenerator
} // namespace atlas
