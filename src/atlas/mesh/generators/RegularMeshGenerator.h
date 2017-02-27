
#ifndef atlas_mesh_generators_RegularMeshGenerator_h
#define atlas_mesh_generators_RegularMeshGenerator_h

#include "atlas/mesh/generators/MeshGenerator.h"
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

//----------------------------------------------------------------------------------------------------------------------

class RegularMeshGenerator : public MeshGenerator {

public:

    RegularMeshGenerator(const eckit::Parametrisation& = util::NoConfig() );

    virtual void generate(const atlas::grid::Grid&, const grid::GridDistribution&, Mesh&) const;
    virtual void generate(const atlas::grid::Grid&, Mesh&) const;

    using MeshGenerator::generate;

private:

    virtual void hash(eckit::MD5&) const;

    void configure_defaults();

    void generate_mesh(
      const atlas::grid::Regular&,
      const std::vector<int>& parts,
      Mesh& m ) const;

    void generate_global_element_numbering(
      Mesh& mesh ) const;

private:

    util::Metadata options;

};

//----------------------------------------------------------------------------------------------------------------------

} // namespace generators
} // namespace mesh
} // namespace atlas

#endif // atlas_mesh_generators_Structured_h
