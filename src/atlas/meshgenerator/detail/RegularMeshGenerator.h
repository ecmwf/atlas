
#pragma once

#include "atlas/meshgenerator/MeshGenerator.h"
#include "atlas/meshgenerator/detail/MeshGeneratorImpl.h"
#include "atlas/util/Config.h"
#include "atlas/util/Metadata.h"

namespace eckit {
class Parametrisation;
}

namespace atlas {
class RegularGrid;
class Mesh;
}  // namespace atlas

namespace atlas {
namespace grid {
class Distribution;
}  // namespace grid
}  // namespace atlas

namespace atlas {
namespace meshgenerator {

//----------------------------------------------------------------------------------------------------------------------

class RegularMeshGenerator : public MeshGenerator::Implementation {
public:
    RegularMeshGenerator( const eckit::Parametrisation& = util::NoConfig() );

    virtual void generate( const Grid&, const grid::Distribution&, Mesh& ) const override;
    virtual void generate( const Grid&, Mesh& ) const override;

    using MeshGenerator::Implementation::generate;

private:
    virtual void hash( eckit::Hash& ) const override;

    void configure_defaults();

    void generate_mesh( const RegularGrid&, const std::vector<int>& parts, Mesh& m ) const;

private:
    util::Metadata options;
};

//----------------------------------------------------------------------------------------------------------------------

}  // namespace meshgenerator
}  // namespace atlas
