/*
 * (C) Crown Copyright 2021 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include "atlas/array.h"
#include "atlas/meshgenerator/MeshGenerator.h"
#include "atlas/meshgenerator/detail/MeshGeneratorImpl.h"
#include "atlas/util/Config.h"
#include "atlas/util/Metadata.h"

#ifndef DOXYGEN_SHOULD_SKIP_THIS
namespace eckit {
class Parametrisation;
}

namespace atlas {
class CubedSphereGrid;
class Mesh;
template <typename T>
class vector;
}  // namespace atlas

namespace atlas {
namespace grid {
class Distribution;
}  // namespace grid
}  // namespace atlas
#endif

namespace atlas {
namespace meshgenerator {

//--------------------------------------------------------------------------------------------------

class CubedSphereDualMeshGenerator : public MeshGenerator::Implementation {
public:
    CubedSphereDualMeshGenerator(const eckit::Parametrisation& = util::NoConfig());

    virtual void generate(const Grid&, const grid::Distribution&, Mesh&) const override;
    virtual void generate(const Grid&, Mesh&) const override;

    using MeshGenerator::Implementation::generate;

    static std::string static_type() { return "cubedsphere_dual"; }
    std::string type() const override { return static_type(); }

private:
    virtual void hash(eckit::Hash&) const override;

    void configure_defaults();

    void generate_mesh(const CubedSphereGrid&, const grid::Distribution&, Mesh&) const;

    void set_metadata(Mesh&) const;

private:
    util::Metadata options;
};

//--------------------------------------------------------------------------------------------------

}  // namespace meshgenerator
}  // namespace atlas
