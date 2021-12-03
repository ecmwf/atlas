/*
 * (C) Copyright 2020 UCAR
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

class NodalCubedSphereMeshGenerator : public MeshGenerator::Implementation {
public:
    NodalCubedSphereMeshGenerator(const eckit::Parametrisation& = util::NoConfig());

    virtual void generate(const Grid&, const grid::Distribution&, Mesh&) const override;
    virtual void generate(const Grid&, Mesh&) const override;

    using MeshGenerator::Implementation::generate;

    static std::string static_type() { return "nodal-cubedsphere"; }
    std::string type() const override { return static_type(); }

private:
    virtual void hash(eckit::Hash&) const override;

    void configure_defaults();

private:
    util::Metadata options;
};

//--------------------------------------------------------------------------------------------------

}  // namespace meshgenerator
}  // namespace atlas
