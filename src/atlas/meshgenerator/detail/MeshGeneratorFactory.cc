/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <string>

#include "atlas/meshgenerator/detail/CubedSphereDualMeshGenerator.h"
#include "atlas/meshgenerator/detail/CubedSphereMeshGenerator.h"
#include "atlas/meshgenerator/detail/DelaunayMeshGenerator.h"
#include "atlas/meshgenerator/detail/HealpixMeshGenerator.h"
#include "atlas/meshgenerator/detail/MeshGeneratorFactory.h"
#include "atlas/meshgenerator/detail/NodalCubedSphereMeshGenerator.h"
#include "atlas/meshgenerator/detail/StructuredMeshGenerator.h"

using atlas::Mesh;

namespace atlas {
namespace meshgenerator {

//----------------------------------------------------------------------------------------------------------------------

void force_link() {
    static struct Link {
        Link() {
            MeshGeneratorBuilder<meshgenerator::CubedSphereMeshGenerator>();
            MeshGeneratorBuilder<meshgenerator::CubedSphereDualMeshGenerator>();
            MeshGeneratorBuilder<meshgenerator::StructuredMeshGenerator>();
            MeshGeneratorBuilder<meshgenerator::DelaunayMeshGenerator>();
            MeshGeneratorBuilder<meshgenerator::HealpixMeshGenerator>();
            MeshGeneratorBuilder<meshgenerator::NodalCubedSphereMeshGenerator>();
        }
    } link;
}

//----------------------------------------------------------------------------------------------------------------------

const MeshGenerator::Implementation* MeshGeneratorFactory::build(const std::string& builder) {
    return build(builder, util::NoConfig());
}

const MeshGenerator::Implementation* MeshGeneratorFactory::build(const std::string& builder,
                                                                 const eckit::Parametrisation& param) {
    force_link();
    auto factory = get(builder);
    return factory->make(param);
}

//----------------------------------------------------------------------------------------------------------------------

}  // namespace meshgenerator
}  // namespace atlas
