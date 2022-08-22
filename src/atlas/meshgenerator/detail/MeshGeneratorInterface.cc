/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "atlas/meshgenerator/detail/MeshGeneratorInterface.h"
#include "atlas/grid/Distribution.h"
#include "atlas/grid/Grid.h"
#include "atlas/grid/Partitioner.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/meshgenerator.h"
#include "atlas/meshgenerator/detail/MeshGeneratorImpl.h"
#include "atlas/runtime/Exception.h"

using atlas::Mesh;

namespace atlas {
namespace meshgenerator {

//----------------------------------------------------------------------------------------------------------------------

extern "C" {

void atlas__MeshGenerator__delete(MeshGenerator::Implementation* This) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialisd atlas_MeshGenerator");
    delete This;
}

const MeshGenerator::Implementation* atlas__MeshGenerator__create_noconfig(const char* name) {
    const MeshGenerator::Implementation* meshgenerator(nullptr);
    {
        MeshGenerator m(std::string{name});
        meshgenerator = m.get();
        meshgenerator->attach();
    }
    meshgenerator->detach();
    return meshgenerator;
}

const MeshGenerator::Implementation* atlas__MeshGenerator__create(const char* name,
                                                                  const eckit::Parametrisation* config) {
    const MeshGenerator::Implementation* meshgenerator(nullptr);
    ATLAS_ASSERT(config);
    {
        MeshGenerator m(std::string(name), *config);
        meshgenerator = m.get();
        meshgenerator->attach();
    }
    meshgenerator->detach();
    return meshgenerator;
}

Mesh::Implementation* atlas__MeshGenerator__generate__grid_griddist(
    const MeshGenerator::Implementation* This, const Grid::Implementation* grid,
    const grid::Distribution::Implementation* distribution) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialisd atlas_MeshGenerator");
    ATLAS_ASSERT(grid != nullptr, "Cannot access uninitialisd atlas_Grid");
    ATLAS_ASSERT(distribution != nullptr, "Cannot access uninitialisd atlas_GridDistribution");

    Mesh::Implementation* m;
    {
        Mesh mesh = This->generate(Grid(grid), grid::Distribution(distribution));
        mesh.get()->attach();
        m = mesh.get();
    }
    m->detach();
    return m;
}

Mesh::Implementation* atlas__MeshGenerator__generate__grid(const MeshGenerator::Implementation* This,
                                                           const Grid::Implementation* grid) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialisd atlas_MeshGenerator");
    ATLAS_ASSERT(grid != nullptr, "Cannot access uninitialisd atlas_Grid");
    Mesh::Implementation* m;
    {
        Mesh mesh = This->generate(Grid(grid));
        ;
        mesh.get()->attach();
        m = mesh.get();
    }
    m->detach();
    return m;
}

Mesh::Implementation* atlas__MeshGenerator__generate__grid_partitioner(
    const MeshGenerator::Implementation* This, const Grid::Implementation* grid,
    const grid::Partitioner::Implementation* partitioner) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_MeshGenerator");
    ATLAS_ASSERT(grid != nullptr, "Cannot access uninitialised atlas_Grid");
    ATLAS_ASSERT(partitioner != nullptr, "Cannot access uninitialised atlas_Partitioner");

    Mesh::Implementation* m;
    {
        Mesh mesh = This->generate(Grid(grid), grid::Partitioner(partitioner));
        mesh.get()->attach();
        m = mesh.get();
    }
    m->detach();
    return m;
}
}  // extern "C"

//----------------------------------------------------------------------------------------------------------------------

}  // namespace meshgenerator

//----------------------------------------------------------------------------------------------------------------------

}  // namespace atlas
