/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "atlas/mesh/Mesh.h"
#include "atlas/grid/Grid.h"
#include "atlas/grid/Partitioner.h"
#include "atlas/meshgenerator/MeshGenerator.h"

namespace atlas {

//----------------------------------------------------------------------------------------------------------------------

Mesh::Mesh(): Handle(new Implementation()) {}

Mesh::Mesh(const Grid& grid):
    Handle([&]() {
        auto meshgenerator = MeshGenerator{grid.meshgenerator()};
        auto mesh          = meshgenerator.generate(grid, grid::Partitioner(grid.partitioner()));
        mesh.get()->attach();
        return mesh.get();
    }()) {
    get()->detach();
}

Mesh::Mesh(const Grid& grid, const grid::Partitioner& partitioner):
    Handle([&]() {
        auto meshgenerator = MeshGenerator{grid.meshgenerator()};
        auto mesh          = meshgenerator.generate(grid, partitioner);
        mesh.get()->attach();
        return mesh.get();
    }()) {
    get()->detach();
}


Mesh::Mesh(eckit::Stream& stream): Handle(new Implementation(stream)) {}

//----------------------------------------------------------------------------------------------------------------------

}  // namespace atlas
