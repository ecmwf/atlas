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

#include "atlas/grid/Distribution.h"
#include "atlas/grid/Grid.h"
#include "atlas/grid/Partitioner.h"
#include "atlas/mesh/Nodes.h"
#include "atlas/meshgenerator/MeshGenerator.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/runtime/Exception.h"

namespace atlas {

//----------------------------------------------------------------------------------------------------------------------

Mesh::Mesh(): Handle(new Implementation()) {}

Mesh::Mesh(const Grid& grid, const eckit::Configuration& config):
    Handle([&]() {
        if (config.has("mpi_comm")) {
            mpi::push(config.getString("mpi_comm"));
        }
        auto cfg           = grid.meshgenerator() | util::Config(config);
        auto meshgenerator = MeshGenerator{grid.meshgenerator() | config};
        auto mesh          = meshgenerator.generate(grid, grid::Partitioner(grid.partitioner() | config));
        if (config.has("mpi_comm")) {
            mpi::pop();
        }
        mesh.get()->attach();
        return mesh.get();
    }()) {
    get()->detach();
}

Mesh::Mesh(const Grid& grid, const grid::Partitioner& partitioner, const eckit::Configuration& config):
    Handle([&]() {
        auto mpi_comm = partitioner.mpi_comm();
        if (config.has("mpi_comm")) {
            mpi_comm = config.getString("mpi_comm");
            ATLAS_ASSERT(mpi_comm == partitioner.mpi_comm());
        }
        mpi::Scope mpi_scope(mpi_comm);
        auto meshgenerator = MeshGenerator{grid.meshgenerator() | config};
        auto mesh          = meshgenerator.generate(grid, partitioner);
        mesh.get()->attach();
        return mesh.get();
    }()) {
    get()->detach();
}

Mesh::Mesh(const Grid& grid, const grid::Distribution& distribution, const eckit::Configuration& config):
    Handle([&]() {
        auto mpi_comm = mpi::comm().name();
        if (config.has("mpi_comm")) {
            mpi_comm = config.getString("mpi_comm");
        }
        mpi::Scope mpi_scope(mpi_comm);
        auto meshgenerator = MeshGenerator{grid.meshgenerator() | config};
        auto mesh          = meshgenerator.generate(grid, distribution);
        mesh.get()->attach();
        return mesh.get();
    }()) {
    get()->detach();
}



Mesh::Mesh(eckit::Stream& stream): Handle(new Implementation(stream)) {}

Mesh::operator bool() const {
    return get()->nodes().size() > 0;
}

//----------------------------------------------------------------------------------------------------------------------

}  // namespace atlas
