/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "PartitionedMesh.h"

#include "atlas/grid/Partitioner.h"
#include "atlas/output/Gmsh.h"
#include "atlas/runtime/Log.h"
#include "atlas/runtime/Trace.h"

namespace atlas {
namespace interpolation {

PartitionedMesh::PartitionedMesh(const std::string& partitioner, const std::string& generator,
                                 bool generatorTriangulate, double generatorAngle, bool patchPole):
    optionPartitioner_(partitioner), optionGenerator_(generator) {
    generatorParams_.set("three_dimensional", false);
    generatorParams_.set("patch_pole", patchPole);
    generatorParams_.set("include_pole", false);
    generatorParams_.set("triangulate", generatorTriangulate);
    generatorParams_.set("angle", generatorAngle);
    generatorParams_.set("fixup", true);
}

void PartitionedMesh::writeGmsh(const std::string& fileName, const FieldSet& fields) {
    util::Config output_config;
    //output_config.set( "coordinates", std::string( "xyz" ) );
    output_config.set("ghost", true);

    output::Gmsh out(fileName, output_config);
    out.write(mesh_);

    if (not fields.empty()) {
        out.write(fields);
    }
}

void PartitionedMesh::partition(const Grid& grid) {
    ATLAS_TRACE("PartitionedMesh::partition()");

    auto meshgen_config = grid.meshgenerator();
    meshgen_config.set(generatorParams_);
    if (optionGenerator_ != "default") {
        meshgen_config.set("type", optionGenerator_);
    }
    MeshGenerator meshgen(meshgen_config);

    auto partitioner_config = grid.partitioner();
    if (optionPartitioner_ != "default") {
        partitioner_config.set("type", optionPartitioner_);
    }
    partitioner_ = Partitioner(partitioner_config);

    mesh_ = meshgen.generate(grid, partitioner_.partition(grid));
}

void PartitionedMesh::partition(const Grid& grid, const PartitionedMesh& other) {
    ATLAS_TRACE("PartitionedMesh::partition(other)");

    partitioner_ = grid::MatchingMeshPartitioner(other.mesh_, util::Config("type", optionPartitioner_));

    auto meshgen_config = grid.meshgenerator();
    meshgen_config.set(generatorParams_);
    if (optionGenerator_ != "default") {
        meshgen_config.set("type", optionGenerator_);
    }
    MeshGenerator meshgen(meshgen_config);

    mesh_ = meshgen.generate(grid, partitioner_.partition(grid));
}

}  // namespace interpolation
}  // namespace atlas
