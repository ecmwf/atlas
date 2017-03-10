/*
 * (C) Copyright 1996-2017 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */


#include "PartitionedMesh.h"

#include <typeinfo>
#include "eckit/log/Timer.h"
#include "atlas/functionspace/NodeColumns.h"
#include "atlas/grid/Distribution.h"
#include "atlas/grid/Partitioner.h"
#include "atlas/grid/detail/partitioner/Partitioner.h"
#include "atlas/output/Gmsh.h"
#include "atlas/runtime/LibAtlas.h"
#include "atlas/runtime/Log.h"


namespace atlas {
namespace interpolation {

PartitionedMesh::PartitionedMesh(
        const std::string& partitioner,
        const std::string& generator,
        bool generatorTriangulate,
        double generatorAngle ) :
    optionPartitioner_(partitioner),
    optionGenerator_(generator) {

    generatorParams_.set("three_dimensional", false );
    generatorParams_.set("patch_pole",        true  );
    generatorParams_.set("include_pole",      false );
    generatorParams_.set("triangulate",        generatorTriangulate );
    generatorParams_.set("angle",              generatorAngle );
}


void PartitionedMesh::writeGmsh(const std::string& fileName, const field::FieldSet* fields) {
    ASSERT(mesh_);

    util::Config output_config;
    // output_config.set("coordinates", std::string("xyz"));
    output_config.set("ghost", true);

    output::Gmsh out(fileName, output_config);
    out.write(*mesh_);

    if (fields) {
        out.write(*fields);
    }
}


void PartitionedMesh::partition(const grid::Grid& grid) {
    eckit::TraceTimer<LibAtlas> tim("PartitionedMesh::partition()");

    partitioner_ = Partitioner(optionPartitioner_, grid);


    Generator::Ptr meshgen(meshgenerator::MeshGeneratorFactory::build(optionGenerator_, generatorParams_));
    mesh_.reset(meshgen->generate(grid, partitioner_.distribution()));
}


void PartitionedMesh::partition(const grid::Grid& grid, const PartitionedMesh& other) {
    eckit::TraceTimer<LibAtlas> tim("PartitionedMesh::partition(other)");

    partitioner_ = Partitioner( grid::MatchedPartitionerFactory::build(optionPartitioner_,grid, *other.mesh_));

    Generator::Ptr meshgen(meshgenerator::MeshGeneratorFactory::build(optionGenerator_, generatorParams_));
    mesh_.reset(meshgen->generate(grid, partitioner_.distribution() ));
}


}  // namespace interpolation
}  // namespace atlas

