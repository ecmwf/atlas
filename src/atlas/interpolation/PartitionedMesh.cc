/*
 * (C) Copyright 1996-2017 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */


#include "atlas/interpolation/PartitionedMesh.h"

#include <typeinfo>
#include "eckit/log/Timer.h"
#include "atlas/functionspace/NodeColumns.h"
#include "atlas/grid/GridDistribution.h"
#include "atlas/grid/partitioners/PartitionerFromPrePartitionedMesh.h"
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
    output_config.set("coordinates", std::string("xyz"));

    output::Gmsh out(fileName, output_config);
    out.write(*mesh_);

    if (fields) {
        out.write(*fields);
    }
}


void PartitionedMesh::partition(const grid::Grid& grid) {
    eckit::TraceTimer<LibAtlas> tim("PartitionedMesh::partition()");

    partitioner_.reset(grid::partitioners::PartitionerFactory::build(optionPartitioner_, grid));
    ASSERT(partitioner_);

    grid::GridDistribution::Ptr dist;
    dist.reset(partitioner_->distribution());
    ASSERT(dist);

    Generator::Ptr meshgen(mesh::generators::MeshGeneratorFactory::build(optionGenerator_, generatorParams_));
    mesh_.reset(meshgen->generate(grid, *dist));
    ASSERT(mesh_);
}


void PartitionedMesh::partition(const grid::Grid& grid, const PartitionedMesh& other) {
    eckit::TraceTimer<LibAtlas> tim("PartitionedMesh::partition(other)");

    ASSERT(other.mesh_);
    ASSERT(other.partitioner_);

    Partitioner::Ptr partitioner_(grid::partitioners::PartitionerFactory::build(optionPartitioner_, grid));
    ASSERT(partitioner_);

    grid::GridDistribution::Ptr dist;
    try {
        typedef grid::partitioners::PartitionerFromPrePartitionedMesh PrePartitioner;
        PrePartitioner& partner = dynamic_cast< PrePartitioner& >(*partitioner_);

        partner.setup(other.mesh_, other.partitioner_->grid().domain());
        dist.reset(partner.distribution());

    } catch (std::bad_cast&) {
        throw eckit::UserError("Partitioner has to be a PartitionerFromPrePartitionedMesh");
    }
    ASSERT(dist);

    Generator::Ptr meshgen(mesh::generators::MeshGeneratorFactory::build(optionGenerator_, generatorParams_));
    mesh_.reset(meshgen->generate(grid, *dist));
    ASSERT(mesh_);
}


}  // namespace interpolation
}  // namespace atlas

