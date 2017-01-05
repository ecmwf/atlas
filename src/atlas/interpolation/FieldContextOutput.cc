/*
 * (C) Copyright 1996-2017 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */


#include "atlas/interpolation/FieldContextOutput.h"

//include "eckit/linalg/Vector.h"
//include "atlas/atlas.h"
#include "atlas/functionspace/NodeColumns.h"
#include "atlas/grid/GridDistribution.h"
#include "atlas/grid/partitioners/Partitioner.h"
#include "atlas/grid/Structured.h"
//include "atlas/internals/AtlasTool.h"
#include "atlas/mesh/generators/Structured.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/mesh/Nodes.h"
#include "atlas/output/Gmsh.h"
//include "atlas/runtime/Log.h"

#include "atlas/grid/partitioners/PartitionerFromPrePartitionedMesh.h"
//include "atlas/interpolation/Interpolation.h"


namespace atlas {
namespace interpolation {


FieldContextOutput::FieldContextOutput(
        const std::string& gridname,
        const std::string& partitioner,
        const std::string& meshGenerator,
        const mesh::Mesh::Ptr prePartitionedMesh,
        const grid::Domain& prePartitionedDomain,
        bool meshGeneratorTriangulate,
        double meshGeneratorAngle ) :
    FieldContext(gridname, partitioner, meshGenerator, meshGeneratorTriangulate, meshGeneratorAngle) {
    using grid::partitioners::Partitioner;
    using mesh::generators::MeshGenerator;

    grid_.reset(grid::Structured::create(optionGridname_));
    ASSERT(grid_);

    grid::partitioners::Partitioner::Ptr partitioner_(grid::partitioners::PartitionerFactory::build(optionPartitioner_, *grid_));
    ASSERT(partitioner_);

    try {
        grid::partitioners::PartitionerFromPrePartitionedMesh& partner = dynamic_cast< grid::partitioners::PartitionerFromPrePartitionedMesh& >(*partitioner_);
        partner.setup(prePartitionedMesh, prePartitionedDomain);
        distribution_.reset(partner.distribution());
    } catch (std::bad_cast&) {
        throw eckit::UserError("Partitioner has to be a PartitionerFromPrePartitionedMesh");
    }
    ASSERT(distribution_);

    MeshGenerator::Ptr meshgen(mesh::generators::MeshGeneratorFactory::build(optionMeshGenerator_, meshGeneratorParams_));
    mesh_.reset(meshgen->generate(*grid_, *distribution_));
}


void FieldContextOutput::write(const std::string& fileName) {
    util::Config output_config;
    output_config.set("coordinates", std::string("xyz"));
    output::Gmsh out(fileName, output_config);
    out.write(*mesh_);
    if (fieldSet_) {
        for (size_t i = 0; i < fieldSet_->size(); ++i) {
            out.write(fieldSet_->operator[](i));
        }
    }
}


}  // interpolation
}  // atlas
