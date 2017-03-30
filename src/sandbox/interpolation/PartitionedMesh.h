/*
 * (C) Copyright 1996-2017 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#pragma once

#include "atlas/grid/Partitioner.h"
#include "atlas/grid/Grid.h"
#include "atlas/meshgenerator/StructuredMeshGenerator.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/mesh/Nodes.h"


namespace atlas {
namespace field { class FieldSet; }
}


namespace atlas {
namespace interpolation {


struct PartitionedMesh {

    typedef grid::Partitioner Partitioner;
    typedef meshgenerator::MeshGenerator Generator;
    typedef mesh::Mesh                      Mesh;

    PartitionedMesh(
            const std::string& partitioner,
            const std::string& generator,
            bool meshGeneratorTriangulate = false,
            double meshGeneratorAngle = 0 );

    virtual ~PartitionedMesh() {}

    const Mesh& mesh() const { return mesh_; }
    Mesh&       mesh()       { return mesh_; }

    void writeGmsh(const std::string& fileName, const field::FieldSet* fields = NULL);

    void partition(const grid::Grid&);
    void partition(const grid::Grid&, const PartitionedMesh&);

protected:

    const std::string optionPartitioner_;
    const std::string optionGenerator_;

    Generator::Parameters generatorParams_;
    Partitioner partitioner_;
    Mesh mesh_;

};


}  // namespace interpolation
}  // namespace atlas
