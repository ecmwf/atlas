/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#pragma once

#include "atlas/field.h"
#include "atlas/grid.h"
#include "atlas/mesh.h"
#include "atlas/meshgenerator.h"

namespace atlas {
namespace interpolation {

struct PartitionedMesh {
    typedef grid::Partitioner Partitioner;

    PartitionedMesh(const std::string& partitioner, const std::string& generator, bool meshGeneratorTriangulate = false,
                    double meshGeneratorAngle = 0, bool patchPole = true);

    virtual ~PartitionedMesh() {}

    const Mesh& mesh() const { return mesh_; }
    Mesh& mesh() { return mesh_; }

    void writeGmsh(const std::string& fileName, const FieldSet& fields = FieldSet());

    void partition(const Grid&);
    void partition(const Grid&, const PartitionedMesh&);

protected:
    const std::string optionPartitioner_;
    const std::string optionGenerator_;

    MeshGenerator::Parameters generatorParams_;
    Partitioner partitioner_;
    Mesh mesh_;
};

}  // namespace interpolation
}  // namespace atlas
