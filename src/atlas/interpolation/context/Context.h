/*
 * (C) Copyright 1996-2017 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */


#ifndef atlas_interpolation_context_Context_h
#define atlas_interpolation_context_Context_h

//include "eckit/linalg/LinearAlgebra.h"
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

//#include "atlas/grid/partitioners/PartitionerFromPrePartitionedMesh.h"
//include "atlas/interpolation/Interpolation.h"


namespace atlas {
namespace interpolation {
namespace context {


struct Context {

    Context(
            const std::string& gridname,
            const std::string& partitioner,
            const std::string& meshGenerator,
            bool meshGeneratorTriangulate = false,
            double meshGeneratorAngle = 0 );

    virtual ~Context() {}

    size_t meshHaloSize(size_t haloSize) { meshHaloSize_ = haloSize; }
    size_t meshHaloSize() const { return meshHaloSize_; }

    const grid::Structured& grid()  const { return *grid_; }
    const mesh::Mesh& mesh()        const { return *mesh_; }
    const mesh::Mesh::Ptr meshPtr() const { return  mesh_; }
    const grid::partitioners::Partitioner& partitioner() const { return *partitioner_; }
    const field::FieldSet& fieldSet()                    const { return *fieldSet_; }
    const functionspace::NodeColumns& functionSpace()    const { return *functionSpace_; }

    grid::Structured& grid()  { return *grid_; }
    mesh::Mesh& mesh()        { return *mesh_; }
    mesh::Mesh::Ptr meshPtr() { return  mesh_; }
    grid::partitioners::Partitioner& partitioner() { return *partitioner_; }
    field::FieldSet& fieldSet()                    { return *fieldSet_; }
    functionspace::NodeColumns& functionSpace()    { return *functionSpace_; }

protected:

    const std::string& optionGridname_;
    const std::string& optionPartitioner_;
    const std::string& optionMeshGenerator_;

    size_t meshHaloSize_;
    mesh::generators::MeshGenerator::Parameters meshGeneratorParams_;

    grid::Structured::Ptr grid_;
    grid::partitioners::Partitioner::Ptr partitioner_;
    mesh::Mesh::Ptr mesh_;
    grid::GridDistribution::Ptr distribution_;

    functionspace::NodeColumns::Ptr functionSpace_;
    field::FieldSet::Ptr fieldSet_;

};


}  // context
}  // interpolation
}  // atlas


#endif
