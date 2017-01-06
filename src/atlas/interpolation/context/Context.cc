/*
 * (C) Copyright 1996-2017 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */


#include "atlas/interpolation/context/Context.h"

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
namespace context {


Context::Context(
        const std::string& gridname,
        const std::string& partitioner,
        const std::string& meshGenerator,
        bool meshGeneratorTriangulate,
        double meshGeneratorAngle ) :
    optionGridname_(gridname),
    optionPartitioner_(partitioner),
    optionMeshGenerator_(meshGenerator),
    meshHaloSize_(0) {

    meshGeneratorParams_.set("three_dimensional", false);
    meshGeneratorParams_.set("patch_pole",        true);
    meshGeneratorParams_.set("include_pole",      false);
    meshGeneratorParams_.set("triangulate", meshGeneratorTriangulate);
    meshGeneratorParams_.set("angle",       meshGeneratorAngle);
}


}  // context
}  // interpolation
}  // atlas
