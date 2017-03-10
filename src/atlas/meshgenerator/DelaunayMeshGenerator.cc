/*
 * (C) Copyright 1996-2017 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "eckit/utils/MD5.h"
#include "atlas/meshgenerator/DelaunayMeshGenerator.h"
#include "atlas/grid/Distribution.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/mesh/Nodes.h"
#include "atlas/field/Field.h"
#include "atlas/mesh/actions/ExtendNodesGlobal.h"
#include "atlas/mesh/actions/BuildXYZField.h"
#include "atlas/mesh/actions/BuildConvexHull3D.h"
#include "atlas/mesh/HybridElements.h"
#include "atlas/array/ArrayView.h"
#include "atlas/runtime/Log.h"

using atlas::mesh::Mesh;

namespace atlas {
namespace meshgenerator {

//----------------------------------------------------------------------------------------------------------------------

DelaunayMeshGenerator::DelaunayMeshGenerator()
{
}

DelaunayMeshGenerator::DelaunayMeshGenerator(const eckit::Parametrisation& p)
{
}


DelaunayMeshGenerator::~DelaunayMeshGenerator() {
}

void DelaunayMeshGenerator::hash(eckit::MD5& md5) const
{
    md5.add("Delaunay");

    // no other settings
}

void DelaunayMeshGenerator::generate(const grid::Grid& grid, const grid::Distribution& dist, mesh::Mesh& mesh) const
{
  if( dist.nb_partitions() > 1 )
  {
    Log::warning() << "Delaunay triangulation does not support a GridDistribution"
                             "with more than 1 partition"
                          << std::endl;
    NOTIMP;
    /// TODO: Read mesh on 1 MPI task, and distribute according to GridDistribution
    /// HINT: use atlas/actions/DistributeMesh
  }
  else
  {
    generate(grid, mesh);
  }
}

void DelaunayMeshGenerator::generate(const grid::Grid& g, mesh::Mesh& mesh) const
{
  mesh.createNodes(g);

  array::ArrayView<gidx_t,1> gidx( mesh.nodes().global_index() );
  for( size_t jnode=0; jnode<mesh.nodes().size(); ++ jnode ) {
    gidx(jnode) = jnode+1;
  }

  mesh::actions::BuildXYZField()(mesh);
  mesh::actions::ExtendNodesGlobal()(g, mesh);    ///< does nothing if global domain
  mesh::actions::BuildConvexHull3D()(mesh);
}

namespace {
static MeshGeneratorBuilder< DelaunayMeshGenerator > __delaunay("delaunay");
}

//----------------------------------------------------------------------------------------------------------------------

} // namespace meshgenerator
} // namespace atlas

