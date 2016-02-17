/*
 * (C) Copyright 1996-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include <cmath>

#include "eckit/exception/Exceptions.h"
#include "atlas/Mesh.h"
#include "atlas/FunctionSpace.h"
#include "atlas/mesh/Nodes.h"
#include "atlas/mesh/HybridElements.h"
#include "atlas/runtime/ErrorHandling.h"

#include "atlas/functionspace/EdgeBasedFiniteVolume.h"
#include "atlas/actions/BuildEdges.h"
#include "atlas/actions/BuildParallelFields.h"
#include "atlas/actions/BuildDualMesh.h"
#include "atlas/util/ArrayView.h"
#include "atlas/util/IndexView.h"
#include "atlas/util/Bitflags.h"
#include "atlas/atlas_omp.h"
#include "atlas/Parameters.h"

// =======================================================

using namespace atlas::actions;

namespace atlas {
namespace functionspace {

namespace {
  Halo EdgeBasedFiniteVolume_halo(const eckit::Parametrisation &params)
  {
    size_t halo_size(1);
    params.get("halo",halo_size);
    return Halo(halo_size);
  }
}

EdgeBasedFiniteVolume::EdgeBasedFiniteVolume( Mesh &mesh )
 : Nodes(mesh,Halo(1)), edges_(mesh.edges())
{
  radius_ = Earth::radiusInMeters();
  setup();
}

EdgeBasedFiniteVolume::EdgeBasedFiniteVolume( Mesh &mesh, const Halo &halo )
 : Nodes(mesh,halo), edges_(mesh.edges())
{
  radius_ = Earth::radiusInMeters();
  setup();
}

EdgeBasedFiniteVolume::EdgeBasedFiniteVolume( Mesh &mesh, const eckit::Parametrisation &params )
  : Nodes(mesh,EdgeBasedFiniteVolume_halo(params),params), edges_(mesh.edges())
{
  radius_ = Earth::radiusInMeters();
  params.get("radius",radius_);
  setup();
}

void EdgeBasedFiniteVolume::setup()
{
  if( edges_.size() == 0 )
  {
    build_edges(mesh());
    build_pole_edges(mesh());
    build_edges_parallel_fields( mesh() );
    build_median_dual_mesh(mesh());
    build_node_to_edge_connectivity(mesh());

    const size_t nnodes = nodes().size();

    // Compute sign
    {
      const ArrayView<int,1> is_pole_edge   ( edges_.field("is_pole_edge") );

      const Connectivity &node_edge_connectivity = nodes().edge_connectivity();
      const Connectivity &edge_node_connectivity = edges_.node_connectivity();
      nodes().add( Field::create<double>("node2edge_sign",make_shape(nnodes,node_edge_connectivity.maxcols()) ) );
      ArrayView<double,2> node2edge_sign( nodes().field("node2edge_sign") );

      atlas_omp_parallel_for( int jnode=0; jnode<nnodes; ++jnode )
      {
        for(size_t jedge = 0; jedge < node_edge_connectivity.cols(jnode); ++jedge)
        {
          size_t iedge = node_edge_connectivity(jnode,jedge);
          size_t ip1 = edge_node_connectivity(iedge,0);
          if( jnode == ip1 )
            node2edge_sign(jnode,jedge) = 1.;
          else
          {
            node2edge_sign(jnode,jedge) = -1.;
            if( is_pole_edge(iedge) )
              node2edge_sign(jnode,jedge) = 1.;
          }
        }
      }
    }

    // Metrics
    {
      const size_t nedges = edges_.size();
      const ArrayView<double,2> lonlat_deg( nodes().lonlat() );
      ArrayView<double,1> dual_volumes ( nodes().field("dual_volumes") );
      ArrayView<double,2> dual_normals ( edges_.field("dual_normals") );

      const double deg2rad = M_PI/180.;
      atlas_omp_parallel_for( size_t jnode=0; jnode<nnodes; ++jnode )
      {
        double y  = lonlat_deg(jnode,LAT) * deg2rad;
        double hx = radius_*std::cos(y);
        double hy = radius_;
        double G  = hx*hy;
        dual_volumes(jnode) *= std::pow(deg2rad,2) * G;
      }

      atlas_omp_parallel_for( size_t jedge=0; jedge<nedges; ++jedge )
      {
        dual_normals(jedge,LON) *= deg2rad;
        dual_normals(jedge,LAT) *= deg2rad;
      }
    }
  }
  edges_fs_.reset( new functionspace::Edges(mesh()) );
}

// ------------------------------------------------------------------------------------------
extern "C" {

EdgeBasedFiniteVolume* atlas__functionspace__EdgeBasedFiniteVolume__new (Mesh* mesh, const eckit::Parametrisation* params)
{
  EdgeBasedFiniteVolume* functionspace(0);
  ATLAS_ERROR_HANDLING(
    ASSERT(mesh);
    functionspace = new EdgeBasedFiniteVolume(*mesh,*params);
  );
  return functionspace;
}

}
// ------------------------------------------------------------------------------------------

} // namespace functionspace
} // namespace atlas
