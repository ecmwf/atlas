/*
 * (C) Copyright 1996-2017 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include <cmath>
#include "eckit/exception/Exceptions.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/mesh/Nodes.h"
#include "atlas/mesh/HybridElements.h"
#include "atlas/mesh/actions/BuildEdges.h"
#include "atlas/mesh/actions/BuildParallelFields.h"
#include "atlas/mesh/actions/BuildDualMesh.h"
#include "atlas/functionspace/NodeColumns.h"
#include "atlas/functionspace/EdgeColumns.h"
#include "atlas/numerics/fvm/Method.h"
#include "atlas/internals/Parameters.h"
#include "atlas/util/Constants.h"
#include "atlas/internals/Bitflags.h"
#include "atlas/array/ArrayView.h"
#include "atlas/array/MakeView.h"
#include "atlas/parallel/omp/omp.h"
#include "atlas/runtime/ErrorHandling.h"

// =======================================================

using namespace atlas::mesh::actions;

namespace atlas {
namespace numerics {
namespace fvm {

namespace {
  mesh::Halo Method_halo(const eckit::Parametrisation &params)
  {
    size_t halo_size(1);
    params.get("halo",halo_size);
    return mesh::Halo(halo_size);
  }
}

Method::Method( mesh::Mesh &mesh ) :
  mesh_(mesh),
  halo_(mesh),
  nodes_(mesh.nodes()),
  edges_(mesh.edges()),
  radius_(util::Earth::radiusInMeters())
{
  setup();
}

Method::Method( mesh::Mesh &mesh, const mesh::Halo &halo ) :
  mesh_(mesh),
  halo_(halo),
  nodes_(mesh.nodes()),
  edges_(mesh.edges()),
  radius_(util::Earth::radiusInMeters())
{
  setup();
}

Method::Method( mesh::Mesh &mesh, const eckit::Parametrisation &params ) :
  mesh_(mesh),
  halo_(Method_halo(params)),
  nodes_(mesh.nodes()),
  edges_(mesh.edges()),
  radius_(util::Earth::radiusInMeters())
{
  params.get("radius",radius_);
  setup();
}

void Method::setup()
{
  node_columns_.reset( new functionspace::NodeColumns(mesh(),halo_) );
  if( edges_.size() == 0 )
  {
    build_edges(mesh());
    build_pole_edges(mesh());
    build_edges_parallel_fields( mesh() );
    build_median_dual_mesh(mesh());
    build_node_to_edge_connectivity(mesh());

    const size_t nnodes = nodes_.size();

    // Compute sign
    {
      const array::ArrayView<int,1> is_pole_edge = array::make_view<int,1>( edges_.field("is_pole_edge") );

      const mesh::Connectivity &node_edge_connectivity = nodes_.edge_connectivity();
      const mesh::MultiBlockConnectivity &edge_node_connectivity = edges_.node_connectivity();
      if( ! nodes_.has_field("node2edge_sign") )
      {
        nodes_.add(
              field::Field::create<double>("node2edge_sign",
              array::make_shape(nnodes,node_edge_connectivity.maxcols()) ) );
      }
      array::ArrayView<double,2> node2edge_sign = array::make_view<double,2>( nodes_.field("node2edge_sign") );

      atlas_omp_parallel_for( size_t jnode=0; jnode<nnodes; ++jnode )
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
    if (0) {
      const size_t nedges = edges_.size();
      const array::ArrayView<double,2> lonlat_deg = array::make_view<double,2>( nodes_.lonlat() );
      array::ArrayView<double,1> dual_volumes = array::make_view<double,1>( nodes_.field("dual_volumes") );
      array::ArrayView<double,2> dual_normals = array::make_view<double,2>( edges_.field("dual_normals") );

      const double deg2rad = M_PI/180.;
      atlas_omp_parallel_for( size_t jnode=0; jnode<nnodes; ++jnode )
      {
        double y  = lonlat_deg(jnode,internals::LAT) * deg2rad;
        double hx = radius_*std::cos(y);
        double hy = radius_;
        double G  = hx*hy;
        dual_volumes(jnode) *= std::pow(deg2rad,2) * G;
      }

      atlas_omp_parallel_for( size_t jedge=0; jedge<nedges; ++jedge )
      {
        dual_normals(jedge,internals::LON) *= deg2rad;
        dual_normals(jedge,internals::LAT) *= deg2rad;
      }
    }
  }
  edge_columns_.reset( new functionspace::EdgeColumns(mesh()) );
}

// ------------------------------------------------------------------------------------------
extern "C" {
Method* atlas__numerics__fvm__Method__new (mesh::Mesh* mesh, const eckit::Parametrisation* params)
{
  Method* method(0);
  ATLAS_ERROR_HANDLING(
    ASSERT(mesh);
    method = new Method(*mesh,*params);
  );
  return method;
}

functionspace::NodeColumns* atlas__numerics__fvm__Method__functionspace_nodes (Method* This)
{
  ATLAS_ERROR_HANDLING(
        ASSERT(This);
        return &This->node_columns();
  );
  return 0;
}

functionspace::EdgeColumns* atlas__numerics__fvm__Method__functionspace_edges (Method* This)
{
  ATLAS_ERROR_HANDLING(
        ASSERT(This);
        return &This->edge_columns();
  );
  return 0;
}



}
// ------------------------------------------------------------------------------------------

} // namepace fvm
} // namespace numerics
} // namespace atlas
