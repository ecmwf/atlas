/*
 * (C) Copyright 1996-2015 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include <cmath>

#include "eckit/exception/Exceptions.h"
#include "eckit/memory/ScopedPtr.h"
#include "atlas/numerics/nabla/EdgeBasedFiniteVolume.h"
#include "atlas/functionspace/NodesFunctionSpace.h"
#include "atlas/Mesh.h"
#include "atlas/Nodes.h"

#include "atlas/actions/BuildEdges.h"
#include "atlas/actions/BuildParallelFields.h"
#include "atlas/actions/BuildDualMesh.h"
#include "atlas/util/ArrayView.h"
#include "atlas/util/IndexView.h"
#include "atlas/atlas_omp.h"
#include "atlas/Parameters.h"
#include "atlas/util/Bitflags.h"



// =======================================================

using namespace atlas::actions;

namespace atlas {
namespace functionspace {

EdgeBasedFiniteVolumeFunctionSpace::EdgeBasedFiniteVolumeFunctionSpace(Mesh &mesh, const Halo &halo )
 : next::FunctionSpace("EdgeBasedFiniteVolume"),
   mesh_(mesh),
   halo_(halo)
{
  nodes_.reset( new NodesFunctionSpace(mesh,halo) );
  if(mesh_.has_function_space("edges")) {
    edges_ = &mesh_.function_space("edges");
  } else {
    build_edges(mesh_);
    build_pole_edges(mesh_);
    edges_ = &mesh_.function_space("edges");
    build_edges_parallel_fields(*edges_,mesh_.nodes());
    build_median_dual_mesh(mesh_);
    build_node_to_edge_connectivity(mesh_);

    size_t nnodes = mesh_.nodes().size();
    IndexView<int,2> node2edge( mesh_.nodes().field("to_edge") );
    ArrayView<int,1> node2edge_size( mesh_.nodes().field("to_edge_size") );
    mesh_.nodes().add( Field::create<double>("node2edge_sign",make_shape(nnodes,node2edge.shape(1)) ) );
    ArrayView<double,2> node2edge_sign( mesh_.nodes().field("node2edge_sign") );
    IndexView<int,   2> edge2node( edges_->field("nodes") );

    atlas_omp_parallel_for( int jnode=0; jnode<nnodes; ++jnode )
    {
      for(size_t jedge = 0; jedge < node2edge_size(jnode); ++jedge)
      {
        size_t iedge = node2edge(jnode,jedge);
        size_t ip1 = edge2node(iedge,0);
        if( jnode == ip1 )
          node2edge_sign(jnode,jedge) = 1.;
        else
          node2edge_sign(jnode,jedge) = -1.;
      }
    }
  }
}


}
}

using namespace atlas::functionspace;
using atlas::util::Topology;

namespace atlas {
namespace numerics {
namespace nabla {

namespace {
static NablaBuilder< EdgeBasedFiniteVolume > __edgebasedfinitevolume("EdgeBasedFiniteVolume");
}

EdgeBasedFiniteVolume::EdgeBasedFiniteVolume(const next::FunctionSpace &fs, const eckit::Parametrisation &p) :
  Nabla(fs,p)
{
  fvm_ = dynamic_cast<const EdgeBasedFiniteVolumeFunctionSpace *>(&fs);
  if( ! fvm_ )
    throw eckit::BadCast("nabla::EdgeBasedFiniteVolume needs a EdgeBasedFiniteVolumeFunctionSpace",Here());
  eckit::Log::info() << "EdgeBasedFiniteVolume constructed for functionspace " << fvm_->name()
                     << " with " << fvm_->nodes_->nb_nodes_global() << " nodes total" << std::endl;

  setup();

}

EdgeBasedFiniteVolume::~EdgeBasedFiniteVolume()
{
}

void EdgeBasedFiniteVolume::setup()
{
  atlas::FunctionSpace &edges = fvm_->nodes_->mesh().function_space("edges");
  atlas::Nodes const   &nodes = fvm_->nodes_->nodes();

  IndexView<int,2> edge2node;

  ArrayView<double,2> lonlat_deg;
  ArrayView<double,1> V;
  ArrayView<double,2> S;

  IndexView<int,   2> node2edge;
  ArrayView<int,   1> node2edge_size;
  ArrayView<double,2> node2edge_sign;
  ArrayView<int,   1> edge_is_pole;

  size_t nnodes;
  size_t nedges;

  nnodes = nodes.size();
  nedges = edges.shape(0);

  edge2node  = IndexView<int,   2> ( edges.field("nodes") );
  lonlat_deg = ArrayView<double,2> ( nodes.lonlat() );
  V      = ArrayView<double,1> ( nodes.field("dual_volumes") );
  S      = ArrayView<double,2> ( edges.field("dual_normals") );

  double radius = Earth::radiusInMeters();
  double deg2rad = M_PI/180.;
  atlas_omp_parallel_for( size_t jnode=0; jnode<nnodes; ++jnode )
  {
    double y  = lonlat_deg(jnode,LAT) * deg2rad;
    double hx = radius*std::cos(y);
    double hy = radius;
    double G  = hx*hy;
    V(jnode) *= std::pow(deg2rad,2) * G;
  }
  atlas_omp_parallel_for( size_t jedge=0; jedge<nedges; ++jedge )
  {
    S(jedge,LON) *= deg2rad;
    S(jedge,LAT) *= deg2rad;
  }

  edge_is_pole   = ArrayView<int,1> ( edges.field("is_pole_edge") );
  node2edge      = IndexView<int,2> ( nodes.field("to_edge") );
  node2edge_size = ArrayView<int,1> ( nodes.field("to_edge_size") );
  node2edge_sign = ArrayView<double,2> ( nodes.field("node2edge_sign") );

  atlas_omp_parallel_for( int jnode=0; jnode<nnodes; ++jnode )
  {
    for(size_t jedge = 0; jedge < node2edge_size(jnode); ++jedge)
    {
      size_t iedge = node2edge(jnode,jedge);
      size_t ip1 = edge2node(iedge,0);
      if( jnode == ip1 )
        node2edge_sign(jnode,jedge) = 1.;
      else
        node2edge_sign(jnode,jedge) = -1.;
    }
  }

  // Filter pole_edges out of all edges
  std::vector<size_t> tmp(nedges);
  size_t c(0);
  for(size_t jedge = 0; jedge < nedges; ++jedge)
  {
    if( edge_is_pole(jedge) )
      tmp[c++] = jedge;
  }
  pole_edges_.clear();
  pole_edges_.reserve(c);
  for( size_t jedge=0; jedge<c; ++jedge )
    pole_edges_.push_back(tmp[jedge]);
}


void EdgeBasedFiniteVolume::gradient(const Field& _field, Field& _grad)
{
  atlas::FunctionSpace &edges = fvm_->nodes_->mesh().function_space("edges");
  atlas::Nodes const   &nodes = fvm_->nodes_->nodes();

  IndexView<int,2> edge2node;

  ArrayView<double,2> lonlat_deg;
  ArrayView<double,1> V;
  ArrayView<double,2> S;

  IndexView<int,   2> node2edge;
  ArrayView<int,   1> node2edge_size;
  ArrayView<double,2> node2edge_sign;
  ArrayView<int,   1> edge_is_pole;

  size_t nnodes;
  size_t nedges;
  size_t nlev=_field.levels();
  if( _grad.levels() != nlev )
    throw eckit::AssertionFailed("gradient field should have same number of levels",Here());

  nnodes = nodes.size();
  nedges = edges.shape(0);

  ArrayView<double,2> field( _field.data<double>(), make_shape(nnodes,nlev) );
  ArrayView<double,3> grad( _grad.data<double>(), make_shape(nnodes,nlev,2) );

  edge2node  = IndexView<int,   2> ( edges.field("nodes") );
  lonlat_deg = ArrayView<double,2> ( nodes.lonlat() );
  V      = ArrayView<double,1> ( nodes.field("dual_volumes") );
  S      = ArrayView<double,2> ( edges.field("dual_normals") );
  edge_is_pole   = ArrayView<int,1> ( edges.field("is_pole_edge") );
  node2edge      = IndexView<int,2> ( nodes.field("to_edge") );
  node2edge_size = ArrayView<int,1> ( nodes.field("to_edge_size") );
  node2edge_sign = ArrayView<double,2> ( nodes.field("node2edge_sign") );

  eckit::SharedPtr<Array> avgS_arr( Array::create<double>(nedges,nlev,2) );
  ArrayView<double,3> avgS(*avgS_arr);

  atlas_omp_parallel_for( int jedge=0; jedge<nedges; ++jedge )
  {
    int ip1 = edge2node(jedge,0);
    int ip2 = edge2node(jedge,1);

    for(size_t jlev = 0; jlev < nlev; ++jlev)
    {
      double avg = ( field(ip1,jlev) + field(ip2,jlev) ) * 0.5;
      avgS(jedge,jlev,LON) = S(jedge,LON)*avg;
      avgS(jedge,jlev,LAT) = S(jedge,LAT)*avg;
    }
  }

  atlas_omp_parallel_for( int jnode=0; jnode<nnodes; ++jnode )
  {
    for(size_t jlev = 0; jlev < nlev; ++jlev )
    {
      grad(jnode,jlev,LON) = 0.;
      grad(jnode,jlev,LAT) = 0.;
    }
    for( int jedge=0; jedge<node2edge_size(jnode); ++jedge )
    {
      int iedge = node2edge(jnode,jedge);
      double add = node2edge_sign(jnode,jedge);
      for(size_t jlev = 0; jlev < nlev; ++jlev)
      {
        grad(jnode,jlev,LON) += add*avgS(iedge,jlev,LON);
        grad(jnode,jlev,LAT) += add*avgS(iedge,jlev,LAT);
      }
    }
    for(size_t jlev = 0; jlev < nlev; ++jlev)
    {
      grad(jnode,jlev,LON) /= V(jnode);
      grad(jnode,jlev,LAT) /= V(jnode);
    }
  }
  // special treatment for the north & south pole cell faces
  // Sx == 0 at pole, and Sy has same sign at both sides of pole
  for(size_t jedge = 0; jedge < pole_edges_.size(); ++jedge)
  {
    int iedge = pole_edges_[jedge];
    int ip2 = edge2node(iedge,1);
    // correct for wrong Y-derivatives in previous loop
    for(size_t jlev = 0; jlev < nlev; ++jlev)
      grad(ip2,jlev,LAT) += 2.*avgS(iedge,jlev,LAT)/V(ip2);
  }

  // halo-exchange
  fvm_->nodes_->haloExchange(_grad);

}



} // namespace nabla
} // namespace numerics
} // namespace atlas
