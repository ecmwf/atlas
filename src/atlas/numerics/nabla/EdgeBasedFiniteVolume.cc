/*
 * (C) Copyright 1996-2015 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "eckit/exception/Exceptions.h"
#include "eckit/config/Parametrisation.h"
#include "atlas/numerics/nabla/EdgeBasedFiniteVolume.h"
#include "atlas/functionspace/EdgeBasedFiniteVolume.h"
#include "atlas/Mesh.h"
#include "atlas/Nodes.h"
#include "atlas/Field.h"

#include "atlas/util/ArrayView.h"
#include "atlas/util/IndexView.h"
#include "atlas/atlas_omp.h"
#include "atlas/Parameters.h"
#include "atlas/util/Bitflags.h"
#include "atlas/io/Gmsh.h"



// =======================================================

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
  fvm_ = dynamic_cast<const functionspace::EdgeBasedFiniteVolume *>(&fs);
  if( ! fvm_ )
    throw eckit::BadCast("nabla::EdgeBasedFiniteVolume needs a EdgeBasedFiniteVolumeFunctionSpace",Here());
  eckit::Log::info() << "EdgeBasedFiniteVolume constructed for functionspace " << fvm_->name()
                     << " with " << fvm_->nb_nodes_global() << " nodes total" << std::endl;

  setup();

}

EdgeBasedFiniteVolume::~EdgeBasedFiniteVolume()
{
}

void EdgeBasedFiniteVolume::setup()
{
  atlas::FunctionSpace &edges = fvm_->mesh().function_space("edges");
  atlas::Nodes const   &nodes = fvm_->nodes();

  const size_t nnodes = nodes.size();
  const size_t nedges = edges.shape(0);

  const ArrayView<int,1> edge_is_pole ( edges.field("is_pole_edge") );

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


void EdgeBasedFiniteVolume::gradient(const Field& scalar_field, Field& grad_field) const
{
  const double radius = fvm_->radius();
  const double deg2rad = M_PI/180.;

  atlas::FunctionSpace &edges = fvm_->mesh().function_space("edges");
  atlas::Nodes const   &nodes = fvm_->nodes();

  const size_t nnodes = nodes.size();
  const size_t nedges = edges.shape(0);
  const size_t nlev = scalar_field.levels();
  if( grad_field.levels() != nlev )
    throw eckit::AssertionFailed("gradient field should have same number of levels",Here());


  const ArrayView<double,2> scalar ( scalar_field.data<double>(), make_shape(nnodes,nlev)   );
        ArrayView<double,3> grad   ( grad_field.  data<double>(), make_shape(nnodes,nlev,2) );

  const ArrayView<double,2> lonlat_deg     ( nodes.lonlat() );
  const IndexView<int,   2> edge2node      ( edges.field("nodes") );
  const ArrayView<double,1> V              ( nodes.field("dual_volumes") );
  const ArrayView<double,2> S              ( edges.field("dual_normals") );
  const ArrayView<int,   1> edge_is_pole   ( edges.field("is_pole_edge") );
  const IndexView<int,   2> node2edge      ( nodes.field("to_edge") );
  const ArrayView<int,   1> node2edge_size ( nodes.field("to_edge_size") );
  const ArrayView<double,2> node2edge_sign ( nodes.field("node2edge_sign") );

  ArrayT<double> avgS_arr( nedges,nlev,2 );
  ArrayView<double,3> avgS(avgS_arr);

  atlas_omp_parallel_for( size_t jedge=0; jedge<nedges; ++jedge )
  {
    int ip1 = edge2node(jedge,0);
    int ip2 = edge2node(jedge,1);

    for(size_t jlev = 0; jlev < nlev; ++jlev)
    {
      double avg = ( scalar(ip1,jlev) + scalar(ip2,jlev) ) * 0.5;
      avgS(jedge,jlev,LON) = S(jedge,LON)*avg;
      avgS(jedge,jlev,LAT) = S(jedge,LAT)*avg;
    }
  }

  atlas_omp_parallel_for( size_t jnode=0; jnode<nnodes; ++jnode )
  {
    for(size_t jlev = 0; jlev < nlev; ++jlev )
    {
      grad(jnode,jlev,LON) = 0.;
      grad(jnode,jlev,LAT) = 0.;
    }
    for( size_t jedge=0; jedge<node2edge_size(jnode); ++jedge )
    {
      int iedge = node2edge(jnode,jedge);
      double add = node2edge_sign(jnode,jedge);
      for(size_t jlev = 0; jlev < nlev; ++jlev)
      {
        grad(jnode,jlev,LON) += add*avgS(iedge,jlev,LON);
        grad(jnode,jlev,LAT) += add*avgS(iedge,jlev,LAT);
      }
    }
    double y  = lonlat_deg(jnode,LAT) * deg2rad;
    double hx = radius*std::cos(y);
    double hy = radius;
    for(size_t jlev = 0; jlev < nlev; ++jlev)
    {
      grad(jnode,jlev,LON) *= hy/V(jnode);
      grad(jnode,jlev,LAT) *= hx/V(jnode);
    }
  }

  // special treatment for the north & south pole cell faces
  // Sx == 0 at pole, and Sy has same sign at both sides of pole
  for(size_t jedge = 0; jedge < pole_edges_.size(); ++jedge)
  {
    int iedge = pole_edges_[jedge];
    int ip2 = edge2node(iedge,1);
    double y  = lonlat_deg(ip2,LAT) * deg2rad;
    double hx = radius*std::cos(y);
    // correct for wrong Y-derivatives in previous loop
    for(size_t jlev = 0; jlev < nlev; ++jlev)
      grad(ip2,jlev,LAT) += 2.*avgS(iedge,jlev,LAT)*hx/V(ip2);
  }
}

// ================================================================================

void EdgeBasedFiniteVolume::divergence(const Field& vector_field, Field& div_field) const
{
  const double radius = fvm_->radius();
  const double deg2rad = M_PI/180.;

  atlas::FunctionSpace &edges = fvm_->mesh().function_space("edges");
  atlas::Nodes const   &nodes = fvm_->nodes();

  const size_t nnodes = nodes.size();
  const size_t nedges = edges.shape(0);
  const size_t nlev = vector_field.levels();
  if( div_field.levels() != nlev )
    throw eckit::AssertionFailed("gradient field should have same number of levels",Here());

  const ArrayView<double,3> vector ( vector_field.data<double>(), make_shape(nnodes,nlev,2));
        ArrayView<double,2> div    ( div_field   .data<double>(), make_shape(nnodes,nlev)  );

  const ArrayView<double,2> lonlat_deg     ( nodes.lonlat() );
  const IndexView<int,   2> edge2node      ( edges.field("nodes") );
  const ArrayView<double,1> V              ( nodes.field("dual_volumes") );
  const ArrayView<double,2> S              ( edges.field("dual_normals") );
  const ArrayView<int,   1> edge_is_pole   ( edges.field("is_pole_edge") );
  const IndexView<int,   2> node2edge      ( nodes.field("to_edge") );
  const ArrayView<int,   1> node2edge_size ( nodes.field("to_edge_size") );
  const ArrayView<double,2> node2edge_sign ( nodes.field("node2edge_sign") );

  ArrayT<double> avgS_arr( nedges,nlev,2 );
  ArrayView<double,3> avgS(avgS_arr);

  atlas_omp_parallel_for( size_t jedge=0; jedge<nedges; ++jedge )
  {
    int ip1 = edge2node(jedge,0);
    int ip2 = edge2node(jedge,1);
    double y1  = lonlat_deg(ip1,LAT) * deg2rad;
    double y2  = lonlat_deg(ip2,LAT) * deg2rad;
    double cosy1 = std::cos(y1);
    double cosy2 = std::cos(y2);

    for(size_t jlev = 0; jlev < nlev; ++jlev)
    {
      double avg[2] = {
        (vector(ip1,jlev,LON) + vector(ip2,jlev,LON) ) * 0.5,
        (cosy1*vector(ip1,jlev,LAT) + cosy2*vector(ip2,jlev,LAT) ) * 0.5
      };
      avgS(jedge,jlev,LON) = S(jedge,LON)*avg[LON];
      avgS(jedge,jlev,LAT) = S(jedge,LAT)*avg[LAT];
      // We don't need the cross terms for divergence
    }
  }


  ArrayT<double> dudx_arr( nedges,nlev );
  ArrayT<double> dvdy_arr( nedges,nlev );
  ArrayView<double,2> dudx(dudx_arr);
  ArrayView<double,2> dvdy(dvdy_arr);

//  eckit::SharedPtr<Field> dudx_field ( fvm_->createField<double>("dudx",nlev) );
//  eckit::SharedPtr<Field> dvdy_field ( fvm_->createField<double>("dvdy",nlev) );
//  ArrayView<double,2> dudx(*dudx_field);
//  ArrayView<double,2> dvdy(*dvdy_field);

  atlas_omp_parallel_for( size_t jnode=0; jnode<nnodes; ++jnode )
  {
    for(size_t jlev = 0; jlev < nlev; ++jlev )
    {
      dudx(jnode,jlev) = 0.;
      dvdy(jnode,jlev) = 0.;
    }
    for( size_t jedge=0; jedge<node2edge_size(jnode); ++jedge )
    {
      int iedge = node2edge(jnode,jedge);
      double add = node2edge_sign(jnode,jedge);
      for(size_t jlev = 0; jlev < nlev; ++jlev)
      {
        dudx(jnode,jlev) += add*avgS(iedge,jlev,LON);
        dvdy(jnode,jlev) += add*avgS(iedge,jlev,LAT);
      }
    }
    double y  = lonlat_deg(jnode,LAT) * deg2rad;
    double hy = radius;
    for(size_t jlev = 0; jlev < nlev; ++jlev)
    {
      dudx(jnode,jlev) *= hy/V(jnode);
      dvdy(jnode,jlev) *= hy/V(jnode);
      div(jnode,jlev) = dudx(jnode,jlev) + dvdy(jnode,jlev);
    }
  }

  //io::Gmsh().write(*mesh,grid->shortName()+".msh");
//  fvm_->haloExchange(*dudx_field);
//  fvm_->haloExchange(*dvdy_field);
//  io::Gmsh().write(*dudx_field,"gg.N64_fields.msh",std::ios::app);
//  io::Gmsh().write(*dvdy_field,"gg.N64_fields.msh",std::ios::app);
}




} // namespace nabla
} // namespace numerics
} // namespace atlas
