/*
 * (C) Copyright 1996-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "eckit/exception/Exceptions.h"
#include "eckit/config/Parametrisation.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/mesh/Nodes.h"
#include "atlas/mesh/HybridElements.h"
#include "atlas/field/Field.h"
#include "atlas/numerics/fvm/Nabla.h"
#include "atlas/numerics/fvm/Method.h"
#include "atlas/internals/Parameters.h"
#include "atlas/internals/Bitflags.h"
#include "atlas/util/array/ArrayView.h"
#include "atlas/util/array/IndexView.h"
#include "atlas/util/parallel/atlas_omp.h"
#include "atlas/util/io/Gmsh.h"

// =======================================================

using atlas::internals::Topology;

namespace atlas {
namespace numerics {
namespace fvm {

namespace {
static NablaBuilder< Nabla > __fvm_nabla("fvm");
}

Nabla::Nabla(const numerics::Method &method, const eckit::Parametrisation &p) :
  atlas::numerics::Nabla(method,p)
{
  fvm_ = dynamic_cast<const fvm::Method *>(&method);
  if( ! fvm_ )
    throw eckit::BadCast("atlas::numerics::fvm::Nabla needs a atlas::numerics::fvm::Method",Here());
  Log::info() << "Nabla constructed for method " << fvm_->name()
                     << " with " << fvm_->functionspace_nodes().nb_nodes_global() << " nodes total" << std::endl;

  setup();

}

Nabla::~Nabla()
{
}

void Nabla::setup()
{
  const mesh::Edges &edges = fvm_->mesh().edges();

  const size_t nedges = edges.size();

  const util::array::ArrayView<int,1> edge_is_pole ( edges.field("is_pole_edge") );

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


void Nabla::gradient(const field::Field& scalar_field, field::Field& grad_field) const
{
  const double radius = fvm_->radius();
  const double deg2rad = M_PI/180.;

  const mesh::Edges &edges = fvm_->mesh().edges();
  const mesh::Nodes &nodes = fvm_->mesh().nodes();

  const size_t nnodes = nodes.size();
  const size_t nedges = edges.size();
  const size_t nlev = scalar_field.levels();
  if( grad_field.levels() != nlev )
    throw eckit::AssertionFailed("gradient field should have same number of levels",Here());


  const util::array::ArrayView<double,2> scalar ( scalar_field.data<double>(), util::array::make_shape(nnodes,nlev)   );
        util::array::ArrayView<double,3> grad   ( grad_field.  data<double>(), util::array::make_shape(nnodes,nlev,2) );

  const util::array::ArrayView<double,2> lonlat_deg     ( nodes.lonlat() );
  const util::array::ArrayView<double,1> V              ( nodes.field("dual_volumes") );
  const util::array::ArrayView<double,2> S              ( edges.field("dual_normals") );
  const util::array::ArrayView<int,   1> edge_is_pole   ( edges.field("is_pole_edge") );
  const util::array::ArrayView<double,2> node2edge_sign ( nodes.field("node2edge_sign") );

  const mesh::Connectivity& node2edge = nodes.edge_connectivity();
  const mesh::Connectivity& edge2node = edges.node_connectivity();

  util::array::ArrayT<double> avgS_arr( nedges,nlev,2 );
  util::array::ArrayView<double,3> avgS(avgS_arr);

  atlas_omp_parallel
  {
    atlas_omp_for( size_t jedge=0; jedge<nedges; ++jedge )
    {
      int ip1 = edge2node(jedge,0);
      int ip2 = edge2node(jedge,1);

      #pragma ivdep
      for(size_t jlev = 0; jlev < nlev; ++jlev)
      {
        double avg = ( scalar(ip1,jlev) + scalar(ip2,jlev) ) * 0.5;
        avgS(jedge,jlev,internals::LON) = S(jedge,internals::LON)*avg;
        avgS(jedge,jlev,internals::LAT) = S(jedge,internals::LAT)*avg;
      }
    }

    atlas_omp_for( size_t jnode=0; jnode<nnodes; ++jnode )
    {
      #pragma ivdep
      for(size_t jlev = 0; jlev < nlev; ++jlev )
      {
        grad(jnode,jlev,internals::LON) = 0.;
        grad(jnode,jlev,internals::LAT) = 0.;
      }
      for( size_t jedge=0; jedge<node2edge.cols(jnode); ++jedge )
      {
        const int iedge = node2edge(jnode,jedge);
        const double add = node2edge_sign(jnode,jedge);
        #pragma ivdep
        for(size_t jlev = 0; jlev < nlev; ++jlev)
        {
          grad(jnode,jlev,internals::LON) += add*avgS(iedge,jlev,internals::LON);
          grad(jnode,jlev,internals::LAT) += add*avgS(iedge,jlev,internals::LAT);
        }
      }
      const double y  = lonlat_deg(jnode,internals::LAT) * deg2rad;
      const double metric_x = radius/V(jnode);
      const double metric_y = metric_x*std::cos(y);
      #pragma ivdep
      for(size_t jlev = 0; jlev < nlev; ++jlev)
      {
        grad(jnode,jlev,internals::LON) *= metric_x;
        grad(jnode,jlev,internals::LAT) *= metric_y;
      }
    }
  }
}

// ================================================================================

void Nabla::divergence(const field::Field& vector_field, field::Field& div_field) const
{
  const double radius = fvm_->radius();
  const double deg2rad = M_PI/180.;

  const mesh::Edges &edges = fvm_->mesh().edges();
  const mesh::Nodes &nodes = fvm_->mesh().nodes();

  const size_t nnodes = nodes.size();
  const size_t nedges = edges.size();
  const size_t nlev = vector_field.levels();
  if( div_field.levels() != nlev )
    throw eckit::AssertionFailed("divergence field should have same number of levels",Here());

  const util::array::ArrayView<double,3> vector ( vector_field.data<double>(), util::array::make_shape(nnodes,nlev,2));
        util::array::ArrayView<double,2> div    ( div_field   .data<double>(), util::array::make_shape(nnodes,nlev)  );

  const util::array::ArrayView<double,2> lonlat_deg     ( nodes.lonlat() );
  const util::array::ArrayView<double,1> V              ( nodes.field("dual_volumes") );
  const util::array::ArrayView<double,2> S              ( edges.field("dual_normals") );
  const util::array::ArrayView<int,   1> edge_is_pole   ( edges.field("is_pole_edge") );
  const util::array::ArrayView<double,2> node2edge_sign ( nodes.field("node2edge_sign") );
  const mesh::Connectivity& node2edge = nodes.edge_connectivity();
  const mesh::Connectivity& edge2node = edges.node_connectivity();

  util::array::ArrayT<double> avgS_arr( nedges,nlev,2 );
  util::array::ArrayView<double,3> avgS(avgS_arr);

  atlas_omp_parallel
  {
    atlas_omp_for( size_t jedge=0; jedge<nedges; ++jedge )
    {
      int ip1 = edge2node(jedge,0);
      int ip2 = edge2node(jedge,1);
      double y1  = lonlat_deg(ip1,internals::LAT) * deg2rad;
      double y2  = lonlat_deg(ip2,internals::LAT) * deg2rad;
      double cosy1 = std::cos(y1);
      double cosy2 = std::cos(y2);

      double pbc = 1.-edge_is_pole(jedge);

      for(size_t jlev = 0; jlev < nlev; ++jlev)
      {
        double avg[2] = {
          (      vector(ip1,jlev,internals::LON) +       vector(ip2,jlev,internals::LON) ) * 0.5,
          (cosy1*vector(ip1,jlev,internals::LAT) + cosy2*vector(ip2,jlev,internals::LAT) ) * 0.5 * pbc // (force cos(y)=0 at pole)
        };
        avgS(jedge,jlev,internals::LON) = S(jedge,internals::LON)*avg[internals::LON];  // 0 at pole by construction of S
        avgS(jedge,jlev,internals::LAT) = S(jedge,internals::LAT)*avg[internals::LAT];  // 0 at pole by construction of pbc
        // We don't need the cross terms for divergence, i.e.  S(jedge,internals::LON)*avg[internals::LAT] / S(jedge,internals::LAT)*avg[internals::LON]
      }
    }

    atlas_omp_for( size_t jnode=0; jnode<nnodes; ++jnode )
    {
      for(size_t jlev = 0; jlev < nlev; ++jlev )
      {
        div(jnode,jlev) = 0.;
      }
      for( size_t jedge=0; jedge<node2edge.cols(jnode); ++jedge )
      {
        int iedge = node2edge(jnode,jedge);
        double add = node2edge_sign(jnode,jedge);
        for(size_t jlev = 0; jlev < nlev; ++jlev)
        {
          div(jnode,jlev) += add*(avgS(iedge,jlev,internals::LON)+avgS(iedge,jlev,internals::LAT));
        }
      }
      double metric = radius/V(jnode);
      for(size_t jlev = 0; jlev < nlev; ++jlev)
      {
        div(jnode,jlev) *= metric;
      }
    }
  }
}


void Nabla::curl(const field::Field& vector_field, field::Field& curl_field) const
{
  const double radius = fvm_->radius();
  const double deg2rad = M_PI/180.;

  const mesh::Edges &edges = fvm_->mesh().edges();
  const mesh::Nodes &nodes = fvm_->mesh().nodes();

  const size_t nnodes = nodes.size();
  const size_t nedges = edges.size();
  const size_t nlev = vector_field.levels();
  if( curl_field.levels() != nlev )
    throw eckit::AssertionFailed("curl field should have same number of levels",Here());

  const util::array::ArrayView<double,3> vector ( vector_field.data<double>(), util::array::make_shape(nnodes,nlev,2));
        util::array::ArrayView<double,2> curl   ( curl_field  .data<double>(), util::array::make_shape(nnodes,nlev)  );

  const util::array::ArrayView<double,2> lonlat_deg     ( nodes.lonlat() );
  const util::array::ArrayView<double,1> V              ( nodes.field("dual_volumes") );
  const util::array::ArrayView<double,2> S              ( edges.field("dual_normals") );
  const util::array::ArrayView<int,   1> edge_is_pole   ( edges.field("is_pole_edge") );
  const util::array::ArrayView<double,2> node2edge_sign ( nodes.field("node2edge_sign") );

  const mesh::Connectivity& node2edge = nodes.edge_connectivity();
  const mesh::Connectivity& edge2node = edges.node_connectivity();

  util::array::ArrayT<double> avgS_arr( nedges,nlev,2 );
  util::array::ArrayView<double,3> avgS(avgS_arr);

  atlas_omp_parallel
  {
    atlas_omp_for( size_t jedge=0; jedge<nedges; ++jedge )
    {
      int ip1 = edge2node(jedge,0);
      int ip2 = edge2node(jedge,1);
      double y1  = lonlat_deg(ip1,internals::LAT) * deg2rad;
      double y2  = lonlat_deg(ip2,internals::LAT) * deg2rad;
      double rcosy1 = radius*std::cos(y1);
      double rcosy2 = radius*std::cos(y2);

      double pbc = 1-edge_is_pole(jedge);

      for(size_t jlev = 0; jlev < nlev; ++jlev)
      {
        double avg[2] = {
          (rcosy1*vector(ip1,jlev,internals::LON) + rcosy2*vector(ip2,jlev,internals::LON) ) * 0.5 * pbc, // (force R*cos(y)=0 at pole)
          (radius*vector(ip1,jlev,internals::LAT) + radius*vector(ip2,jlev,internals::LAT) ) * 0.5
        };
        avgS(jedge,jlev,internals::LON) = S(jedge,internals::LAT)*avg[internals::LON]; // 0 at pole by construction of pbc
        avgS(jedge,jlev,internals::LAT) = S(jedge,internals::LON)*avg[internals::LAT]; // 0 at pole by construction of S
        // We don't need the non-cross terms for curl, i.e.  S(jedge,internals::LON)*avg[internals::LON] / S(jedge,internals::LAT)*avg[internals::LAT]
      }
    }

    atlas_omp_for( size_t jnode=0; jnode<nnodes; ++jnode )
    {
      for(size_t jlev = 0; jlev < nlev; ++jlev )
      {
        curl(jnode,jlev) = 0.;
      }
      for( size_t jedge=0; jedge<node2edge.cols(jnode); ++jedge )
      {
        int iedge = node2edge(jnode,jedge);
        double add = node2edge_sign(jnode,jedge);
        for(size_t jlev = 0; jlev < nlev; ++jlev)
        {
          curl(jnode,jlev) += add*(avgS(iedge,jlev,internals::LAT)-avgS(iedge,jlev,internals::LON));
        }
      }
      double metric = 1./V(jnode);
      for(size_t jlev = 0; jlev < nlev; ++jlev)
      {
        curl(jnode,jlev) *= metric;
      }
    }
  }
}


void Nabla::laplacian(const field::Field& scalar, field::Field& lapl) const
{
  eckit::SharedPtr<field::Field> grad (
       fvm_->functionspace_nodes().createField<double>("grad",
       scalar.levels(),util::array::make_shape(2)) );
  gradient(scalar,*grad);
  if( fvm_->functionspace_nodes().halo().size() < 2 )
    fvm_->functionspace_nodes().haloExchange(*grad);
  divergence(*grad,lapl);
}



} // namespace fvm
} // namespace numerics
} // namespace atlas
