/*
 * (C) Copyright 1996-2017 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "atlas/grid/Grid.h"

#include <limits>
#include <vector>
#include "eckit/config/Parametrisation.h"
#include "eckit/exception/Exceptions.h"

#include "atlas/grid/Grid.h"
#include "atlas/grid/Spacing.h"
#include "atlas/grid/Domain.h"
#include "atlas/grid/Projection.h"
#include "atlas/grid/detail/grid/Structured.h"
#include "atlas/grid/detail/spacing/LinearSpacing.h"
#include "atlas/grid/detail/spacing/CustomSpacing.h"
#include "atlas/util/Config.h"
#include "atlas/grid/detail/grid/types/Gaussian.h"

namespace atlas {
namespace grid {

Grid::Grid():
    grid_( nullptr ) {
}

Grid::Grid(const Grid& grid):
    grid_( grid.grid_ ) {
}

Grid::Grid( const Grid::grid_t *grid ):
    grid_( grid ) {
}

Grid::Grid( const std::string& shortname ) {
    grid_ = Grid::grid_t::create( shortname );
}

Grid::Grid( const Config& p ) {
    grid_ = Grid::grid_t::create(p);
}


StructuredGrid::grid_t* create_structured( const Grid::Config& config ) {

  Projection projection;
  Spacing    yspace;
  Domain     domain;

  Grid::Config config_proj;
  if( config.get("projection",config_proj) )
    projection = Projection(config_proj);

  Grid::Config config_yspace;
  if( not config.get("yspace",config_yspace) )
    throw eckit::BadParameter("yspace missing in configuration");
  yspace = Spacing(config_yspace);

  const size_t ny = yspace.size();

  StructuredGrid::grid_t::XSpace *X = new StructuredGrid::grid_t::XSpace(ny);

  double dom_xmin =  std::numeric_limits<double>::max();
  double dom_xmax = -std::numeric_limits<double>::max();

  std::vector<Grid::Config> config_xspace_list;
  if( config.get("xspace[]",config_xspace_list) ) {

    ASSERT( config_xspace_list.size() == ny );
    std::string xspace_type;

    for( size_t j=0; j<ny; ++j ) {
      config_xspace_list[j].get("type",xspace_type);
      ASSERT( xspace_type == "linear" );
      spacing::LinearSpacing::Params xspace( config_xspace_list[j] );
      X->xmin.push_back(xspace.start);
      X->xmax.push_back(xspace.end);
      X->nx.push_back(xspace.N);
      X->dx.push_back(xspace.step);
      dom_xmin = std::min(dom_xmin,xspace.start);
      dom_xmax = std::max(dom_xmax,xspace.end);
    }

  } else {

    Grid::Config config_xspace;
    if( not config.get("xspace",config_xspace) )
      throw eckit::BadParameter("xspace missing in configuration");

    std::string xspace_type;
    config_xspace.get("type",xspace_type);
    ASSERT( xspace_type == "linear" );

    std::vector<long>   v_N;
    std::vector<double> v_start;
    std::vector<double> v_end;
    std::vector<double> v_length;
    config_xspace.get("N[]",      v_N     );
    config_xspace.get("start[]",  v_start );
    config_xspace.get("end[]",    v_end   );
    config_xspace.get("length[]", v_length);

    if( not v_N.     empty() ) ASSERT(v_N.     size() == ny);
    if( not v_start. empty() ) ASSERT(v_start. size() == ny);
    if( not v_end.   empty() ) ASSERT(v_end.   size() == ny);
    if( not v_length.empty() ) ASSERT(v_length.size() == ny);

    for( size_t j=0; j<ny; ++j ) {
      if( not v_N.     empty() ) config_xspace.set("N",     v_N[j]);
      if( not v_start. empty() ) config_xspace.set("start", v_start[j]);
      if( not v_end.   empty() ) config_xspace.set("end",   v_end[j]);
      if( not v_length.empty() ) config_xspace.set("length",v_length[j]);
      spacing::LinearSpacing::Params xspace( config_xspace );
      X->xmin.push_back(xspace.start);
      X->xmax.push_back(xspace.end);
      X->nx.push_back(xspace.N);
      X->dx.push_back(xspace.step);
      dom_xmin = std::min(dom_xmin,xspace.start);
      dom_xmax = std::max(dom_xmax,xspace.end);
    }
  }

  Grid::Config config_domain;
  if( config.get("domain",config_domain) )
    domain = Domain(config_domain);
  else {
    config_domain.set("type","rectangular");
    config_domain.set("ymin",yspace.min());
    config_domain.set("ymax",yspace.max());
    config_domain.set("xmin",dom_xmin);
    config_domain.set("xmax",dom_xmax);
    config_domain.set("units",projection.units());
    domain = Domain(config_domain);
  }

  return new StructuredGrid::grid_t(projection, X, yspace, domain );
}

const StructuredGrid::grid_t* structured_grid( const Grid::grid_t *grid ) {
    const StructuredGrid::grid_t* g( dynamic_cast<const StructuredGrid::grid_t*>(grid) );
    return g;
}

StructuredGrid::StructuredGrid():
    Grid(),
    grid_( nullptr ) {
}

StructuredGrid::StructuredGrid( const Grid& grid ):
    Grid( grid ),
    grid_( structured_grid(get()) ) {
}

StructuredGrid::StructuredGrid( const Grid::grid_t* grid ):
    Grid( grid ),
    grid_( structured_grid(get()) ) {
}

StructuredGrid::StructuredGrid( const std::string& grid ):
    Grid( grid ),
    grid_( structured_grid(get()) ) {
}

StructuredGrid::StructuredGrid( const Config& p ):
    Grid( create_structured(p) ),
    grid_( structured_grid(get()) ) {
}

const RegularGrid::grid_t* regular_grid( const Grid::grid_t* grid ) {
    const RegularGrid::grid_t* g( dynamic_cast<const RegularGrid::grid_t*>(grid) );
    if( g && g->reduced() ) {
        return nullptr;
    }
    return g;
}

RegularGrid::grid_t* RegularGrid::create( const Config& ) {
  NOTIMP;
  // return nullptr;
}

RegularGrid::RegularGrid():
    StructuredGrid(),
    grid_( nullptr ),
    nx_(0) {
}

RegularGrid::RegularGrid( const Grid& grid ):
    StructuredGrid( grid ),
    grid_( regular_grid(get()) ) {
    if( grid_ ) nx_ = StructuredGrid::nx().front();
}

RegularGrid::RegularGrid( const detail::grid::Grid *grid ):
    StructuredGrid(grid),
    grid_( regular_grid(get()) ) {
    if( grid_ ) nx_ = StructuredGrid::nx().front();
}

RegularGrid::RegularGrid( const std::string& grid ):
    StructuredGrid(grid),
    grid_( regular_grid(get()) ) {
    if( grid_ ) nx_ = StructuredGrid::nx().front();
}


RegularGrid::RegularGrid( const Config& p ):
    StructuredGrid(create(p)),
    grid_( regular_grid(get()) ) {
    if( grid_ ) nx_ = StructuredGrid::nx().front();
}

ReducedGaussianGrid::ReducedGaussianGrid( const std::vector<long>& nx ):
    ReducedGaussianGrid::Grid( reduced_gaussian(nx) ) {
}

ReducedGaussianGrid::ReducedGaussianGrid( const std::initializer_list<long>& nx ):
    ReducedGaussianGrid( std::vector<long>(nx) ) {
}


} // namespace Grid
} // namespace atlas
