#include "atlas/grid/ptr/Grid.h"

#include <limits>
#include <vector>
#include "eckit/config/Parametrisation.h"
#include "eckit/exception/Exceptions.h"

#include "atlas/grid/Grid.h"
#include "atlas/grid/Structured.h"
#include "atlas/grid/Regular.h"
#include "atlas/grid/ptr/Spacing.h"
#include "atlas/grid/ptr/Domain.h"
#include "atlas/grid/ptr/Projection.h"
#include "atlas/grid/spacing/LinearSpacing.h"
#include "atlas/grid/spacing/CustomSpacing.h"
#include "atlas/util/Config.h"

namespace atlas {
namespace grid {
namespace ptr {


Grid::Grid():
    grid_( nullptr ),
    projection_( nullptr ),
    domain_( nullptr ) {
}

Grid::Grid(const Grid& grid):
    grid_( grid.grid_ ),
    projection_( grid.projection_ ),
    domain_(grid.domain_) {
}

Grid::Grid( const atlas::grid::Grid *grid ):
    grid_( grid ),
    projection_( nullptr ),
    domain_( nullptr ) {
    if( grid_ ) {
      projection_ = Projection( &grid_->projection() );
      domain_ = Domain( &grid_->domain() );
    }
}

Grid::Grid( const std::string& shortname ) {
    grid_ = atlas::grid::Grid::create( shortname );
    projection_ = Projection( &grid_->projection() );
    domain_ = Domain( &grid_->domain() );
}



atlas::grid::Structured* Structured::create( const eckit::Parametrisation &p ) {
  auto config = dynamic_cast<const util::Config&>(p);

  Projection projection;
  Spacing    yspace;
  Domain     domain;

  util::Config config_proj;
  if( config.get("projection",config_proj) )
    projection = Projection(config_proj);

  util::Config config_yspace;
  if( not config.get("yspace",config_yspace) )
    throw eckit::BadParameter("yspace missing in configuration");
  yspace = Spacing(config_yspace);

  const size_t ny = yspace.size();

  atlas::grid::Structured::XSpace *X = new atlas::grid::Structured::XSpace(ny);

  double dom_xmin =  std::numeric_limits<double>::max();
  double dom_xmax = -std::numeric_limits<double>::max();

  std::vector<util::Config> config_xspace_list;
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

    util::Config config_xspace;
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

  util::Config config_domain;
  if( dynamic_cast<const util::Config&>(config).get("domain",config_domain) )
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

  return new atlas::grid::Structured(projection, X, yspace, domain );
}

Structured::Structured():
    Grid(),
    grid_( nullptr ) {
}

Structured::Structured( const Grid& grid ):
    Grid( grid ),
    grid_(dynamic_cast<const atlas::grid::Structured*>(raw())) {
}

Structured::Structured( const atlas::grid::Grid *grid ):
    Grid( grid ),
    grid_(dynamic_cast<const atlas::grid::Structured*>(raw())) {
}

Structured::Structured( const eckit::Parametrisation& p ):
    Grid(create(p)),
    grid_(dynamic_cast<const atlas::grid::Structured*>(raw())) {
}


const atlas::grid::Structured* regular_grid( const atlas::grid::Grid *grid ) {
    const atlas::grid::Structured* g( dynamic_cast<const atlas::grid::Structured*>(grid) );
    if( g && g->reduced() ) {
        return nullptr;
    }
    return g;
}

atlas::grid::Structured* Regular::create( const eckit::Parametrisation& ) {
  NOTIMP;
}

Regular::Regular( const Grid& grid ):
    Structured( grid ),
    grid_(regular_grid(raw())) {
    if( grid_ ) nx_ = Structured::nx().front();
}

Regular::Regular( const atlas::grid::Grid *grid ):
    Structured(grid),
    grid_( regular_grid(raw()) ) {
    if( grid_ ) nx_ = Structured::nx().front();
}

Regular::Regular( const eckit::Parametrisation& p ):
    Structured(create(p)),
    grid_( regular_grid(raw()) ) {
    if( grid_ ) nx_ = Structured::nx().front();
}


} // namespace ptr
} // namespace Grid
} // namespace atlas
