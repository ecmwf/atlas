#include <memory>

#include "atlas/grid/regular/RegularRegional.h"
#include "atlas/grid/domain/RectangularDomain.h"
#include "atlas/grid/spacing/LinearSpacing.h"
#include "atlas/internals/Debug.h"

namespace atlas {
namespace grid {
namespace regular {

//-----------------------------------------------------------------------------

register_BuilderT1(Grid,RegularRegional,RegularRegional::grid_type_str());

std::string RegularRegional::grid_type_str() {
    return "regular_regional";
}

std::string RegularRegional::className() {
    return "atlas.grid.regular.RegularRegional";
}

std::string RegularRegional::shortName() const {
    std::ostringstream s;
    s << "RR"<< nlonmin() << "x" << nlat();
    return s.str();
}

//-----------------------------------------------------------------------------

RegularRegional::ParseUniformCentred::ParseUniformCentred(const eckit::Parametrisation& config) {
  valid = config.get("nx",nx)
      &&  config.get("ny",ny)
      &&  config.get("centre",centre_lonlat)
      &&  config.get("dx",dx)
      &&  config.get("dy",dy);
  endpoint_x = true;
  endpoint_y = true;
  config.get("endpoint_x",endpoint_x);
  config.get("endpoint_y",endpoint_y);
}

void RegularRegional::ParseUniformCentred::apply(RegularRegional& g) const {
  double p[] = {centre_lonlat[0],centre_lonlat[1]};
  g.projection_->lonlat2coords(p);
  
  NOTIMP;
  // g.spacing_x_.reset( new spacing::LinearSpacing(p[0],dx,nx,endpoint_x) );
  // g.spacing_y_.reset( new spacing::LinearSpacing(p[1],dy,ny,endpoint_y) );

  double LX = 2.*(p[0] - g.spacing_x_->front());
  double LY = 2.*(p[1] - g.spacing_y_->front());
  g.domain_.reset( new domain::RectangularDomain(
      { p[0]-0.5*LX , p[0]+0.5*LX },
      { p[1]-0.5*LY , p[1]+0.5*LY }));
}

//-----------------------------------------------------------------------------

RegularRegional::ParseBounds::ParseBounds(const eckit::Parametrisation& config) {
  valid = config.get("nx",nx)
      &&  config.get("ny",ny)
      &&  config.get("xmin",xmin)
      &&  config.get("xmax",xmax)
      &&  config.get("ymin",ymin)
      &&  config.get("ymax",ymax);
  endpoint_x = true;
  endpoint_y = true;
  config.get("endpoint_x",endpoint_x);
  config.get("endpoint_y",endpoint_y);
}

void RegularRegional::ParseBounds::apply(RegularRegional& g) const {
  g.spacing_x_.reset( new spacing::LinearSpacing( xmin, xmax, nx, endpoint_x ) );
  g.spacing_y_.reset( new spacing::LinearSpacing( ymin, ymax, ny, endpoint_y ) );
  g.domain_.reset( new domain::RectangularDomain( {xmin, xmax}, {ymin, ymax} ) );
}

//-----------------------------------------------------------------------------

RegularRegional::ParseLonLatBounds::ParseLonLatBounds(const eckit::Parametrisation& config) {
  valid = config.get("nx",nx)
      &&  config.get("ny",ny)
      &&  config.get("north",north) // unrotated
      &&  config.get("south",south) // unrotated
      &&  config.get("east",east)   // unrotated
      &&  config.get("west",west);  // unrotated
  endpoint_x = true;
  endpoint_y = true;
  config.get("endpoint_x",endpoint_x);
  config.get("endpoint_y",endpoint_y);
}


//template<typename T>
//bool is_projection(const projection::Projection& p) {
//  return dynamic_cast<const T*>(&p);
//}

//void RegularRegional::ParseLonLatBounds::apply(RegularRegional& g) const {
//  if( not ( is_projection<projection::LonLatProjection>       (g.projection() )||
//            is_projection<projection::RotatedLonLatProjection>(g.projection()))) {
//    throw eckit::BadParameter("Projection is not compatible with lon/lat bounds",Here());
//  }
//  g.spacing_x_.reset( new spacing::LinearSpacing( {west,east}, nx, endpoint_x ) );
//  g.spacing_y_.reset( new spacing::LinearSpacing( {south,north}, ny, endpoint_y ) );
//  g.domain_.reset( new domain::RectangularDomain( {west,east}, {south,north} ) );
//}


void RegularRegional::ParseLonLatBounds::apply(RegularRegional& g) const {
  double sw[] = {west,south};
  double ne[] = {east,north};
  g.projection_->lonlat2coords(sw);
  g.projection_->lonlat2coords(ne);
  g.spacing_x_.reset( new spacing::LinearSpacing( sw[0],ne[0], nx, endpoint_x ) );
  g.spacing_y_.reset( new spacing::LinearSpacing( sw[1],ne[1], ny, endpoint_y ) );
  g.domain_.reset( new domain::RectangularDomain( {sw[0],ne[0]}, {sw[1],ne[1]} ) );
}

//-----------------------------------------------------------------------------

void RegularRegional::setup(const util::Config& config) {

    // read projection subconfiguration
    util::Config config_proj;
    if ( not config.get("projection",config_proj) ) {
      config_proj.set("type","lonlat");
    }
    projection_.reset( projection::Projection::create(config_proj) );

    using Parsers = std::vector< std::unique_ptr<Parse> >;
    Parsers parsers;
    parsers.push_back( std::unique_ptr<ParseUniformCentred>( new ParseUniformCentred(config) ) );
    parsers.push_back( std::unique_ptr<ParseLonLatBounds>( new ParseLonLatBounds(config) ) );
    parsers.push_back( std::unique_ptr<ParseBounds>( new ParseBounds(config) ) );

    bool configured(false);
    for( Parsers::const_iterator it=parsers.begin(); it!=parsers.end(); ++it ) {
      const Parse& p = **it;
      if( not configured && p.valid ) {
        p.apply(*this);
        configured = true;
      }
    }

    if( not configured ) {
      throw eckit::BadParameter("Could not configure RegularRegional grid", Here());
    }

    // setup regular grid
    Regular::setup();

}

RegularRegional::RegularRegional(const util::Config& config) :
    Regular() {
    setup(config);
}

RegularRegional::RegularRegional() :
    Regular() {
}

//-----------------------------------------------------------------------------

}  // namespace regular
}  // namespace grid
}  // namespace atlas

