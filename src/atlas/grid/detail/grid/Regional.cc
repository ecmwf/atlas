#include "Regional.h"

#include "atlas/grid/detail/domain/RectangularDomain.h"
#include "atlas/grid/detail/domain/ZonalBandDomain.h"
#include "atlas/grid/detail/spacing/LinearSpacing.h"
#include "atlas/grid/detail/grid/GridBuilder.h"

using atlas::grid::domain::RectangularDomain;
using atlas::grid::domain::ZonalBandDomain;
using atlas::grid::spacing::LinearSpacing;
using XSpace = atlas::grid::StructuredGrid::grid_t::XSpace;
using YSpace = atlas::grid::StructuredGrid::YSpace;

namespace atlas {
namespace grid {
namespace { // anonymous

struct ConfigParser {

  struct Parsed {
    Parsed() {}
    Parsed( std::initializer_list<double> interval ) : min(*interval.begin()), max(*(interval.begin()+1)) {}
    double min;
    double max;
    long N;
    double step;
    bool endpoint = {true};
  };
  bool valid = {false};
  Parsed x;
  Parsed y;

  double step(double min,double max, long N, bool endpoint=true) {
    double l = max-min;
    if( endpoint && N>1 )
      return l/double(N-1);
    else
      return l/double(N);
  }

  static bool parse( const Projection&, const Grid::Config&, Parsed& x, Parsed& y );

  template <typename Parser>
  static bool parse( const Projection&, const Grid::Config&, Parsed& x, Parsed& y );

};

struct Parse_llc_step : ConfigParser {
  Parse_llc_step( const Projection& p, const Grid::Config& config ) {
      std::vector<double> centre_lonlat;

      valid = config.get("nx",x.N)
          &&  config.get("ny",y.N)
          &&  config.get("dx",x.step)
          &&  config.get("dy",y.step)
          &&  config.get("lonlat(centre)",centre_lonlat);

      if( not valid) return;

      double centre[] = {centre_lonlat[0],centre_lonlat[1]};
      p.lonlat2xy(centre);

      double lx = x.step * double(x.N-1);
      double ly = y.step * double(y.N-1);

      x.min = centre[0] - 0.5 * lx;
      x.max = centre[0] + 0.5 * lx;
      y.min = centre[1] - 0.5 * ly;
      y.max = centre[1] + 0.5 * ly;

  }
};

struct Parse_bounds_xy : ConfigParser {
  Parse_bounds_xy( const Projection& p, const Grid::Config& config ) {
    valid = config.get("nx",x.N)
        &&  config.get("ny",y.N)
        &&  config.get("xmin",x.min)
        &&  config.get("xmax",x.max)
        &&  config.get("ymin",y.min)
        &&  config.get("ymax",y.max);

    if( not valid ) return;

    x.step = step(x.min,x.max,x.N);
    y.step = step(y.min,y.max,y.N);
  }
};

struct Parse_bounds_lonlat : ConfigParser {
  Parse_bounds_lonlat( const Projection& p, const Grid::Config& config ) {
    valid = config.get("nx",x.N)
        &&  config.get("ny",y.N)
        &&  config.get("north",y.max)  // unrotated!
        &&  config.get("south",y.min)  // unrotated!
        &&  config.get("east", x.max)  // unrotated!
        &&  config.get("west", x.min); // unrotated!

    // This version only works with a "lonlat" or "rotated_lonlat" projection!!!
    if( valid ) {
      bool valid_projection = p || p.type() == "rotated_lonlat";
      if( not valid_projection ) {
        throw eckit::BadParameter("This configuration requires that the projection is \"lonlat\" or \"rotated_lonlat\"",Here());
      }
    }

    if( not valid ) return;

    x.step = step(x.min,x.max,x.N);
    y.step = step(y.min,y.max,y.N);
  }
};

struct Parse_ll00_ll11 : ConfigParser {
  Parse_ll00_ll11( const Projection& p, const Grid::Config& config ) {
    std::vector<double> sw;
    std::vector<double> ne;
    valid = config.get("nx",x.N)
        &&  config.get("ny",y.N)
        &&  config.get("lonlat(xmin,ymin)",sw)   // includes rotation
        &&  config.get("lonlat(xmax,ymax)",ne);  // includes rotation

    if( not valid ) return;

    p.lonlat2xy(sw.data());
    p.lonlat2xy(ne.data());
    x.min = sw[0];    x.max = ne[0];
    y.min = sw[1];    y.max = ne[1];

    x.step = step(x.min,x.max,x.N);
    y.step = step(y.min,y.max,y.N);
  }
};

struct Parse_ll00_step : ConfigParser {
  Parse_ll00_step( const Projection& p, const Grid::Config& config ) {
    std::vector<double> sw;
    valid = config.get("nx",x.N)
        &&  config.get("ny",y.N)
        &&  config.get("dx",x.step)
        &&  config.get("dy",y.step)
        &&  config.get("lonlat(xmin,ymin)",sw); // includes rotation

    if( not valid ) return;

    p.lonlat2xy(sw.data());
    x.min = sw[0];
    y.min = sw[1];

    x.max = x.min + x.step * (x.N-1);
    y.max = y.min + y.step * (y.N-1);
  }
};


template <typename Parser>
bool ConfigParser::parse( const Projection& projection, const Grid::Config& config, Parsed& x, Parsed& y ) {
  Parser p(projection,config);
  if( p.valid ) {
    x = p.x;
    y = p.y;
    return true; // success
  }
  return false; // failure
}

bool ConfigParser::parse( const Projection& projection, const Grid::Config& config, Parsed& x, Parsed& y ) {

  // bounding box using 4 variables  (any projection allowed)
  if( ConfigParser::parse< Parse_bounds_xy >( projection, config, x, y) ) return true;

  // centre of domain and increments  (any projection allowed)
  if( ConfigParser::parse< Parse_llc_step >( projection, config, x, y) ) return true;

  // bottom-left of domain and increments (any projection allowed)
  if( ConfigParser::parse< Parse_ll00_step >( projection, config, x, y) ) return true;

  // bounding box using two points defined in lonlat (any projection allowed)
  if( ConfigParser::parse< Parse_ll00_ll11 >( projection, config, x, y) ) return true;

// From here on, projection must be (rotated) lonlat

  // bounding box using 4 variables (south west north east)
  if( ConfigParser::parse< Parse_bounds_lonlat >( projection, config, x, y) ) return true;

  return false;
}


static class regional : public GridBuilder {

public:

  regional() : GridBuilder("regional") {}

  virtual void print(std::ostream& os) const {
    //os << std::left << std::setw(20) << "O<gauss>" << "Octahedral Gaussian grid";
  }

  virtual const Grid::grid_t* create( const std::string& name, const Grid::Config& config ) const {
    eckit::NotImplemented( "There are no named regional grids implemented.", Here() );
    return nullptr;
  }

  virtual const Grid::grid_t* create( const Grid::Config& config ) const {


    // read projection subconfiguration
    Projection projection;
    {
      util::Config config_proj;
      if ( config.get("projection",config_proj) ) {
        projection = Projection( config_proj );
      }
    }

    // Read grid configuration
    ConfigParser::Parsed x,y;
    if( not ConfigParser::parse(projection,config,x,y) ) {
        throw eckit::BadParameter("Could not parse configuration for RegularRegional grid", Here());
    }

    // Deduce the domain
    Domain domain ( new RectangularDomain( {x.min,x.max}, {y.min,y.max}, projection.units() ) );

    YSpace yspace( new LinearSpacing(y.min,y.max,y.N,y.endpoint) );

    bool with_endpoint = true;
    XSpace* xspace( new XSpace( {x.min,x.max}, std::vector<long>(y.N,x.N), with_endpoint ) );

    return new StructuredGrid::grid_t( projection, xspace, yspace, domain );

  }

} regional_;




static class zonal_band : public GridBuilder {

public:

  zonal_band() : GridBuilder("zonal_band") {}

  virtual void print(std::ostream& os) const {
    //os << std::left << std::setw(20) << "O<gauss>" << "Octahedral Gaussian grid";
  }

  virtual const Grid::grid_t* create( const std::string& name, const Grid::Config& config ) const {
    eckit::NotImplemented( "There are no named zonal_band grids implemented.", Here() );
    return nullptr;
  }

  virtual const Grid::grid_t* create( const Grid::Config& config ) const {


    // read projection subconfiguration
    Projection projection;
    {
      util::Config config_proj;
      if ( config.get("projection",config_proj) ) {
        projection = Projection( config_proj );
      }
    }

    ASSERT( projection.units() == "degrees" );

    // Read grid configuration
    ConfigParser::Parsed y;
    long nx;
    if( not config.get("nx",nx ) ) throw eckit::BadParameter("Parameter 'nx' missing in configuration", Here());
    if( not config.get("ny",y.N) ) throw eckit::BadParameter("Parameter 'ny' missing in configuration", Here());
    if( not (config.get("ymin",y.min) or config.get("south",y.min)) ) y.min = -90.;
    if( not (config.get("ymax",y.max) or config.get("north",y.max)) ) y.max =  90.;

    // Deduce the domain
    Domain domain ( new ZonalBandDomain( {y.min,y.max} ) );

    YSpace yspace( new LinearSpacing(y.min,y.max,y.N,true) );

    XSpace* xspace( new XSpace( {0.,360.}, std::vector<long>(y.N,nx), false ) );

    return new StructuredGrid::grid_t( projection, xspace, yspace, domain );

  }

} zonal_band_;




} // anonymous
} // namespace grid
} // namespace atlas
