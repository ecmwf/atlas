#include <memory>

#include "atlas/grid/regular/RegularRegional.h"
#include "atlas/grid/domain/RectangularDomain.h"
#include "atlas/grid/spacing/LinearSpacing.h"
#include "atlas/grid/projection/LonLatProjection.h"
#include "atlas/internals/Debug.h"

using atlas::grid::domain::RectangularDomain;
using atlas::grid::projection::Projection;
using atlas::grid::projection::RotatedLonLatProjection;
using atlas::grid::spacing::LinearSpacing;
using atlas::grid::spacing::Spacing;
using eckit::Parametrisation;

namespace atlas {
namespace grid {
namespace regular {

namespace {

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

    static bool parse( const Grid& grid, const Parametrisation& config, Parsed& x, Parsed& y );
    
    template <typename Parser>
    static bool parse( const Grid& grid, const Parametrisation& config, Parsed& x, Parsed& y );

  };

  struct Parse_llc_step : ConfigParser {
    Parse_llc_step( const Projection& p, const Parametrisation& config ) {
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
    Parse_bounds_xy( const Projection& p, const Parametrisation& config ) {
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
    Parse_bounds_lonlat( const Projection& p, const Parametrisation& config ) {
      valid = config.get("nx",x.N)
          &&  config.get("ny",y.N)
          &&  config.get("north",y.max)  // unrotated!
          &&  config.get("south",y.min)  // unrotated!
          &&  config.get("east", x.max)  // unrotated!
          &&  config.get("west", x.min); // unrotated!

      // This version only works with a "lonlat" or "rotated_lonlat" projection!!!
      if( valid ) {
        bool valid_projection = p || dynamic_cast<const RotatedLonLatProjection*>( &p );
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
    Parse_ll00_ll11( const Projection& p, const Parametrisation& config ) {
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
    Parse_ll00_step( const Projection& p, const Parametrisation& config ) {
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
  bool ConfigParser::parse( const Grid& g, const Parametrisation& config, Parsed& x, Parsed& y ) {
    Parser p(g.projection(),config);
    if( p.valid ) {
      x = p.x;
      y = p.y;
      return true; // success
    }
    return false; // failure
  }
  
  bool ConfigParser::parse( const Grid& g, const Parametrisation& config, Parsed& x, Parsed& y ) {

    // bounding box using 4 variables  (any projection allowed)
    if( ConfigParser::parse< Parse_bounds_xy >( g, config, x, y) ) return true;

    // centre of domain and increments  (any projection allowed)
    if( ConfigParser::parse< Parse_llc_step >( g, config, x, y) ) return true;

    // bottom-left of domain and increments (any projection allowed)
    if( ConfigParser::parse< Parse_ll00_step >( g, config, x, y) ) return true;

    // bounding box using two points defined in lonlat (any projection allowed)
    if( ConfigParser::parse< Parse_ll00_ll11 >( g, config, x, y) ) return true;

// From here on, projection must be (rotated) lonlat

    // bounding box using 4 variables (south west north east)
    if( ConfigParser::parse< Parse_bounds_lonlat >( g, config, x, y) ) return true;

    return false;
  }

}




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

void RegularRegional::setup(const util::Config& config) {

    // read projection subconfiguration
    {
      util::Config config_proj;
      if ( not config.get("projection",config_proj) ) {
        config_proj.set("type","lonlat");
      }
      projection_.reset( projection::Projection::create(config_proj) );
    } 

    // Read grid configuration
    ConfigParser::Parsed x,y;
    if( not ConfigParser::parse(*this,config,x,y) ) {
        throw eckit::BadParameter("Could not parse configuration for RegularRegional grid", Here());
    }
    
    // Deduce the domain
    domain_.reset( new RectangularDomain( {x.min,x.max}, {y.min,y.max}, projection_->units() ) );

    // Delegate further setup
    spacing::Spacing* yspace = new LinearSpacing(y.min,y.max,y.N,y.endpoint);
    Structured::setup(yspace,x.N,x.min,x.max,x.step);
}

RegularRegional::RegularRegional(const util::Config& config) :
    Regular() {
    setup(config);
}

//-----------------------------------------------------------------------------

}  // namespace regular
}  // namespace grid
}  // namespace atlas

