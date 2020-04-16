/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "Regional.h"

#include "atlas/grid/StructuredGrid.h"
#include "atlas/grid/detail/grid/GridBuilder.h"
#include "atlas/runtime/Exception.h"
#include "atlas/runtime/Log.h"

using atlas::grid::LinearSpacing;
using XSpace = atlas::StructuredGrid::XSpace;
using YSpace = atlas::StructuredGrid::YSpace;

namespace atlas {
namespace grid {
namespace {  // anonymous

static Domain domain( const Grid::Config& grid ) {
    Grid::Config config;
    if ( grid.get( "domain", config ) ) {
        return Domain( config );
    }
    return Domain();
}

struct ConfigParser {
    struct Parsed {
        Parsed() = default;
        Parsed( std::initializer_list<double> interval ) : min( *interval.begin() ), max( *( interval.begin() + 1 ) ) {}
        double min;
        double max;
        long N;
        double step;
        bool endpoint = {true};
    };
    bool valid = {false};
    Parsed x;
    Parsed y;

    double step( double min, double max, long N, bool endpoint = true ) {
        double l = max - min;
        if ( endpoint && N > 1 ) {
            return l / double( N - 1 );
        }
        else {
            return l / double( N );
        }
    }

    static bool parse( const Projection&, const Grid::Config&, Parsed& x, Parsed& y );

    template <typename Parser>
    static bool parse( const Projection&, const Grid::Config&, Parsed& x, Parsed& y );
};

struct Parse_llc_step : ConfigParser {
    Parse_llc_step( const Projection& p, const Grid::Config& config ) {
        std::vector<double> centre_lonlat;

        valid = config.get( "nx", x.N ) && config.get( "ny", y.N ) && config.get( "dx", x.step ) &&
                config.get( "dy", y.step ) && config.get( "lonlat(centre)", centre_lonlat );

        if ( not valid ) {
            return;
        }

        double centre[] = {centre_lonlat[0], centre_lonlat[1]};
        p.lonlat2xy( centre );

        double lx = x.step * double( x.N - 1 );
        double ly = y.step * double( y.N - 1 );

        x.min = centre[0] - 0.5 * lx;
        x.max = centre[0] + 0.5 * lx;
        y.min = centre[1] - 0.5 * ly;
        y.max = centre[1] + 0.5 * ly;
    }
};

struct Parse_bounds_xy : ConfigParser {
    Parse_bounds_xy( const Projection& p, const Grid::Config& config ) {
        valid = config.get( "nx", x.N ) && config.get( "ny", y.N ) && config.get( "xmin", x.min ) &&
                config.get( "xmax", x.max ) && config.get( "ymin", y.min ) && config.get( "ymax", y.max );

        if ( not valid ) {
            return;
        }

        x.step = step( x.min, x.max, x.N );
        y.step = step( y.min, y.max, y.N );
    }
};

struct Parse_bounds_lonlat : ConfigParser {
    Parse_bounds_lonlat( const Projection& p, const Grid::Config& config ) {
        valid = config.get( "nx", x.N ) && config.get( "ny", y.N ) && config.get( "north", y.max )  // unrotated!
                && config.get( "south", y.min )                                                     // unrotated!
                && config.get( "east", x.max )                                                      // unrotated!
                && config.get( "west", x.min );                                                     // unrotated!

        // This version only works with a "lonlat" or "rotated_lonlat" projection!!!
        if ( valid ) {
            bool valid_projection = not p or p.type() == "rotated_lonlat";
            if ( not valid_projection ) {
                std::stringstream errmsg;
                errmsg << "This configuration requires that the projection is "
                          "\"lonlat\" or \"rotated_lonlat\". Received: "
                       << p.type();
                errmsg << "\n"
                          "p.bool() = "
                       << bool( p );
                throw_Exception( errmsg.str(), Here() );
            }
        }

        if ( not valid ) {
            return;
        }

        x.step = step( x.min, x.max, x.N );
        y.step = step( y.min, y.max, y.N );
    }
};

struct Parse_ll00_ll11 : ConfigParser {
    Parse_ll00_ll11( const Projection& p, const Grid::Config& config ) {
        std::vector<double> sw;
        std::vector<double> ne;
        valid = config.get( "nx", x.N ) && config.get( "ny", y.N ) &&
                config.get( "lonlat(xmin,ymin)", sw )      // includes rotation
                && config.get( "lonlat(xmax,ymax)", ne );  // includes rotation

        if ( not valid ) {
            return;
        }

        p.lonlat2xy( sw.data() );
        p.lonlat2xy( ne.data() );
        x.min = sw[0];
        x.max = ne[0];
        y.min = sw[1];
        y.max = ne[1];

        x.step = step( x.min, x.max, x.N );
        y.step = step( y.min, y.max, y.N );
    }
};

struct Parse_xy01_step : ConfigParser {
    Parse_xy01_step( const Projection&, const Grid::Config& config ) {
        valid = config.get( "nx", x.N ) && config.get( "ny", y.N ) && config.get( "dx", x.step ) &&
                config.get( "dy", y.step ) && config.get( "xmin", x.min ) &&
                config.get( "ymax", y.max );  // includes rotation

        if ( not valid ) {
            return;
        }

        x.max = x.min + x.step * ( x.N - 1 );
        y.min = y.max - y.step * ( y.N - 1 );
    }
};

struct Parse_ll00_step : ConfigParser {
    Parse_ll00_step( const Projection& p, const Grid::Config& config ) {
        std::vector<double> sw;
        valid = config.get( "nx", x.N ) && config.get( "ny", y.N ) && config.get( "dx", x.step ) &&
                config.get( "dy", y.step ) && config.get( "lonlat(xmin,ymin)", sw );  // includes rotation

        if ( not valid ) {
            return;
        }

        p.lonlat2xy( sw.data() );
        x.min = sw[0];
        y.min = sw[1];

        x.max = x.min + x.step * ( x.N - 1 );
        y.max = y.min + y.step * ( y.N - 1 );
    }
};

struct Parse_xy00_step : ConfigParser {
    Parse_xy00_step( const Projection&, const Grid::Config& config ) {
        valid = config.get( "nx", x.N ) && config.get( "ny", y.N ) && config.get( "dx", x.step ) &&
                config.get( "dy", y.step ) && config.get( "xmin", x.min ) && config.get( "ymin", y.min );

        if ( not valid ) {
            return;
        }

        x.max = x.min + x.step * ( x.N - 1 );
        y.max = y.min + y.step * ( y.N - 1 );
    }
};


struct Parse_ll01_step : ConfigParser {  // This resembles GRIB input for "lambert_azimutal_equal_area"
    Parse_ll01_step( const Projection& p, const Grid::Config& config ) {
        std::vector<double> nw;
        valid = config.get( "nx", x.N ) && config.get( "ny", y.N ) && config.get( "dx", x.step ) &&
                config.get( "dy", y.step ) && config.get( "lonlat(xmin,ymax)", nw );  // includes rotation

        if ( not valid ) {
            return;
        }

        p.lonlat2xy( nw.data() );
        x.min = nw[0];
        y.max = nw[1];

        x.max = x.min + x.step * ( x.N - 1 );
        y.min = y.max - y.step * ( y.N - 1 );
    }
};


template <typename Parser>
bool ConfigParser::parse( const Projection& projection, const Grid::Config& config, Parsed& x, Parsed& y ) {
    Parser p( projection, config );
    if ( p.valid ) {
        x = p.x;
        y = p.y;
        return true;  // success
    }
    return false;  // failure
}

bool ConfigParser::parse( const Projection& projection, const Grid::Config& config, Parsed& x, Parsed& y ) {
    // bounding box using 4 variables  (any projection allowed)
    if ( ConfigParser::parse<Parse_bounds_xy>( projection, config, x, y ) ) {
        return true;
    }

    // centre of domain and increments  (any projection allowed)
    if ( ConfigParser::parse<Parse_llc_step>( projection, config, x, y ) ) {
        return true;
    }

    // top-left of domain and increments (any projection allowed)
    if ( ConfigParser::parse<Parse_ll01_step>( projection, config, x, y ) ) {
        return true;
    }
    if ( ConfigParser::parse<Parse_xy01_step>( projection, config, x, y ) ) {
        return true;
    }

    // bottom-left of domain and increments (any projection allowed)
    if ( ConfigParser::parse<Parse_ll00_step>( projection, config, x, y ) ) {
        return true;
    }
    if ( ConfigParser::parse<Parse_xy00_step>( projection, config, x, y ) ) {
        return true;
    }

    // bounding box using two points defined in lonlat (any projection allowed)
    if ( ConfigParser::parse<Parse_ll00_ll11>( projection, config, x, y ) ) {
        return true;
    }

    //-------- From here on, projection must be (rotated) lonlat --------//

    // bounding box using 4 variables (south west north east)
    if ( ConfigParser::parse<Parse_bounds_lonlat>( projection, config, x, y ) ) {
        return true;
    }

    return false;
}

static class regional : public GridBuilder {
public:
    regional() : GridBuilder( "regional" ) {}

    void print( std::ostream& os ) const override {
        // os << std::left << std::setw(20) << "O<gauss>" << "Octahedral Gaussian
        // grid";
    }

    const Grid::Implementation* create( const std::string& name, const Grid::Config& config ) const override {
        throw_NotImplemented( "There are no named regional grids implemented.", Here() );
        return nullptr;
    }

    const Grid::Implementation* create( const Grid::Config& config ) const override {
        // read projection subconfiguration
        Projection projection;
        {
            util::Config config_proj;
            if ( config.get( "projection", config_proj ) ) {
                projection = Projection( config_proj );
            }
        }

        // Read grid configuration
        ConfigParser::Parsed x, y;
        if ( not ConfigParser::parse( projection, config, x, y ) ) {
            throw_Exception( "Could not parse configuration for RegularRegional grid", Here() );
        }

        YSpace yspace = config.getInt( "y_numbering", -1 ) < 0 ? LinearSpacing( y.max, y.min, y.N, y.endpoint )
                                                               : LinearSpacing( y.min, y.max, y.N, y.endpoint );

        bool with_endpoint = true;
        XSpace xspace( {x.min, x.max}, std::vector<long>( y.N, x.N ), with_endpoint );
        return new StructuredGrid::grid_t( xspace, yspace, projection, domain( config ) );
    }

    void force_link() {}

} regional_;

static class zonal_band : public GridBuilder {
public:
    zonal_band() : GridBuilder( "zonal_band" ) {}

    void print( std::ostream& os ) const override {
        // os << std::left << std::setw(20) << "O<gauss>" << "Octahedral Gaussian
        // grid";
    }

    const Grid::Implementation* create( const std::string& name, const Grid::Config& config ) const override {
        throw_NotImplemented( "There are no named zonal_band grids implemented.", Here() );
        return nullptr;
    }

    const Grid::Implementation* create( const Grid::Config& config ) const override {
        // read projection subconfiguration
        Projection projection;
        {
            util::Config config_proj;
            if ( config.get( "projection", config_proj ) ) {
                projection = Projection( config_proj );
            }
        }

        ATLAS_ASSERT( projection.units() == "degrees" );

        // Read grid configuration
        ConfigParser::Parsed y;
        long nx;
        if ( not config.get( "nx", nx ) ) {
            throw_Exception( "Parameter 'nx' missing in configuration", Here() );
        }
        if ( not config.get( "ny", y.N ) ) {
            throw_Exception( "Parameter 'ny' missing in configuration", Here() );
        }
        if ( not( config.get( "ymin", y.min ) or config.get( "south", y.min ) ) ) {
            y.min = -90.;
        }
        if ( not( config.get( "ymax", y.max ) or config.get( "north", y.max ) ) ) {
            y.max = 90.;
        }

        YSpace yspace = config.getInt( "y_numbering", -1 ) < 0 ? LinearSpacing( y.max, y.min, y.N, true )
                                                               : LinearSpacing( y.min, y.max, y.N, true );

        XSpace xspace( {0., 360.}, std::vector<long>( y.N, nx ), false );

        return new StructuredGrid::grid_t( xspace, yspace, projection, domain( config ) );
    }

    void force_link() {}

} zonal_band_;

}  // namespace

namespace detail {
namespace grid {

void force_link_Regional() {
    regional_.force_link();
    zonal_band_.force_link();
}


static 
StructuredGrid::grid_t* lambertregional (long Nx, long Ny, double XMinInMeters, double YMinInMetres, double DxInMetres, double DyInMetres, 
                                         double LoVInDegrees, double LaDInDegrees, double Latin1InDegrees, double Latin2InDegrees)
{
  using namespace atlas::util;

  std::vector<Spacing> spacings (Ny);

  for (int i = 0; i < Ny; i++)
    spacings[i] = Spacing (Config ("type", "linear") | Config ("N", Nx) 
                         | Config ("start", XMinInMeters) | Config ("end", XMinInMeters + (Nx - 1) * DxInMetres));

  StructuredGrid::XSpace xspace (spacings);
  StructuredGrid::YSpace yspace (Config ("type", "linear") | Config ("N", Ny) | Config ("start", YMinInMetres) | Config ("end", YMinInMetres + (Ny - 1 ) * DyInMetres));
  Projection proj (Config ("type", "lambert_conformal_conic") | Config ("longitude0", LoVInDegrees)
                 | Config ("latitude0", LaDInDegrees) | Config ("latitude1", Latin1InDegrees) 
                 | Config ("latitude2", Latin2InDegrees));

  return new StructuredGrid::grid_t (xspace, yspace, proj, Domain ());
}

static 
StructuredGrid::grid_t* latlonregional (long Nx, long Ny, double xmin, double xmax, double ymin, double ymax)
{
  using namespace atlas::util;

  std::vector<Spacing> spacings (Ny);

  for (int i = 0; i < Ny; i++)
    spacings[i] = Spacing (Config ("type", "linear") | Config ("N", Nx) 
                         | Config ("start", xmin) | Config ("end", xmax));

  StructuredGrid::XSpace xspace (spacings);
  StructuredGrid::YSpace yspace (Config ("type", "linear") | Config ("N", Ny) | Config ("start", ymin) | Config ("end", ymax));
  Projection proj (Config ("type", "lonlat"));

  return new StructuredGrid::grid_t (xspace, yspace, proj, Domain ());
}


extern "C" {


StructuredGrid::grid_t* atlas__grid__LambertRegional_int ( int nx, int ny, double xmin, double ymin, double dx, double dy, 
                                                           double longitude0, double latitude0, double latitude1, double latitude2 ) {
    return lambertregional( nx, ny, xmin, ymin, dx, dy, longitude0, latitude0, latitude1, latitude2);
}

StructuredGrid::grid_t* atlas__grid__LambertRegional_long ( long nx, long ny, double xmin, double ymin, double dx, double dy, 
                                                            double longitude0, double latitude0, double latitude1, double latitude2 ) {
    return lambertregional( nx, ny, xmin, ymin, dx, dy, longitude0, latitude0, latitude1, latitude2);
}

StructuredGrid::grid_t* atlas__grid__LatLonRegional_int ( int nx, int ny, double xmin, double xmax, double ymin, double ymax) {
    return latlonregional (nx, ny, xmin, xmax, ymin, ymax);
}

StructuredGrid::grid_t* atlas__grid__LatLonRegional_long ( long nx, long ny, double xmin, double xmax, double ymin, double ymax) {
    return latlonregional (nx, ny, xmin, xmax, ymin, ymax);
}


}



}  // namespace grid
}  // namespace detail

}  // namespace grid
}  // namespace atlas
