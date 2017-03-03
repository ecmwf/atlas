#include "LonLat.h"

#include "eckit/utils/Translator.h"

namespace atlas {
namespace grid {
namespace { // anonymous

static eckit::Translator<std::string,int> to_int;


struct Shift {

    enum Bits {
        NONE = 0,
        LAT  = (1<<1),
        LON  = (1<<2)
    };

    Shift(int bits=NONE) : bits_(bits) {
    }

    Shift(bool shift_lon, bool shift_lat) : bits_((shift_lon? LON:NONE) | (shift_lat? LAT:NONE)) {
    }

    bool operator()(int bits) const {
        return (bits_ & bits) == bits;
    }

    const int bits_;

};

using XSpace = detail::grid::Structured::XSpace;

StructuredGrid::grid_t* create_lonlat(long nlon, long nlat, Shift shift) {

    bool shifted_x = shift(Shift::LON);
    bool shifted_y = shift(Shift::LAT);
    
    double start_x = (shifted_x ? 0.5 : 0.0)*360.0/double(nlon);
    std::array<double,2> interval_x = { start_x, start_x+360. };
    bool no_endpoint = false;
    XSpace* xspace = new XSpace( interval_x, std::vector<long>(nlat,nlon), no_endpoint );

    // spacing is uniform in y
    Grid::Config config_spacing;
    config_spacing.set("type","linear");
    config_spacing.set("start", 90.0-(shifted_y ? 90.0/double(nlat) : 0.0) );
    config_spacing.set("end",-90.0+(shifted_y ? 90.0/double(nlat) : 0.0) );
    config_spacing.set("N",nlat);
    Spacing yspace(config_spacing);

    // domain is global
    Grid::Config config_domain;
    config_domain.set("type","global");
    Domain domain(config_domain);
    
    return new StructuredGrid::grid_t( Projection(), xspace, yspace, domain );
}

StructuredGrid::grid_t* create_lonlat( const Grid::Config& config, Shift shift ) {
  
  bool shifted_y = shift(Shift::LAT);
  
  long N, nx, ny;
  // dimensions
  if ( config.get("N",N) ) {
    nx = 4*N;
    ny = shifted_y ? 2*N : 2*N+1;
  } else if( config.get("nx",nx)
          && config.get("ny",ny) ) {
  } else {
      throw eckit::BadParameter("Configuration requires either N, or (nx,ny)",Here());
  }
  
  return create_lonlat(nx,ny,shift);
  
}

//---------------------------------------------------------------------------------------------------------------------

static class regular_lonlat : public GridCreator {

public:

  regular_lonlat(): GridCreator( "regular_lonlat", {
    "^[Ll]([0-9]+)x([0-9]+)$",
    "^[Ll]([0-9]+)$"           } ){}

  virtual void print(std::ostream& os) const {
    os << std::left << std::setw(20) << "L<nx>x<ny> / L<gauss>" << "Regular longitude-latitude grid";
  }

  virtual const Grid::grid_t* create( const std::string& name ) const {
    int id;
    std::vector<std::string> matches;
    if( match( name, matches, id ) ) {
      
      util::Config grid;
      grid.set("type", type());

      if( id == 0 ) {
        grid.set( "nx", to_int(matches[0]) );
        grid.set( "ny", to_int(matches[1]) );
        return create( grid );
      }

      if( id == 1 ) {
        grid.set( "N", to_int(matches[0]) );
        return create( grid );
      }

    }
    return nullptr;
  }
  
  virtual const Grid::grid_t* create( const Grid::Config& config ) const {
    return create_lonlat( config, Shift(false,false) );
  }

} regular_lonlat_;

//---------------------------------------------------------------------------------------------------------------------

static class shifted_lonlat : public GridCreator {

public:

  shifted_lonlat(): GridCreator( "shifted_lonlat", {
    "^[Ss]([0-9]+)x([0-9]+)$",
    "^[Ss]([0-9]+)$"           } ){}

  virtual void print(std::ostream& os) const {
    os << std::left << std::setw(20) << "S<nx>x<ny> / S<gauss>" << "Shifted longitude-latitude grid";
  }

  virtual const Grid::grid_t* create( const std::string& name ) const {
    int id;
    std::vector<std::string> matches;
    if( match( name, matches, id ) ) {
      
      util::Config grid;
      grid.set("type", type());

      if( id == 0 ) {
        grid.set( "nx", to_int(matches[0]) );
        grid.set( "ny", to_int(matches[1]) );
        return create( grid );
      }

      if( id == 1 ) {
        grid.set( "N", to_int(matches[0]) );
        return create( grid );
      }

    }
    return nullptr;
  }
  
  virtual const Grid::grid_t* create( const Grid::Config& config ) const {
    return create_lonlat( config, Shift(true,true) );
  }

} shifted_lonlat_;

//---------------------------------------------------------------------------------------------------------------------

static class shifted_lon : public GridCreator {

public:

  shifted_lon(): GridCreator( "shifted_lon", {
    "^[Ss][Ll][Oo][Nn]([0-9]+)x([0-9]+)$",
    "^[Ss][Ll][Oo][Nn]([0-9]+)$"           } ){}

  virtual void print(std::ostream& os) const {
    os << std::left << std::setw(20) << "Slon<nx>x<ny> / Slon<gauss>" << "Shifted longitude grid";
  }

  virtual const Grid::grid_t* create( const std::string& name ) const {
    int id;
    std::vector<std::string> matches;
    if( match( name, matches, id ) ) {
      
      util::Config grid;
      grid.set("type", type());

      if( id == 0 ) {
        grid.set( "nx", to_int(matches[0]) );
        grid.set( "ny", to_int(matches[1]) );
        return create( grid );
      }

      if( id == 1 ) {
        grid.set( "N", to_int(matches[0]) );
        return create( grid );
      }

    }
    return nullptr;
  }
  
  virtual const Grid::grid_t* create( const Grid::Config& config ) const {
    return create_lonlat( config, Shift(true,false) );
  }

} shifted_lon_;


//---------------------------------------------------------------------------------------------------------------------

static class shifted_lat : public GridCreator {

public:

  shifted_lat(): GridCreator( "shifted_lat", {
    "^[Ss][Ll][Aa][Tt]([0-9]+)x([0-9]+)$",
    "^[Ss][Ll][Aa][Tt]([0-9]+)$"           } ){}

  virtual void print(std::ostream& os) const {
    os << std::left << std::setw(20) << "Slat<nx>x<ny> / Slat<gauss>" << "Shifted latitude grid";    
  }

  virtual const Grid::grid_t* create( const std::string& name ) const {
    int id;
    std::vector<std::string> matches;
    if( match( name, matches, id ) ) {
      
      util::Config grid;
      grid.set("type", type());

      if( id == 0 ) {
        grid.set( "nx", to_int(matches[0]) );
        grid.set( "ny", to_int(matches[1]) );
        return create( grid );
      }

      if( id == 1 ) {
        grid.set( "N", to_int(matches[0]) );
        return create( grid );
      }

    }
    return nullptr;
  }
  
  virtual const Grid::grid_t* create( const Grid::Config& config ) const {
    return create_lonlat( config, Shift(false,true) );
  }

} shifted_lat_;

//---------------------------------------------------------------------------------------------------------------------

}  // anonymous namespace
}  // namespace grid
}  // namespace atlas
