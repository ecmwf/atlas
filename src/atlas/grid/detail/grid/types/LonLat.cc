#include "LonLat.h"

#include "eckit/utils/Translator.h"

#include "atlas/grid/detail/grid/regular/RegularLonLat.h"
#include "atlas/grid/detail/grid/regular/ShiftedLonLat.h"
#include "atlas/grid/detail/grid/regular/ShiftedLon.h"
#include "atlas/grid/detail/grid/regular/ShiftedLat.h"

namespace atlas {
namespace grid {
namespace { // anonymous

static eckit::Translator<std::string,int> to_int;

//---------------------------------------------------------------------------------------------------------------------

static class regular_lonlat : public GridCreator {

public:

  regular_lonlat() : GridCreator("^[Ll]([0-9]+)$") {}

  virtual void print(std::ostream& os) const {
    os << std::left << std::setw(20) << "L<gauss>" << "Regular longitude-latitude grid";
  }

  virtual const Grid::grid_t* create( const std::string& name ) const {
    std::vector<std::string> matches;
    if( match( name, matches ) ) {
      util::Config grid;
      int N = to_int(matches[0]);
      grid.set("type", detail::grid::regular::RegularLonLat::static_type());
      grid.set("N",N);
      return GridCreator::create( grid );
    }
    return nullptr;
  }

} regular_lonlat_;

//---------------------------------------------------------------------------------------------------------------------

static class shifted_lonlat : public GridCreator {

public:

  shifted_lonlat() : GridCreator("^[Ss]([0-9]+)$") {}

  virtual void print(std::ostream& os) const {
    os << std::left << std::setw(20) << "S<gauss>" << "Shifted longitude-latitude grid";
  }

  virtual const Grid::grid_t* create( const std::string& name ) const {
    std::vector<std::string> matches;
    if( match( name, matches ) ) {
      util::Config grid;
      int N = to_int(matches[0]);
      grid.set("type", detail::grid::regular::ShiftedLonLat::static_type());
      grid.set("N",N);
      return GridCreator::create( grid );
    }
    return nullptr;
  }

} shifted_lonlat_;

//---------------------------------------------------------------------------------------------------------------------

static class shifted_lon : public GridCreator {

public:

  shifted_lon() : GridCreator("^[Ss][Ll][Oo][Nn]([0-9]+)$") {}

  virtual void print(std::ostream& os) const {
    os << std::left << std::setw(20) << "Slon<gauss>" << "Shifted longitude grid";
  }

  virtual const Grid::grid_t* create( const std::string& name ) const {
    std::vector<std::string> matches;
    if( match( name, matches ) ) {
      util::Config grid;
      int N = to_int(matches[0]);
      grid.set("type", detail::grid::regular::ShiftedLon::static_type());
      grid.set("N",N);
      return GridCreator::create( grid );
    }
    return nullptr;
  }

} shifted_lon_;


//---------------------------------------------------------------------------------------------------------------------

static class shifted_lat : public GridCreator {

public:

  shifted_lat() : GridCreator("^[Ss][Ll][Aa][Tt]([0-9]+)$") {}

  virtual void print(std::ostream& os) const {
    os << std::left << std::setw(20) << "Slat<gauss>" << "Shifted latitude grid";
  }

  virtual const Grid::grid_t* create( const std::string& name ) const {
    std::vector<std::string> matches;
    if( match( name, matches ) ) {
      util::Config grid;
      int N = to_int(matches[0]);
      grid.set("type", detail::grid::regular::ShiftedLat::static_type());
      grid.set("N",N);
      return GridCreator::create( grid );
    }
    return nullptr;
  }

} shifted_lat_;

//---------------------------------------------------------------------------------------------------------------------

static class regular_lonlat_x : public GridCreator {

public:

  regular_lonlat_x() : GridCreator("^[Ll]([0-9]+)x([0-9]+)$") {}

  virtual void print(std::ostream& os) const {
    os << std::left << std::setw(20) << "L<nx>x<ny>" << "Regular longitude-latitude grid";
  }


  virtual const Grid::grid_t* create( const std::string& name ) const {
    std::vector<std::string> matches;
    if( match( name, matches ) ) {
      util::Config grid;
      int nx = to_int(matches[0]);
      int ny = to_int(matches[1]);
      grid.set("type", detail::grid::regular::RegularLonLat::static_type());
      grid.set("nx",nx);
      grid.set("ny",ny);
      return GridCreator::create( grid );
    }
    return nullptr;
  }

} regular_lonlat_x_;

//---------------------------------------------------------------------------------------------------------------------

static class shifted_lonlat_x : public GridCreator {

public:

  shifted_lonlat_x() : GridCreator("^[Ss]([0-9]+)x([0-9]+)$") {}

  virtual void print(std::ostream& os) const {
    os << std::left << std::setw(20) << "S<nx>x<ny>" << "Shifted longitude-latitude grid";
  }

  virtual const Grid::grid_t* create( const std::string& name ) const {
    std::vector<std::string> matches;
    if( match( name, matches ) ) {
      util::Config grid;
      int nx = to_int(matches[0]);
      int ny = to_int(matches[1]);
      grid.set("type", detail::grid::regular::ShiftedLonLat::static_type());
      grid.set("nx",nx);
      grid.set("ny",ny);
      return GridCreator::create( grid );
    }
    return nullptr;
  }

} shifted_lonlat_x_;


//---------------------------------------------------------------------------------------------------------------------

static class shifted_lon_x : public GridCreator {

public:

  shifted_lon_x() : GridCreator("^[Ss][Ll][Oo][Nn]([0-9]+)x([0-9]+)$") {}

  virtual void print(std::ostream& os) const {
    os << std::left << std::setw(20) << "Slon<nx>x<ny>" << "Shifted longitude grid";
  }

  virtual const Grid::grid_t* create( const std::string& name ) const {
    std::vector<std::string> matches;
    if( match( name, matches ) ) {
      util::Config grid;
      int nx = to_int(matches[0]);
      int ny = to_int(matches[1]);
      grid.set("type", detail::grid::regular::ShiftedLon::static_type());
      grid.set("nx",nx);
      grid.set("ny",ny);
      return GridCreator::create( grid );
    }
    return nullptr;
  }

} shifted_lon_x_;


//---------------------------------------------------------------------------------------------------------------------

static class shifted_lat_x : public GridCreator {

public:

  shifted_lat_x() : GridCreator("^[Ss][Ll][Aa][Tt]([0-9]+)x([0-9]+)$") {}

  virtual void print(std::ostream& os) const {
    os << std::left << std::setw(20) << "Slat<nx>x<ny>" << "Shifted latitude grid";
  }

  virtual const Grid::grid_t* create( const std::string& name ) const {
    std::vector<std::string> matches;
    if( match( name, matches ) ) {
      util::Config grid;
      int nx = to_int(matches[0]);
      int ny = to_int(matches[1]);
      grid.set("type", detail::grid::regular::ShiftedLat::static_type());
      grid.set("nx",nx);
      grid.set("ny",ny);
      return GridCreator::create( grid );
    }
    return nullptr;
  }

} shifted_lat_x_;

}  // anonymous namespace
}  // namespace grid
}  // namespace atlas
