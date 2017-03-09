#include "Gaussian.h"

#include <cmath>
#include <limits>

#include "eckit/utils/Translator.h"

#include "atlas/grid/detail/pl/classic_gaussian/PointsPerLatitude.h"
#include "atlas/grid/detail/grid/GridBuilder.h"

namespace atlas {
namespace grid {
namespace { // anonymous

static eckit::Translator<std::string,int> to_int;

static Projection projection( const Grid::Config& grid ) {

    Grid::Config config;
    if( grid.get("projection",config) ) {
      return Projection(config);
    }
    return Projection();
}

static Domain domain( const Grid::Config& grid ) {

    Grid::Config config;
    if( grid.get("domain",config) ) {
      return Domain(config);
    }
    return Domain();
}

static Spacing yspace( const Grid::Config& grid ) {

    long N;
    grid.get("N",N);

    Grid::Config config;
    config.set("type","gaussian");
    config.set("start", 90.0);
    config.set("end",  -90.0);
    config.set("N",2*N);
    return Spacing(config);

}

static StructuredGrid::grid_t::XSpace* xspace( const std::vector<long>& nx ) {
  StructuredGrid::grid_t::XSpace *_xspace = new StructuredGrid::grid_t::XSpace(nx.size());
  _xspace->nx = nx;
  _xspace->nxmax = 0;
  for( size_t j=0; j<_xspace->ny; ++j ) {
    _xspace->xmin.push_back( 0. );
    _xspace->xmax.push_back( 360. );
    _xspace->dx.push_back( 360./_xspace->nx[j] );
    _xspace->nxmin = std::min( _xspace->nxmin, size_t(_xspace->nx[j]) );
    _xspace->nxmax = std::max( _xspace->nxmax, size_t(_xspace->nx[j]) );
  }
  return _xspace;
}

//---------------------------------------------------------------------------------------------------------------------

static class classic_gaussian : public GridBuilder {

public:

  classic_gaussian() : GridBuilder("classic_gaussian",{"^[Nn]([0-9]+)$"}) {}

  virtual void print(std::ostream& os) const {
    os << std::left << std::setw(20) << "N<gauss>" << "Classic Gaussian grid";
  }

  virtual const Grid::grid_t* create( const std::string& name, const Grid::Config& config ) const {
    int id;
    std::vector<std::string> matches;
    if( match( name, matches, id ) ) {
      util::Config grid(config);
      int N = to_int(matches[0]);
      grid.set("type", type() );
      grid.set("N",N);
      return create( grid );
    }
    return nullptr;
  }

  virtual const Grid::grid_t* create( const Grid::Config& config ) const {
    long N;
    config.get("N",N);
    std::vector<long> nx(2*N);
    detail::pl::classic_gaussian::points_per_latitude_npole_spole( N, nx.data() );
    return new StructuredGrid::grid_t( "N"+std::to_string(N), projection(config), xspace(nx) , yspace(config), domain(config) );
  }

} classic_gaussian_;

//---------------------------------------------------------------------------------------------------------------------

static class octahedral_gaussian : GridBuilder {

public:

  octahedral_gaussian(): GridBuilder("octahedral_gaussian",{"^[Oo]([0-9]+)$"}) {}

  virtual void print(std::ostream& os) const {
    os << std::left << std::setw(20) << "O<gauss>" << "Octahedral Gaussian grid";
  }

  virtual const Grid::grid_t* create( const std::string& name, const Grid::Config& config ) const {
    int id;
    std::vector<std::string> matches;
    if( match( name, matches, id ) ) {
      util::Config grid(config);
      int N = to_int(matches[0]);
      grid.set("type", type());
      grid.set("N",N);
      return create( grid );
    }
    return nullptr;
  }

  virtual const Grid::grid_t* create( const Grid::Config& config ) const {
    long N;
    config.get("N",N);

    long start = 20;
    config.get("nx[0]",start);

    std::vector<long> nx(2*N);
    for(long j=0; j<N; ++j) {
        nx[j] = start + 4*j;
        nx[2*N-1-j] = nx[j];
    }
    return new StructuredGrid::grid_t( "O"+std::to_string(N), projection(config), xspace(nx) , yspace(config), domain(config) );
  }

} octahedral_gaussian_;

//---------------------------------------------------------------------------------------------------------------------

static class regular_gaussian : GridBuilder {

public:

  regular_gaussian(): GridBuilder("regular_gaussian",{"^[Ff]([0-9]+)$"}) {}

  virtual void print(std::ostream& os) const {
    os << std::left << std::setw(20) << "F<gauss>" << "Regular Gaussian grid";
  }

  virtual const Grid::grid_t* create( const std::string& name, const Grid::Config& config ) const {
    int id;
    std::vector<std::string> matches;
    if( match( name, matches, id ) ) {
      util::Config grid(config);
      int N = to_int(matches[0]);
      grid.set("type", type());
      grid.set("N",N);
      return create( grid );
    }
    return nullptr;
  }

  virtual const Grid::grid_t* create( const Grid::Config& config) const {
    long N;
    config.get("N",N);
    std::vector<long> nx(2*N,4*N);
    return new StructuredGrid::grid_t( "F"+std::to_string(N), projection(config), xspace(nx) , yspace(config), domain(config) );
  }


} regular_gaussian_;


//---------------------------------------------------------------------------------------------------------------------

} // anonymous namespace

namespace detail {
namespace grid {

StructuredGrid::grid_t* reduced_gaussian( const std::vector<long>& nx ) {

  Grid::Config yspace;
  yspace.set("type","gaussian");
  yspace.set("start", 90.0);
  yspace.set("end",  -90.0);
  yspace.set("N",nx.size());

  return new StructuredGrid::grid_t( Projection(), xspace(nx) , Spacing(yspace), Domain() );
}


extern "C" {


    StructuredGrid::grid_t* atlas__grid__reduced__ReducedGaussian_int(long N, int half_pl[]) {
      std::vector<long> nx(2*N);
      for( long j=0; j<N; ++j ) {
        nx[j] = half_pl[j];
        nx[2*N-1-j] = nx[j];
      }
      return reduced_gaussian(nx);
    }


    StructuredGrid::grid_t* atlas__grid__reduced__ReducedGaussian_long(long N, long half_pl[]) {
        std::vector<long> nx(2*N);
        for( long j=0; j<N; ++j ) {
          nx[j] = half_pl[j];
          nx[2*N-1-j] = nx[j];
        }
        return reduced_gaussian(nx);
    }

}


} // namespace grid
} // namespace detail

} // namespace grid
} // namespace atlas
