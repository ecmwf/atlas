#include "atlas/grid/detail/grid/reduced/ReducedLonLat.h"

#include "atlas/grid/detail/spacing/Spacing.h"

namespace atlas {
namespace grid {
namespace detail {
namespace grid {
namespace reduced {

register_BuilderT1(Grid,ReducedLonLat,ReducedLonLat::grid_type_str());

std::string ReducedLonLat::grid_type_str() {
    return "reduced_lonlat";
}

std::string ReducedLonLat::className() {
    return "atlas.grid.reduced.ReducedLonLat";
}

std::string ReducedLonLat::shortName() const {
  return "reduced_lonlat";
}


void ReducedLonLat::setup(size_t ny, long pl[]) {

    util::Config config_spacing, config_domain, config_proj;

    // projection is lonlat
    config_proj.set("type","lonlat");
    projection_ = projection::Projection::create(config_proj);

    // domain is global
    config_domain.set("type","global");
    domain_ = Domain(config_domain);

    // determine input for Structured::setup
    std::vector<double> xmin(ny);    // first longitude per latitude
    std::vector<double> xmax(ny);    // last longitude per latitude
    std::vector<double> dx(ny);      // step per latitude
    std::vector<long> nx(ny);      // step per latitude

    // latitudes: linear spacing
    config_spacing.set("type","linear");
    config_spacing.set("start",90.0);
    config_spacing.set("end",-90.0);
    config_spacing.set("N",ny);
    YSpace yspace(config_spacing);

    // loop over latitudes to set bounds
    for (int j=0;j<ny;j++) {
      nx[j]=pl[j];
      xmin[j]=0.;
      xmax[j]=360.;
      dx[j]=360./double(nx[j]);
    }

    // setup Structured grid
    Structured::setup(yspace, nx, xmin, xmax, dx);
    if (ny%2) {
      Structured::N_=0;    // odd number of latitudes
    } else {
      Structured::N_=ny/2;
    }

}

ReducedLonLat::ReducedLonLat(const util::Config& config) :
    Structured()
{
    size_t N, nlat;
    if( ! config.get("nlat",nlat) ) {
      if ( !config.get("N",N) ) {
        throw eckit::BadParameter("N or nlat missing in config",Here());
      } else {
        nlat=2*N;
      }
    }

    std::vector<long> pl;
    if( ! config.get("pl",pl) ) throw eckit::BadParameter("pl missing in config",Here());

    // mirror if necessary
    if ( pl.size() == N ) {
      // mirror to get length 2*N
      pl.resize(nlat);
      for (int jlat=0;jlat<N;jlat++) {
        pl[nlat-1-jlat]=pl[jlat];
      }
    }

    // check length
     if ( pl.size() != nlat )
       throw eckit::BadParameter("pl should have length N or 2*N",Here());

    // setup
    setup(N,pl.data());
}

Grid::Spec ReducedLonLat::spec() const {
    Grid::Spec grid_spec;

    // general specs
    grid_spec=Grid::spec();

    // specs for reduced gaussian
    grid_spec.set("pl", eckit::makeVectorValue(pl()));

    return grid_spec;
}


}  // namespace regular
}  // namespace grid
}  // namespace detail
}  // namespace grid
}  // namespace atlas

