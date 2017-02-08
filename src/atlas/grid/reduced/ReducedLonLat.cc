#include "atlas/grid/reduced/ReducedLonLat.h"

#include "atlas/grid/spacing/Spacing.h"

namespace atlas {
namespace grid {
namespace reduced {

register_BuilderT1(Grid,ReducedLonLat,ReducedLonLat::grid_type_str());

std::string ReducedLonLat::grid_type_str() {
    return "reducedLonLat";
}

std::string ReducedLonLat::className() {
    return "atlas.grid.reduced.ReducedLonLat";
}

std::string ReducedLonLat::shortName() const {
  return "reducedLonLat";
}


void ReducedLonLat::setup(size_t ny, long pl[]) {

    util::Config config_spacing, config_domain, config_proj;

    // projection is lonlat
    config_proj.set("projectionType","lonlat");
    projection_.reset( projection::Projection::create(config_proj) );

    // domain is global
    config_domain.set("domainType","global");
    domain_.reset( domain::Domain::create(config_domain) );

    // determine input for Structured::setup
    std::vector<double> xmin(ny);    // first longitude per latitude
    std::vector<double> xmax(ny);    // last longitude per latitude
    std::vector<double> y(ny);      // latitudes

    // latitudes: gaussian spacing
    config_spacing.set("spacingType","uniform");
    config_spacing.set("xmin",90.0);
    config_spacing.set("xmax",-90.0);
    config_spacing.set("N",ny);
    eckit::SharedPtr<spacing::Spacing> spacing_y ( spacing::Spacing::create(config_spacing) );
    spacing_y->generate(y);

    // loop over latitudes to set bounds
    xmin.assign(ny,0.0);
    for (int jlat=0;jlat<ny;jlat++) {
      xmax[jlat]=(pl[jlat]-1)*360.0/pl[jlat];
    }

    // setup Structured grid
    Structured::setup(ny,y.data(), pl, xmin.data(), xmax.data());
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

eckit::Properties ReducedLonLat::spec() const {
    eckit::Properties grid_spec;

    // general specs
    grid_spec=Grid::spec();

    // specs for reduced gaussian
    grid_spec.set("pl", eckit::makeVectorValue(pl()));

    return grid_spec;
}


}  // namespace regular
}  // namespace grid
}  // namespace atlas

