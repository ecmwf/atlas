#include "atlas/grid/reduced/ReducedGaussian.h"

#include "atlas/grid/spacing/Spacing.h"

namespace atlas {
namespace grid {
namespace reduced {

register_BuilderT1(Grid,ReducedGaussian,ReducedGaussian::grid_type_str());

std::string ReducedGaussian::grid_type_str() {
    return "reducedGaussian";
}


std::string ReducedGaussian::className() {
    return "atlas.grid.reduced.ReducedGaussian";
}

std::string ReducedGaussian::shortName() const {
  return "reducedGaussian";
}

void ReducedGaussian::setup(const size_t N, const long pl[]) {

    // configs for spacing, domain and projection
    util::Config config_spacing, config_domain, config_proj;

    // number of latitudes
    size_t ny=2*N;

    // domain is global
    config_domain.set("domainType","global");
    domain_.reset( domain::Domain::create(config_domain) );


    // mirror pl around equator
    std::vector<long> pll(ny);
    for (int jlat=0;jlat<N;jlat++) {
      pll[jlat]=pl[jlat];
      pll[ny-1-jlat]=pl[jlat];
    }

    // determine input for Structured::setup
    std::vector<double> xmin(ny);    // first longitude per latitude
    std::vector<double> xmax(ny);    // last longitude per latitude
    std::vector<double> y(ny);      // latitudes

    // latitudes: gaussian spacing
    config_spacing.set("spacingType","gaussian");
    config_spacing.set("xmin",90.0);
    config_spacing.set("xmax",-90.0);
    config_spacing.set("N",ny);
    eckit::SharedPtr<spacing::Spacing> spacing_y ( spacing::Spacing::create(config_spacing) );
    spacing_y->generate(y);

    // loop over latitudes to set bounds
    xmin.assign(ny,0.0);
    for (int jlat=0;jlat<ny;jlat++) {
      xmin[jlat]=0.0;
      xmax[jlat]=(pll[jlat]-1)*360.0/pll[jlat];
    }

    // setup Structured grid
    Structured::setup(ny,y.data(), pll.data(), xmin.data(), xmax.data());
    Structured::N_=N;

}

ReducedGaussian::ReducedGaussian(const util::Config& config) :
    Structured()
{
    size_t N;
    if( ! config.has("N") ) throw eckit::BadParameter("N missing in config",Here());
    config.get("N",N);

    std::vector<long> pl;
    if( ! config.has("pl") ) throw eckit::BadParameter("pl missing in config",Here());
    config.get("pl",pl);

    // check length
     if ( pl.size() != N  )
       throw eckit::BadParameter("pl should have length N",Here());

    // projection is lonlat
    util::Config config_proj;
    config_proj.set("projectionType","lonlat");
    projection_.reset( projection::Projection::create(config_proj) );

    // setup
    setup(N,pl.data());
}

ReducedGaussian::ReducedGaussian(const int N, const long pl[]) {
    // projection is lonlat
    util::Config config_proj;
    config_proj.set("projectionType","lonlat");
    projection_.reset( projection::Projection::create(config_proj) );

    // setup
    setup(N,pl);

}

eckit::Properties ReducedGaussian::spec() const {
    eckit::Properties grid_spec;

    // general specs
    grid_spec=Grid::spec();

    // specs for reduced gaussian
    grid_spec.set("pl", eckit::makeVectorValue(pl()));

    return grid_spec;
}


extern "C" {


    Structured* atlas__grid__reduced__ReducedGaussian_int(size_t N, int pl[]) {
        std::vector<long> pl_vector;
        pl_vector.assign(pl,pl+N);
        return new ReducedGaussian(N,pl_vector.data());
    }


    Structured* atlas__grid__reduced__ReducedGaussian_long(size_t N, long pl[]) {
        return new ReducedGaussian(N,pl);
    }


}



}  // namespace regular
}  // namespace grid
}  // namespace atlas

