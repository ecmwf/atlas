#include "atlas/grid/detail/grid/reduced/ReducedGaussian.h"

#include "atlas/grid/detail/spacing/Spacing.h"

namespace atlas {
namespace grid {
namespace detail {
namespace grid {
namespace reduced {

register_BuilderT1(Grid,ReducedGaussian,ReducedGaussian::grid_type_str());

std::string ReducedGaussian::grid_type_str() {
    return "reduced_gaussian";
}


std::string ReducedGaussian::className() {
    return "atlas.grid.reduced.ReducedGaussian";
}

std::string ReducedGaussian::shortName() const {
  return "reduced_gaussian";
}

void ReducedGaussian::setup(const size_t N, const long pl[]) {

    // configs for spacing, domain and projection
    util::Config config_spacing, config_domain, config_proj;

    // number of latitudes
    size_t ny=2*N;

    // domain is global
    config_domain.set("type","global");
    domain_ = Domain(config_domain);


    // mirror pl around equator
    std::vector<long> nx(ny);
    for (int jlat=0;jlat<N;jlat++) {
      nx[jlat]=pl[jlat];
      nx[ny-1-jlat]=pl[jlat];
    }

    // determine input for Structured::setup
    std::vector<double> xmin(ny);    // first longitude per latitude
    std::vector<double> xmax(ny);    // last longitude per latitude
    std::vector<double> dx(ny);

    // latitudes: gaussian spacing
    config_spacing.set("type","gaussian");
    config_spacing.set("start", 90.0);
    config_spacing.set("end",  -90.0);
    config_spacing.set("N",ny);
    YSpace yspace(config_spacing);

    // loop over latitudes to set bounds
    for (int j=0;j<ny;j++) {
      xmin[j] = 0.0;
      xmax[j] = 360.;
      dx[j]   = 360.0/double(nx[j]);
    }

    // setup Structured grid
    Structured::setup(yspace, nx, xmin, xmax, dx);
    Structured::N_=N;
}

ReducedGaussian::ReducedGaussian() :
    Structured() {
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

     // default projection is lonlat
     util::Config config_proj;
     if( not config.get("projection",config_proj) ) {
       config_proj.set("type","lonlat");
     }
     projection_ = Projection(config_proj);

    // setup
    setup(N,pl.data());
}

ReducedGaussian::ReducedGaussian(const int N, const long pl[]) {
    // projection is lonlat
    util::Config config_proj;
    config_proj.set("type","lonlat");
    projection_ = Projection(config_proj);

    // setup
    setup(N,pl);

}

Grid::Spec ReducedGaussian::spec() const {
    Grid::Spec grid_spec;

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
}  // namespace detail
}  // namespace grid
}  // namespace atlas

