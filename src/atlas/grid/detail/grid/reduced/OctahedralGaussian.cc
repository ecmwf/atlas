#include "atlas/grid/detail/grid/reduced/OctahedralGaussian.h"

#include "atlas/grid/detail/spacing/Spacing.h"

namespace atlas {
namespace grid {
namespace detail {
namespace grid {
namespace reduced {

register_BuilderT1(Grid,OctahedralGaussian,OctahedralGaussian::grid_type_str());

std::string OctahedralGaussian::grid_type_str() {
    return "octahedral_gaussian";
}

std::string OctahedralGaussian::className() {
    return "atlas.grid.reduced.OctahedralGaussian";
}

std::string OctahedralGaussian::shortName() const {
   std::ostringstream s;
   s << "O" << nlat()/2;
   return s.str();
}

std::vector<long> OctahedralGaussian::computePL(const size_t N, const size_t start) {
    std::vector<long> pl(N);
    for(size_t jlat=0; jlat < N; ++jlat) {
        pl[jlat] = start + 4*jlat;
    }
    return pl;
}

void OctahedralGaussian::setup(const size_t N, const size_t start) {

    // number of longitudes per latitude
    std::vector<long>   pl(N);
    pl=computePL(N,start);

    // setup from reducedGaussian
    ReducedGaussian::setup(N,pl.data());
}

OctahedralGaussian::OctahedralGaussian(const util::Config& config) :
    ReducedGaussian()
{
    size_t N;

    // get N from config
    if ( !config.get("N",N) ) throw eckit::BadParameter("OctahedralGaussian requires N",Here());

    // default projection is lonlat
    util::Config config_proj;
    if( not config.get("projection",config_proj) ) {
      config_proj.set("type","lonlat");
    }
    projection_ = Projection(config_proj);

    size_t start=20;
    config.get("nx[0]",start);

    // setup
    setup(N,start);
}

eckit::Properties OctahedralGaussian::spec() const {
    eckit::Properties grid_spec;

    // general specs
    grid_spec=Grid::spec();

    // specs for this grid are in shortname

    return grid_spec;
}


}  // namespace regular
}  // namespace grid
}  // namespace detail
}  // namespace grid
}  // namespace atlas

