#include "atlas/grid/reduced/OctahedralGaussian.h"

#include "atlas/grid/spacing/Spacing.h"

namespace atlas {
namespace grid {
namespace reduced {

register_BuilderT1(Grid,OctahedralGaussian,OctahedralGaussian::grid_type_str());

std::string OctahedralGaussian::grid_type_str() {
    return "octahedralGaussian";
}

std::string OctahedralGaussian::className() {
    return "atlas.grid.reduced.OctahedralGaussian";
}

std::vector<long> OctahedralGaussian::computePL(const size_t N) {
    const size_t start = 20;
    std::vector<long> pl(2*N);
    for(size_t jlat=0; jlat < N; ++jlat) {
        pl[jlat] = start + 4*jlat;
        pl[2*N-1-jlat] = pl[jlat];
    }
    return pl;
}

void OctahedralGaussian::setup(size_t N) {

		util::Config config_spacing, config_domain, config_proj;
		
		// number of longitudes per latitude
		size_t ny=2*N;
		std::vector<long>   pl(ny);
    pl=computePL(N);
    
    // setup from reducedGaussian
		ReducedGaussian::setup(N,pl.data());
		
}

OctahedralGaussian::OctahedralGaussian(const util::Config& config) :
    ReducedGaussian()
{
		size_t N;
		
		// get N from config
		if ( !config.get("N",N) ) throw eckit::BadParameter("OctahedralGaussian requires N",Here());
		
		// setup
		setup(N);
}

eckit::Properties OctahedralGaussian::spec() const {
    eckit::Properties grid_spec;
    grid_spec.set("grid_type",  gridType());
    /*
    grid_spec.set("short_name", shortName());
    grid_spec.set("N",    N());
    grid_spec.set("nlon", nlon());
    grid_spec.set("nlat", nlat());
    grid_spec.set("domain", domain_spec(domain_) );
    */
    return grid_spec;
}


}  // namespace regular
}  // namespace grid
}  // namespace atlas

