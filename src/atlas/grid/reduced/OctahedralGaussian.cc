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

std::string OctahedralGaussian::shortName() const {
   std::ostringstream s;
   s << "O" << nlat()/2;
   return s.str();
}

std::vector<long> OctahedralGaussian::computePL(const size_t N) {
    const size_t start = 20;
    std::vector<long> pl(N);
    for(size_t jlat=0; jlat < N; ++jlat) {
        pl[jlat] = start + 4*jlat;
    }
    return pl;
}

void OctahedralGaussian::setup(size_t N) {

		util::Config config_spacing, config_domain, config_proj;
		
		// number of longitudes per latitude
		std::vector<long>   pl(N);
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

		// projection is lonlat
		util::Config config_proj;
		config_proj.set("projectionType","lonlat");
		projection_=projection::Projection::create(config_proj);
		
		// setup
		setup(N);
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
}  // namespace atlas

