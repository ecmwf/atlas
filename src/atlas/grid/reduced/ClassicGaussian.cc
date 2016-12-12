#include "atlas/grid/reduced/ClassicGaussian.h"

#include "atlas/grid/spacing/Spacing.h"
#include "atlas/grid/reduced/classic/PointsPerLatitude.h"

namespace atlas {
namespace grid {
namespace reduced {

register_BuilderT1(Grid,ClassicGaussian,ClassicGaussian::grid_type_str());

std::string ClassicGaussian::grid_type_str() {
    return "classicGaussian";
}


std::string ClassicGaussian::className() {
    return "atlas.grid.reduced.ClassicGaussian";
}


void ClassicGaussian::setup(size_t N) {

		util::Config config_spacing, config_domain, config_proj;

		// number of longitudes: from predefined sets
		size_t ny=2*N;
		std::vector<long> pl(ny);			// number of longitudes per latitude
		classic::points_per_latitude_npole_spole(N,pl.data());
    
    // setup from reducedGaussian
    ReducedGaussian::setup(N,pl.data());
    
}

ClassicGaussian::ClassicGaussian(const util::Config& config) :
    ReducedGaussian()
{
		size_t N;
		
		// get N from config
		if ( !config.get("N",N) ) throw eckit::BadParameter("ClassicGaussian requires N",Here());

		// projection is lonlat
		util::Config config_proj;
		config_proj.set("projectionType","lonlat");
		projection_=projection::Projection::create(config_proj);
		
		// setup
		setup(N);
}

eckit::Properties ClassicGaussian::spec() const {
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

