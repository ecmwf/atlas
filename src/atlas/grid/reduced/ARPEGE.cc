#include "atlas/grid/reduced/ARPEGE.h"

#include "atlas/grid/projection/Projection.h"

namespace atlas {
namespace grid {
namespace reduced {

register_BuilderT1(Grid,ARPEGE,ARPEGE::grid_type_str());

std::string ARPEGE::grid_type_str() {
    return "ARPEGE";
}


std::string ARPEGE::className() {
    return "atlas.grid.reduced.ARPEGE";
}


void ARPEGE::setup(size_t N) {

		util::Config config_spacing, config_domain, config_proj;

}

ARPEGE::ARPEGE(const util::Config& config) :
    ClassicGaussian()
{
		size_t N;
		double c;
		util::Config config_proj;
		
		// get N from config
		if ( !config.get("N",N) ) throw eckit::BadParameter("ARPEGE requires N",Here());
		
		// set projection
		// get stretching factor; defaults to 1
		if ( !config.get("stretching_factor",c) ) c=1.0;
		config_proj.set("stretching_factor",c);
		std::vector<double> pole(2);
		if ( config.get("pole",pole) ) {
			config_proj.set("projectionType","rotatedSchmidt");
			config_proj.set("pole",pole);
		} else {
			config_proj.set("projectionType","schmidt");
		}
		projection_=projection::Projection::create(config_proj);

		// setup
		ClassicGaussian::setup(N);
}

eckit::Properties ARPEGE::spec() const {
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

