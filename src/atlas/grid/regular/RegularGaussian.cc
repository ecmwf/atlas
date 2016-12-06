#include "atlas/grid/regular/RegularGaussian.h"


namespace atlas {
namespace grid {
namespace regular {

register_BuilderT1(Grid,RegularGaussian,RegularGaussian::grid_type_str());

std::string RegularGaussian::grid_type_str() {
    return "regularGaussian";
}


std::string RegularGaussian::className() {
    return "atlas.grid.regular.RegularGaussian";
}


void RegularGaussian::setup() {

		// setup regular grid
    Regular::setup();
    
}

RegularGaussian::RegularGaussian(const util::Config& config) :
    Regular()
{
		long nlon, nlat, N;
		util::Config config_proj, config_spacing, config_domain;
		
		// projection is lonlat
		config_proj.set("projectionType","lonlat");
		projection_=projection::Projection::create(config_proj);
		
		// dimensions
		if ( config.get("N",N) ) {
			nlon=4*N;nlat=2*N;
		} else {
			throw eckit::BadParameter("RegularGaussian requires parameter N",Here());
		}

		// spacing is uniform in x and y
		config_spacing.set("spacingType","uniform");
		config_spacing.set("xmin",0.0);
		config_spacing.set("xmax",(nlon-1)*360.0/nlon);
		config_spacing.set("N",nlon);
		spacing_x_=spacing::Spacing::create(config_spacing);

		config_spacing.set("spacingType","gaussian");
		config_spacing.set("xmin",90.0);
		config_spacing.set("xmax",-90.0);
		config_spacing.set("N",nlat);
		spacing_y_=spacing::Spacing::create(config_spacing);
		
		// domain is global
		config_domain.set("domainType","global");
		domain_=domain::Domain::create(config_domain);
		
		// perform setup
		setup();
}

eckit::Properties RegularGaussian::spec() const {
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

