#include "atlas/grid/regular/RegularGlobalLonLat.h"


namespace atlas {
namespace grid {
namespace regular {

register_BuilderT1(Grid,RegularGlobalLonLat,RegularGlobalLonLat::grid_type_str());

std::string RegularGlobalLonLat::grid_type_str() {
    return "regularGlobalLonLat";
}


std::string RegularGlobalLonLat::className() {
    return "atlas.grid.regular.RegularGlobalLonLat";
}


void RegularGlobalLonLat::setup() {

		// setup regular grid
    Regular::setup();
    
}

RegularGlobalLonLat::RegularGlobalLonLat(const util::Config& config) :
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
			if ( !config.get("nlon",nlon) || !config.get("nlat",nlat) ) {
				throw eckit::BadParameter("RegularGlobalLonLat requires either N, or (nlon,nlat)",Here());
			}
		}

		// spacing is uniform in x and y
		config_spacing.set("spacingType","uniform");
		config_spacing.set("xmin",0.0);
		config_spacing.set("xmax",(nlon-1)*360.0/nlon);
		config_spacing.set("N",nlon);
		spacing_x_=spacing::Spacing::create(config_spacing);
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

eckit::Properties RegularGlobalLonLat::spec() const {
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

