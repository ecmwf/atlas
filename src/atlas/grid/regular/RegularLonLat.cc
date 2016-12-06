#include "atlas/grid/regular/RegularLonLat.h"


namespace atlas {
namespace grid {
namespace regular {

register_BuilderT1(Grid,RegularLonLat,RegularLonLat::grid_type_str());

std::string RegularLonLat::grid_type_str() {
    return "regularLonLat";
}


std::string RegularLonLat::className() {
    return "atlas.grid.regular.RegularLonLat";
}


void RegularLonLat::setup(long nlon, long nlat) {

		util::Config config_proj, config_spacing, config_domain;

		// projection is lonlat
		config_proj.set("projectionType","lonlat");
		projection_=projection::Projection::create(config_proj);
		
		// spacing is uniform in x
		config_spacing.set("spacingType","uniform");
		config_spacing.set("xmin",(shiftLon_ ? 0.5 : 0.0)*360.0/nlon );
		config_spacing.set("xmax",(shiftLon_ ? nlon-0.5 : nlon-1)*360.0/nlon );
		config_spacing.set("N",nlon);
		spacing_x_=spacing::Spacing::create(config_spacing);
		
		// spacing is uniform in y
		config_spacing.set("spacingType","uniform");
		config_spacing.set("xmin", 90.0-(shiftLat_ ? 90.0/(nlat-1) : 0.0) );
		config_spacing.set("xmax",-90.0+(shiftLat_ ? 90.0/(nlat-1) : 0.0) );
		config_spacing.set("N",nlat+(shiftLat_?-1:0));
		spacing_y_=spacing::Spacing::create(config_spacing);
		
		// domain is global
		config_domain.set("domainType","global");
		domain_=domain::Domain::create(config_domain);
		
		// setup regular grid
    Regular::setup();
    
}

RegularLonLat::RegularLonLat(const util::Config& config) :
    Regular()
{
		long nlon, nlat, N;
		
		// dimensions
		if ( config.get("N",N) ) {
			nlon=4*N;nlat=2*N;
		} else {
			if ( !config.get("nlon",nlon) || !config.get("nlat",nlat) ) {
				throw eckit::BadParameter("RegularLonLat requires either N, or (nlon,nlat)",Here());
			}
		}
		
		// default: no shift
		shiftLon_=false;
		shiftLat_=false;

		// perform setup
		setup(nlon, nlat);
}

RegularLonLat::RegularLonLat() : Regular() {
}

eckit::Properties RegularLonLat::spec() const {
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

