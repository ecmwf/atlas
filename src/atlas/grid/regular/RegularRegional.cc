#include "atlas/grid/regular/RegularRegional.h"
#include "atlas/grid/domain/RectangularDomain.h"


namespace atlas {
namespace grid {
namespace regular {

register_BuilderT1(Grid,RegularRegional,RegularRegional::grid_type_str());

std::string RegularRegional::grid_type_str() {
    return "regularRegional";
}


std::string RegularRegional::className() {
    return "atlas.grid.regular.RegularRegional";
}


void RegularRegional::setup(const util::Config& config) {

		util::Config config_dom, config_spacing_x, config_spacing_y, config_proj;
		
		long nx, ny;
		std::vector<double> bbox(4);
		std::vector<double> sw(2), ne(2), center(2);
		double dx, dy;
		
		// read subconfigurations
		
		// projection
		if ( config.get("projection",config_proj) ) {
		} else {
			// default, error, or hardcoded default?
			config_proj.set("projectionType","lonlat");
		}
		projection_=projection::Projection::create(config_proj);

		// dimensions
		if ( ! config.get("nx",nx) ) throw eckit::BadParameter("nx missing in Params",Here());
		if ( ! config.get("ny",ny) ) throw eckit::BadParameter("ny missing in Params",Here());

		// domain
		if ( config.get("domain",config_dom) ) {
			// domain is specified either by bbox, by sw and ne, or by center and resolution
			if ( !config_dom.get("bbox",bbox) ) {
				if ( config_dom.get("center",center) && config_dom.get("dx",dx) && config_dom.get("dy",dy)	) {
					// coordinates of center
					eckit::geometry::LLPoint2 llcenter(center[0],center[1]);
					eckit::geometry::Point2 xy;
					xy=projection_->lonlat2coords(llcenter);
					// calculate bbox
					bbox[0]=xy[0]-(nx-1)*dx/2;bbox[1]=xy[0]+(nx-1)*dx/2;
					bbox[2]=xy[1]-(ny-1)*dy/2;bbox[3]=xy[1]+(ny-1)*dy/2;
				} else {
					if ( config_dom.get("sw",sw) && config_dom.get("ne",ne) )	{
						// corner coordinates
						eckit::geometry::Point2 xy;
						eckit::geometry::LLPoint2 llsw(sw[0],sw[1]), llne(ne[0],ne[1]);
						// sw corner
						xy=projection_->lonlat2coords(llsw);
						// put in bbox
						bbox[0]=xy[0];bbox[2]=xy[1];
						// ne corner
						xy=projection_->lonlat2coords(llne);
						// put in bbox
						bbox[1]=xy[0];bbox[3]=xy[1];
					} else {
						throw eckit::BadParameter("RegularRegional grid domain should be specified by (i) bbox, (ii) center, dx and dy, or (iii) ne and sw",Here());		
					}
				}
				// put bbox in config_dom
				config_dom.set("bbox",bbox);
			}
			// set domainType if it's missing
			std::string domainType;
			if (!config_dom.get("domainType",domainType)) config_dom.set("domainType","rectangular");
			// create domain from configuration
			domain_=domain::Domain::create(config_dom);
			// check if the domain is rectangular
			domain::RectangularDomain * rd=dynamic_cast<domain::RectangularDomain*>(domain_);
			if (! rd) throw eckit::BadParameter("RegularRegional grid requires a RectangularDomain",Here());
				
		} else {
			// error
			throw eckit::BadParameter("domain is required for a RegularRegional grid",Here());
		}

		// spacing_x
		util::Config config_spacing;
		std::string spacingType;
		
		if ( !config.get("spacing_x",spacingType) ) spacingType="uniform";
		// set configuration of spacing_x
		config_spacing.set("spacingType",spacingType);
		config_spacing.set("xmin",bbox[0]);
		config_spacing.set("xmax",bbox[1]);
		config_spacing.set("N",nx);
		// create spacing
		spacing_x_=spacing::Spacing::create(config_spacing);
		
		// spacing_y
		if ( !config.get("spacing_y",spacingType) ) spacingType="uniform";
		// set configuration of spacing_y
		config_spacing.set("spacingType",spacingType);
		config_spacing.set("xmin",bbox[2]);
		config_spacing.set("xmax",bbox[3]);
		config_spacing.set("N",ny);
		// create spacing
		spacing_y_=spacing::Spacing::create(config_spacing);
		
	
		// setup regular grid
    Regular::setup();
    
}

RegularRegional::RegularRegional(const util::Config& config) :
    Regular()
{
		setup(config);
}

RegularRegional::RegularRegional() :
    Regular()
{
}

eckit::Properties RegularRegional::spec() const {
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

