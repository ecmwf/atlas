#include "atlas/grid/reduced/ReducedGaussian.h"

#include "atlas/grid/spacing/Spacing.h"

namespace atlas {
namespace grid {
namespace reduced {

register_BuilderT1(Grid,ReducedGaussian,ReducedGaussian::grid_type_str());

std::string ReducedGaussian::grid_type_str() {
    return "reducedGaussian";
}


std::string ReducedGaussian::className() {
    return "atlas.grid.reduced.ReducedGaussian";
}

void ReducedGaussian::setup(size_t N, long pl[]) {

		util::Config config_spacing, config_domain, config_proj;

		// domain is global
		config_domain.set("domainType","global");
		domain_=domain::Domain::create(config_domain);

		// determine input for Structured::setup
		size_t ny=2*N;
		std::vector<double> xmin(ny);		// first longitude per latitude
		std::vector<double> xmax(ny);		// last longitude per latitude
		std::vector<double> y(ny);			// latitudes
		
    // latitudes: gaussian spacing
    config_spacing.set("spacingType","gaussian");
    config_spacing.set("xmin",90.0);
    config_spacing.set("xmax",-90.0);
    config_spacing.set("N",ny);
    spacing::Spacing * spacing_y=spacing::Spacing::create(config_spacing);
    spacing_y->generate(y);
    
    // loop over latitudes to set bounds
    xmin.assign(ny,0.0);
    for (int jlat=0;jlat<ny;jlat++) {
    	xmin[jlat]=0.0;
    	xmax[jlat]=(pl[jlat]-1)*360.0/pl[jlat];
    }
		
		// setup Structured grid
		Structured::setup(ny,y.data(), pl, xmin.data(), xmax.data());
		Structured::N_=N;

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
    
    // mirror if necessary
    if ( pl.size() == N ) {
    	// mirror to get length 2*N
    	pl.resize(2*N);
    	for (int jlat=0;jlat<N;jlat++) {
    		pl[2*N-1-jlat]=pl[jlat];
    	}
    }
    
    // check length
   	if ( pl.size() != 2*N )
   		throw eckit::BadParameter("pl should have length N or 2*N",Here());

		// projection is lonlat
		util::Config config_proj;
		config_proj.set("projectionType","lonlat");
		projection_=projection::Projection::create(config_proj);
    
    // setup
    setup(N,pl.data());
}

eckit::Properties ReducedGaussian::spec() const {
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

