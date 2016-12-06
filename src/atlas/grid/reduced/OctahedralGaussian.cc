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

		// projection is lonlat
		config_proj.set("projectionType","lonlat");
		projection_=projection::Projection::create(config_proj);
		
		// domain is global
		config_domain.set("domainType","global");
		domain_=domain::Domain::create(config_domain);

		// determine input for Structured::setup
		size_t ny=2*N;
		std::vector<double> xmin(ny);		// first longitude per latitude
		std::vector<double> xmax(ny);		// last longitude per latitude
		std::vector<long>   pl(ny);			// number of longitudes per latitude
		std::vector<double> y(ny);			// latitudes
		
		// number of longitudes per latitude
    pl=computePL(N);
    
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
		Structured::setup(ny,y.data(), pl.data(), xmin.data(), xmax.data());
		Structured::N_=N;

}

OctahedralGaussian::OctahedralGaussian(const util::Config& config) :
    Structured()
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

