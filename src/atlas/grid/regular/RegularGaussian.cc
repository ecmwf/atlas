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

std::string RegularGaussian::shortName() const {
    std::ostringstream s;
    if ( nlonmin() == 2*nlat() && nlat()%2==0 ) {
      s << "F"<< nlat()/2;
    } else {
      s << "F"<< nlonmin() << "x" << nlat();
    }
    return s.str();
}

void RegularGaussian::setup(long N) {

    util::Config config_proj, config_spacing, config_domain;

    // grid dimensions
    long nlat, nlon;
    nlat=2*N;
    nlon=4*N;

    // projection is lonlat
    config_proj.set("projectionType","lonlat");
    projection_=projection::Projection::create(config_proj);

    // spacing is uniform in x, gaussian in y
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

    // setup regular grid
    Regular::setup();

    // set N_ of Structured
    Structured::N_=N;

}

RegularGaussian::RegularGaussian(long N) : Regular() {

  // perform setup
  setup(N);
}

RegularGaussian::RegularGaussian(const util::Config& config) :
    Regular()
{
    long N;

    // dimensions
    if ( ! config.get("N",N) ) {
      throw eckit::BadParameter("RegularGaussian requires parameter N",Here());
    }

    // perform setup
    setup(N);
}


extern "C" {


    Structured* atlas__grid__regular__RegularGaussian(size_t N) {
        return new RegularGaussian(N);
    }


}

}  // namespace regular
}  // namespace grid
}  // namespace atlas

