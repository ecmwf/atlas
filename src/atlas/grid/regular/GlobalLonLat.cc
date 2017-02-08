#include "atlas/grid/regular/GlobalLonLat.h"


namespace atlas {
namespace grid {
namespace regular {

register_BuilderT1(Grid,GlobalLonLat,GlobalLonLat::grid_type_str());

std::string GlobalLonLat::grid_type_str() {
    return "globalLonLat";
}


std::string GlobalLonLat::className() {
    return "atlas.grid.regular.GlobalLonLat";
}

std::string GlobalLonLat::shortName() const {
    std::ostringstream s;
    if ( nlonmin() == 2*nlat() && nlat()%2==0 ) {
      s << "L"<< nlat()/2;
    } else {
      s << "L"<< nlonmin() << "x" << nlat();
    }
    return s.str();
}

void GlobalLonLat::setup(long nlon, long nlat) {

    util::Config config_proj, config_spacing, config_domain;

    // projection is lonlat
    config_proj.set("projectionType","lonlat");
    projection_.reset( projection::Projection::create(config_proj) );

    // spacing is uniform in x
    config_spacing.set("spacingType","uniform");
    config_spacing.set("xmin",(shiftLon_ ? 0.5 : 0.0)*360.0/nlon );
    config_spacing.set("xmax",(shiftLon_ ? nlon-0.5 : nlon-1)*360.0/nlon );
    config_spacing.set("N",nlon);
    spacing_x_.reset( spacing::Spacing::create(config_spacing) );

    // spacing is uniform in y
    config_spacing.set("spacingType","uniform");
    config_spacing.set("xmin", 90.0-(shiftLat_ ? 90.0/(nlat) : 0.0) );
    config_spacing.set("xmax",-90.0+(shiftLat_ ? 90.0/(nlat) : 0.0) );
    config_spacing.set("N",nlat);
    spacing_y_.reset( spacing::Spacing::create(config_spacing) );

    // domain is global
    config_domain.set("domainType","global");
    domain_.reset( domain::Domain::create(config_domain) );

    // setup regular grid
    Regular::setup();

}

GlobalLonLat::GlobalLonLat(const util::Config& config) :
    Regular()
{
    long nlon, nlat, N;

    // dimensions
    if ( config.get("N",N) ) {
      nlon=4*N;nlat=2*N;
    } else {
      if ( !config.get("nlon",nlon) || !config.get("nlat",nlat) ) {
        throw eckit::BadParameter("GlobalLonLat requires either N, or (nlon,nlat)",Here());
      }
    }

    // default: no shift
    shiftLon_=false;
    shiftLat_=false;

    // perform setup
    setup(nlon, nlat);
}

GlobalLonLat::GlobalLonLat(long nlon, long nlat) : Regular() {

    // default: no shift
    shiftLon_=false;
    shiftLat_=false;

    // perform setup
    setup(nlon, nlat);

}

GlobalLonLat::GlobalLonLat(long N) : Regular() {

    long nlon, nlat;

    // dimensions
    nlon=4*N;nlat=2*N;

    // default: no shift
    shiftLon_=false;
    shiftLat_=false;

    // perform setup
    setup(nlon, nlat);

}

GlobalLonLat::GlobalLonLat() : Regular() {
}

eckit::Properties GlobalLonLat::spec() const {
    eckit::Properties grid_spec;

    // general specs
    grid_spec=Grid::spec();

    // no other specs for GlobalLonLat

    return grid_spec;
}

extern "C" {

    Structured* atlas__grid__regular__GlobalLonLat(size_t nlon, size_t nlat) {
        return new GlobalLonLat(nlon,nlat);
    }

}


}  // namespace regular
}  // namespace grid
}  // namespace atlas

