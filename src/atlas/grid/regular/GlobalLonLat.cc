#include "atlas/grid/regular/GlobalLonLat.h"
#include "atlas/grid/projection/LonLatProjection.h"

namespace atlas {
namespace grid {
namespace regular {

void GlobalLonLat::setup(long nlon, long nlat) {

    periodic_x_ = true;
    periodic_y_ = false;

    util::Config config_proj, config_spacing, config_domain;

    // projection is lonlat
    config_proj.set("type","lonlat");
    projection_.reset( projection::Projection::create(config_proj) );

    // spacing is uniform in x
    config_spacing.set("spacingType","uniform");
    config_spacing.set("xmin",(shiftLon_ ? 0.5 : 0.0)*360.0/double(nlon) );
    config_spacing.set("xmax",(shiftLon_ ? nlon-0.5 : nlon-1)*360.0/double(nlon) );
    config_spacing.set("N",nlon);
    spacing_x_.reset( spacing::Spacing::create(config_spacing) );

    // spacing is uniform in y
    config_spacing.set("spacingType","uniform");
    config_spacing.set("xmin", 90.0-(shiftLat_ ? 90.0/double(nlat) : 0.0) );
    config_spacing.set("xmax",-90.0+(shiftLat_ ? 90.0/double(nlat) : 0.0) );
    config_spacing.set("N",nlat);
    spacing_y_.reset( spacing::Spacing::create(config_spacing) );

    // domain is global
    config_domain.set("type","global");
    domain_.reset( domain::Domain::create(config_domain) );

    // setup regular grid
    Regular::setup();

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


}  // namespace regular
}  // namespace grid
}  // namespace atlas

