#include "atlas/grid/regular/GlobalLonLat.h"
#include "atlas/grid/projection/LonLatProjection.h"
#include "atlas/grid/spacing/LinearSpacing.h"

namespace atlas {
namespace grid {
namespace regular {

void GlobalLonLat::setup(long nlon, long nlat) {

    periodic_x_ = true;

    util::Config config_proj, config_spacing, config_domain;

    // projection is lonlat
    config_proj.set("type","lonlat");
    projection_.reset( projection::Projection::create(config_proj) );

    // spacing is uniform in x
    config_spacing.set("type","linear");
    config_spacing.set("start",(shiftLon_ ? 0.5 : 0.0)*360.0/double(nlon) );
    config_spacing.set("end",(shiftLon_ ? nlon-0.5 : nlon-1)*360.0/double(nlon) );
    config_spacing.set("N",nlon);
    spacing::LinearSpacing xspace(config_spacing);

    // spacing is uniform in y
    config_spacing.set("type","linear");
    config_spacing.set("start", 90.0-(shiftLat_ ? 90.0/double(nlat) : 0.0) );
    config_spacing.set("end",-90.0+(shiftLat_ ? 90.0/double(nlat) : 0.0) );
    config_spacing.set("N",nlat);
    spacing::Spacing* yspace = spacing::Spacing::create(config_spacing);

    // domain is global
    config_domain.set("type","global");
    domain_.reset( domain::Domain::create(config_domain) );

    // setup regular grid
    Structured::setup(yspace,xspace.size(),xspace.min(),xspace.max(),xspace.step());

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

