#include "atlas/grid/detail/grid/regular/GlobalLonLat.h"
#include "atlas/grid/detail/projection/LonLatProjection.h"
#include "atlas/grid/detail/spacing/LinearSpacing.h"

namespace atlas {
namespace grid {
namespace detail {
namespace grid {
namespace regular {

void GlobalLonLat::setup(long nlon, long nlat, Shift shift) {

    periodic_x_ = true;
    shifted_x_ = shift(Shift::LON);
    shifted_y_ = shift(Shift::LAT);

    util::Config config_proj, config_spacing, config_domain;

    // projection is lonlat
    config_proj.set("type","lonlat");
    projection_ = Projection(config_proj);

    // spacing is uniform in x
    config_spacing.set("type","linear");
    config_spacing.set("start",(shifted_x_ ? 0.5 : 0.0)*360.0/double(nlon) );
    config_spacing.set("end",(shifted_x_ ? nlon-0.5 : nlon-1)*360.0/double(nlon) );
    config_spacing.set("N",nlon);
    spacing::LinearSpacing xspace(config_spacing);

    // spacing is uniform in y
    config_spacing.set("type","linear");
    config_spacing.set("start", 90.0-(shifted_y_ ? 90.0/double(nlat) : 0.0) );
    config_spacing.set("end",-90.0+(shifted_y_ ? 90.0/double(nlat) : 0.0) );
    config_spacing.set("N",nlat);
    YSpace yspace(config_spacing);

    // domain is global
    config_domain.set("type","global");
    domain_ = Domain(config_domain);

    // setup regular grid
    Structured::setup(yspace,xspace.size(),xspace.min(),xspace.max(),xspace.step());

}


GlobalLonLat::GlobalLonLat() : Regular() {
}

Grid::Spec GlobalLonLat::spec() const {
    Grid::Spec grid_spec;

    // general specs
    grid_spec=Grid::spec();

    // no other specs for GlobalLonLat

    return grid_spec;
}


}  // namespace regular
}  // namespace grid
}  // namespace detail
}  // namespace grid
}  // namespace atlas

