#include "atlas/grid/regular/Regular.h"

#include "atlas/grid/spacing/LinearSpacing.h"

namespace atlas {
namespace grid {
namespace regular {


std::string Regular::grid_type_str() {
    return "regular";
}

std::string Regular::className() {
    return "atlas.grid.regular.Regular";
}

Regular::Regular() :
    Structured() {
}

eckit::Properties Regular::spec() const {
    eckit::Properties grid_spec;

    // general specs
    grid_spec=Grid::spec();

    // specific specs
    grid_spec.set("nlat",nlat());
    grid_spec.set("nlon",nlonmin());

    return grid_spec;
}


}  // namespace regular
}  // namespace grid
}  // namespace atlas

