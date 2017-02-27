#include "atlas/grid/Regular.h"

namespace atlas {
namespace grid {

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

}  // namespace grid
}  // namespace atlas
