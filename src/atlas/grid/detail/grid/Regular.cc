#include "atlas/grid/detail/grid/Regular.h"

namespace atlas {
namespace grid {
namespace detail {
namespace grid {

std::string Regular::static_type() {
    return "regular";
}

std::string Regular::className() {
    return "atlas.grid.regular.Regular";
}

Regular::Regular() :
    Structured() {
}

Grid::Spec Regular::spec() const {
    Grid::Spec grid_spec;

    // general specs
    grid_spec=Grid::spec();

    // specific specs
    grid_spec.set("nlat",nlat());
    grid_spec.set("nlon",nlonmin());

    return grid_spec;
}

}  // namespace grid
}  // namespace detail
}  // namespace grid
}  // namespace atlas
