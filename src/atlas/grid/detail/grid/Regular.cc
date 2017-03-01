#include "atlas/grid/detail/grid/Regular.h"

namespace atlas {
namespace grid {
namespace detail {
namespace grid {

std::string Regular::static_type() {
    return "regular";
}

Regular::Regular() :
    Structured() {
}

Grid::Spec Regular::spec() const {
    Grid::Spec grid_spec;

    // general specs
    grid_spec=Grid::spec();

    // specific specs
    grid_spec.set("nx",nlon());
    grid_spec.set("ny",nlat());

    return grid_spec;
}

}  // namespace grid
}  // namespace detail
}  // namespace grid
}  // namespace atlas
