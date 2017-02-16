#include "atlas/grid/regular/Regular.h"

#include "atlas/grid/spacing/UniformSpacing.h"

namespace atlas {
namespace grid {
namespace regular {

register_BuilderT1(Grid,Regular,Regular::grid_type_str());

std::string Regular::grid_type_str() {
    return "regular";
}


std::string Regular::className() {
    return "atlas.grid.regular.Regular";
}

void Regular::setup() {

  // perform checks

  // UniformSpacing in x-direction? -- For now, Structured assumes equidistant points along each latitude
  if( not dynamic_cast<const spacing::UniformSpacing*>(spacing_x_.get()) ) {
    throw eckit::BadParameter("(For now,) Structured grids require a UniformSpacing in X-direction",Here());
  }

  // calculate input for Structured grid
  size_t nx=spacing_x_->size();
  size_t ny=spacing_y_->size();

  std::vector<double> xmin(ny);
  std::vector<double> xmax(ny);
  std::vector<long>   pl(ny);

  xmin.assign( ny, (*spacing_x_)[0]    );
  xmax.assign( ny, (*spacing_x_)[nx-1] );

  pl.assign( ny, nx );

  Structured::setup(ny,spacing_y_->data(), pl.data(), xmin.data(), xmax.data());
  Structured::N_=0;  // not true for some global lonlat grids ;-(

  //set_typeinfo();
}

Regular::Regular() :
    Structured() {
}

Regular::Regular(const util::Config& p) :
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

