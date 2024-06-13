/*
 * (C) Crown Copyright 2024 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */


#include "atlas/array/MakeView.h"
#include "atlas/grid/CubedSphereGrid.h"
#include "atlas/mesh/actions/GetCubedSphereNodalArea.h"
#include "atlas/util/NormaliseLongitude.h"
#include "atlas/util/Point.h"


namespace atlas {
namespace mesh {
namespace actions {


Field& GetCubedSphereNodalArea::operator()(Mesh& mesh) {

  constexpr auto deg2rad = M_PI / 180.;
  const auto& proj = mesh.projection();
  auto lonlat = array::make_view<double, 2>(mesh.nodes().lonlat());

  auto gcell_area_field = Field("grid_cell_areas",
                                make_datatype<double>(),
                                array::make_shape(mesh.nodes().size()));
  auto gcell_area_fview = array::make_view<double, 1>(gcell_area_field);

  double gcell_area_cs = [&] {
    ATLAS_ASSERT(CubedSphereGrid(mesh.grid()));
    // (grid_res * grid_res) = no. of cells on a tile
    auto grid_res = CubedSphereGrid(mesh.grid()).N();
    // area of a grid cell (cubed-sphere coord. system)
    return M_PI/(2*grid_res) * M_PI/(2*grid_res);
  }();

  for (size_t i = 0; i < gcell_area_fview.size(); ++i) {
    PointLonLat loc = PointLonLat(lonlat(i, atlas::LON), lonlat(i, atlas::LAT));
    double cos_lat = std::cos(deg2rad * loc.lat());
    double grid_jac_det = 1/proj.jacobian(loc).determinant();
    // area of a grid cell (geographic coord. system)
    gcell_area_fview(i) = grid_jac_det * gcell_area_cs * cos_lat;
  }

  mesh.nodes().add(gcell_area_field);

  return mesh.nodes().field("grid_cell_areas");

}

}  // namespace actions
}  // namespace mesh
}  // namespace atlas
