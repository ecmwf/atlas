/*
 * (C) Crown Copyright 2024 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */


#include <string>

#include "atlas/array.h"
#include "atlas/field.h"
#include "atlas/field/FieldSet.h"
#include "atlas/functionspace/CubedSphereColumns.h"
#include "atlas/functionspace/NodeColumns.h"
#include "atlas/grid.h"
#include "atlas/interpolation.h"
#include "atlas/mesh.h"
#include "atlas/meshgenerator/MeshGenerator.h"
#include "atlas/option/Options.h"
#include "atlas/output/Gmsh.h"
#include "atlas/util/Config.h"
#include "atlas/util/CoordinateEnums.h"
#include "atlas/util/function/VortexRollup.h"

#include "tests/AtlasTestEnvironment.h"


namespace atlas {
namespace test {


/// function to generate a synthetic field (vortex field)
Field getVortexField(
  const functionspace::NodeColumns& csfs,
  const std::string& vfield_name,
  const size_t nb_levels) {
  const auto lonlat = array::make_view<double, 2>(csfs.lonlat());

  auto vfield = csfs.createField<double>(option::name(vfield_name) |
                                         option::levels(nb_levels));
  auto vfield_view = array::make_view<double, 2>(vfield);
  double mean_field = 1.;
  for (idx_t l=0; l < vfield.shape(1); ++l) {
    for (idx_t i=0; i < vfield.shape(0); ++i) {
      vfield_view(i, l) = mean_field + util::function::vortex_rollup(
        lonlat(i, atlas::LON), lonlat(i, atlas::LAT), 1.);
    }
  }

  return vfield;
}


/// function to write a field set in a Gmsh file
void makeGmshOutput(const std::string& file_name,
                    const FieldSet& fields) {
  const auto& fs = fields[0].functionspace();
  const auto& mesh = functionspace::NodeColumns(fs).mesh();

  const auto config_gmsh =
    util::Config("coordinates", "xyz") |
    util::Config("ghost", true) |
    util::Config("info", true);
  const auto gmsh = output::Gmsh(file_name, config_gmsh);
  gmsh.write(mesh);
  gmsh.write(fields, fs);
}


/// test to carry out the rigridding from 'high' to 'low' resolution,
/// for a given type of grid (CS-LFR)
///
/// after the regridding, the field set that has been generated is stored
/// in a Gmsh file (data visualization);
///
CASE("rigridding from high to low resolution; grid type: CS-LFR") {

  // source grid (high res.)
  auto csgrid_s = Grid("CS-LFR-100");
  auto csmesh_s = MeshGenerator("cubedsphere_dual").generate(csgrid_s);
  auto csfs_s = functionspace::NodeColumns(csmesh_s);

  // target grid (low res.)
  auto csgrid_t = Grid("CS-LFR-50");
  auto csmesh_t = MeshGenerator("cubedsphere_dual").generate(csgrid_t);
  auto csfs_t = functionspace::NodeColumns(csmesh_t);

  size_t nb_levels = 1;

  auto field_01_s = getVortexField(csfs_s, "field_01_s", nb_levels);

  atlas::FieldSet fs_s;
  fs_s.add(field_01_s);

  auto field_01_t = csfs_t.createField<double>(option::name("field_01_t") |
                                               option::levels(nb_levels));

  FieldSet fs_t;
  fs_t.add(field_01_t);


  const auto scheme = util::Config("type", "binning") |
    util::Config("ancillary_scheme", option::type("cubedsphere-bilinear"));

  Interpolation regrid_high2low(scheme, csfs_s, csfs_t);

  // performing the regridding from high to low resolution
  regrid_high2low.execute(fs_s, fs_t);


  //--

  const std::string fname_s("regridding_h2l_field_s.msh");

  makeGmshOutput(fname_s, fs_s["field_01_s"]);

  const std::string fname_t("regridding_h2l_field_t.msh");

  makeGmshOutput(fname_t, fs_t["field_01_t"]);
}


}  // namespace test
}  // namespace atlas


//--

int main(int argc, char** argv) { return atlas::test::run(argc, argv); }
