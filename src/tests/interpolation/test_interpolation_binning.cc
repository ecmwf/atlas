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
#include "atlas/functionspace/StructuredColumns.h"
#include "atlas/grid.h"
#include "atlas/interpolation.h"
#include "atlas/mesh.h"
#include "atlas/meshgenerator/MeshGenerator.h"
#include "atlas/option/Options.h"
#include "atlas/output/Gmsh.h"
#include "atlas/runtime/Log.h"
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
                    const Mesh& mesh,
                    const FieldSet& fields) {
  const auto& fs = fields[0].functionspace();

  const auto config_gmsh =
    util::Config("coordinates", "xyz") |
    util::Config("ghost", true) |
    util::Config("info", true);
  const auto gmsh = output::Gmsh(file_name, config_gmsh);
  gmsh.write(mesh);
  gmsh.write(fields, fs);
}


/// function to carry out a dot product
double dotProd(const Field& field01, const Field& field02) {
  double dprod{};

  const auto field01_view = array::make_view<double, 2>(field01);
  const auto field02_view = array::make_view<double, 2>(field02);

  for (idx_t l=0; l < field01_view.shape(1); ++l) {
    for (idx_t i=0; i < field01_view.shape(0); ++i) {
      dprod += field01_view(i, l) * field02_view(i, l);
    }
  }
  eckit::mpi::comm().allReduceInPlace(dprod, eckit::mpi::Operation::SUM);

  return dprod;
}


//--


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

  FieldSet fs_s;
  fs_s.add(field_01_s);

  auto field_01_t = csfs_t.createField<double>(option::name("field_01_t") |
                                               option::levels(nb_levels));

  FieldSet fs_t;
  fs_t.add(field_01_t);


  const auto scheme = util::Config("type", "binning") |
                      util::Config("scheme", option::type("cubedsphere-bilinear"));

  Interpolation regrid_high2low(scheme, csfs_s, csfs_t);

  // performing the regridding from high to low resolution
  regrid_high2low.execute(fs_s, fs_t);

  fs_t["field_01_t"].haloExchange();


  //--

  const std::string fname_s("regridding_h2l_cs_s.msh");

  makeGmshOutput(fname_s, csmesh_s, fs_s["field_01_s"]);

  const std::string fname_t("regridding_h2l_cs_t.msh");

  makeGmshOutput(fname_t, csmesh_t, fs_t["field_01_t"]);
}


/// test to carry out the rigridding from 'high' to 'low' resolution,
/// for a given type of grid (O)
///
/// after the regridding, the field set that has been generated is stored
/// in a Gmsh file (data visualization);
///
CASE("rigridding from high to low resolution; grid type: O") {

  // source grid (high res.)
  auto ncgrid_s = Grid("O32");
  auto ncmesh_s = MeshGenerator("structured").generate(ncgrid_s);
  auto ncfs_s = functionspace::StructuredColumns(ncgrid_s, option::halo(3));

  // target grid (low res.)
  auto ncgrid_t = Grid("O16");
  auto ncmesh_t = MeshGenerator("structured").generate(ncgrid_t);
  auto ncfs_t = functionspace::StructuredColumns(ncgrid_t, option::halo(3));

  size_t nb_levels = 1;

  auto field_01_s = getVortexField(ncfs_s, "field_01_s", nb_levels);

  FieldSet fs_s;
  fs_s.add(field_01_s);

  auto field_01_t = ncfs_t.createField<double>(option::name("field_01_t") |
                                               option::levels(nb_levels));

  FieldSet fs_t;
  fs_t.add(field_01_t);


  const auto scheme = util::Config("type", "binning") |
                      util::Config("scheme", option::type("structured-linear2D"));

  Interpolation regrid_high2low(scheme, ncfs_s, ncfs_t);

  // performing the regridding from high to low resolution
  regrid_high2low.execute(fs_s, fs_t);

  fs_t["field_01_t"].haloExchange();


  //--

  const std::string fname_s("regridding_h2l_nc_s.msh");

  makeGmshOutput(fname_s, ncmesh_s, fs_s["field_01_s"]);

  const std::string fname_t("regridding_h2l_nc_t.msh");

  makeGmshOutput(fname_t, ncmesh_t, fs_t["field_01_t"]);
}


/// test to carry out the 'dot-product' test for the rigridding from
/// 'high' to 'low' resolution, for a given type of grid (CS-LFR)
///
CASE("dot-product test for the rigridding from high to low resolution; grid type: CS-LFR") {

  // source grid (high res.)
  auto csgrid_s = Grid("CS-LFR-100");
  auto csmesh_s = MeshGenerator("cubedsphere_dual").generate(csgrid_s);
  auto csfs_s = functionspace::NodeColumns(csmesh_s);

  // target grid (low res.)
  auto csgrid_t = Grid("CS-LFR-50");
  auto csmesh_t = MeshGenerator("cubedsphere_dual").generate(csgrid_t);
  auto csfs_t = functionspace::NodeColumns(csmesh_t);

  size_t nb_levels = 1;

  // source field
  auto field_01_s = getVortexField(csfs_s, "field_01_s", nb_levels);

  FieldSet fs_s;
  fs_s.add(field_01_s);

  // target field
  auto field_01_t = csfs_t.createField<double>(option::name("field_01_t") |
                                               option::levels(nb_levels));

  FieldSet fs_t;
  fs_t.add(field_01_t);

  const auto scheme = util::Config("type", "binning") |
                      util::Config("scheme", option::type("cubedsphere-bilinear")) |
                      util::Config("adjoint", true);

  Interpolation regrid_high2low(scheme, csfs_s, csfs_t);

  // performing the regridding from high to low resolution
  regrid_high2low.execute(fs_s, fs_t);


  fs_t["field_01_t"].haloExchange();

  // target field (adjoint)
  auto field_01_ad_t = csfs_t.createField<double>(option::name("field_01_ad_t") |
                                                  option::levels(nb_levels));
  array::make_view<double, 2>(field_01_ad_t).assign(
    array::make_view<double, 2>(field_01_t));
  field_01_ad_t.adjointHaloExchange();

  FieldSet fs_ad_t;
  fs_ad_t.add(field_01_ad_t);

  // source field (adjoint)
  auto field_01_ad_s = csfs_s.createField<double>(option::name("field_01_ad_s") |
                                                  option::levels(nb_levels));
  array::make_view<double, 2>(field_01_ad_s).assign(0.);

  FieldSet fs_ad_s;
  fs_ad_s.add(field_01_ad_s);

  // performing adjoint operation
  regrid_high2low.execute_adjoint(fs_ad_s, fs_ad_t);


  const auto t_dot_t = dotProd(fs_t["field_01_t"], fs_t["field_01_t"]);
  const auto s_dot_ad_s = dotProd(fs_s["field_01_s"], fs_ad_s["field_01_ad_s"]);

  double scaled_diff = std::abs(t_dot_t - s_dot_ad_s)/std::abs(t_dot_t);

  // carrrying out a dot-product test ...
  Log::info() << "\n- dot-product test:\n"
              << "(Ax) . (Ax) = " << t_dot_t << "; "
              << "x . (A^t A x) = " << s_dot_ad_s << "; "
              << "scaled difference = " << scaled_diff << "\n" << std::endl;

  EXPECT(scaled_diff < 1e-12);
}


}  // namespace test
}  // namespace atlas


//--

int main(int argc, char** argv) { return atlas::test::run(argc, argv); }
