/*
 * (C) Crown Copyright 2026 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */


#include <cmath>
#include <map>
#include <tuple>
#include <utility>
#include <vector>

#include "atlas/array.h"
#include "atlas/array/helpers/ArrayForEach.h"
#include "atlas/field/Field.h"
#include "atlas/field/FieldSet.h"
#include "atlas/functionspace/NodeColumns.h"
#include "atlas/functionspace/PointCloud.h"
#include "atlas/grid/Grid.h"
#include "atlas/interpolation.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/meshgenerator.h"
#include "atlas/option.h"
#include "atlas/util/Config.h"
#include "atlas/util/Constants.h"
#include "atlas/util/function/VortexRollup.h"
#include "atlas/parallel/mpi/mpi.h"

#include "tests/AtlasTestEnvironment.h"



namespace atlas {
namespace test {


struct FieldSpecsFixture {
  static const util::Config& get(const std::string& fixture) {
        static const auto fieldSpecs = std::map<std::string_view, util::Config>{
            {"source", option::name("sfield") |
	               option::variables(2) |
                       option::type("vector")},
            {"target", option::name("tfield") |
	               option::variables(2) |
	               option::type("vector")}};
        return fieldSpecs.at(fixture);
    }
};

// function to define the 'source' function space
auto generateFunctionspaceSource() {
    const auto csgrid = Grid("CS-LFR-48");
    const auto csmesh = MeshGenerator("cubedsphere_dual").generate(csgrid);
    
    return functionspace::NodeColumns(csmesh);
}

// function to define the 'target' function space
auto generateFunctionspaceTarget() {
    const mpi::Comm& comm_w = mpi::comm();
    std::size_t mpi_rank = comm_w.rank();

    // data structure storing a single point
    // (in a cubed sphere with 6 partitions, the point near
    // the North pole is part of partition #4)
    std::vector<PointXY> spoint;
    if (mpi_rank == 4) {
      spoint.push_back({360., 90.});
    }
       
    return functionspace::PointCloud(spoint);
}

std::pair<double, double> vortexHorizontal(double lon, double lat) {
    // set hLon and hLat step size.
    const double hLon = 0.0001;
    const double hLat = 0.0001;

    const double u = (
      util::function::vortex_rollup(lon, lat + 0.5 * hLat, 0.1) -
      util::function::vortex_rollup(lon, lat - 0.5 * hLat, 0.1)) /
      hLat;

    const double v = -(
      util::function::vortex_rollup(lon + 0.5 * hLon, lat, 0.1) -
      util::function::vortex_rollup(lon - 0.5 * hLon, lat, 0.1)) /
      (hLon * std::cos(lat * util::Constants::degreesToRadians()));

    return std::make_pair(u, v);
}

// function to generate a data structure containing the 'source' field
FieldSet generateFieldsSource(FunctionSpace& fsSource) {
    const auto fieldSourceSpecs = FieldSpecsFixture::get("source");
    FieldSet fieldsSource;
    auto fieldSource =
        fieldsSource.add(fsSource.createField<double>(fieldSourceSpecs));
    auto fieldSourceView = array::make_view<double, 2>(fieldSource);

    const auto lonLatSource =
        array::make_view<double, 2>(fsSource.lonlat());

    array::helpers::ArrayForEach<0>::apply(
        std::tie(lonLatSource, fieldSourceView),
        [](auto&& lonLat, auto&& columnSource) {
            const auto setElems = [&](auto&& elemSource) {
                std::tie(elemSource(0), elemSource(1)) =
                    vortexHorizontal(lonLat(0), lonLat(1));
            };
            setElems(columnSource);
        });
    
    return fieldsSource;
}

// function to generate a data structure containing the 'target' field
FieldSet generateFieldsTarget(FunctionSpace& fs) {
    const auto fieldTargetSpecs = FieldSpecsFixture::get("target");
    FieldSet fieldsTarget;
    fieldsTarget.add(fs.createField<double>(fieldTargetSpecs));
    
    return fieldsTarget;
}

double dotProduct(const array::ArrayView<double, 2>& a,
                  const array::ArrayView<double, 2>& b) {
    auto dotProd = 0.;
    array::helpers::arrayForEachDim(
        std::make_integer_sequence<int, 2>{}, std::tie(a, b),
        [&](const double& aElem, const double& bElem) {
          dotProd += aElem * bElem;
        });
    
    return dotProd;
}

CASE("source grid: cubed sphere (CS-LFR-48); "
     "target grid: point cloud (single point); "
     "field distribution: 3D-field, single level, 2D-vector; "
     "procedure: vector interpolation; "
     "parallelization: MPI (6 PEs)") {
  
    auto fsSource = generateFunctionspaceSource();
    auto fsTarget = generateFunctionspaceTarget();

    auto fieldsSource = generateFieldsSource(fsSource);
    auto fieldSourceView = array::make_view<double, 2>(fieldsSource["sfield"]);
    auto fieldsTarget = generateFieldsTarget(fsTarget);
    auto fieldTargetView = array::make_view<double, 2>(fieldsTarget["tfield"]);
    
    const auto ansScheme = option::type("cubedsphere-bilinear") |
                           util::Config("adjoint", true);
    const auto scheme = util::Config("type", "spherical-vector") |
                        util::Config("adjoint", true) |
                        util::Config("scheme", ansScheme);

    Interpolation interp(scheme, fsSource, fsTarget);

    interp.execute(fieldsSource, fieldsTarget);
    fieldsTarget.haloExchange();

    //--
    
    // performing a dot-product test ...

    const auto fieldTargetSpecs = FieldSpecsFixture::get("target");
    auto adjointTarget = fsTarget.createField<double>(fieldTargetSpecs);
    adjointTarget.array().copy(fieldsTarget["tfield"]);
    adjointTarget.adjointHaloExchange();

    const auto fieldSourceSpecs = FieldSpecsFixture::get("source");
    auto adjointSource = fsSource.createField<double>(fieldSourceSpecs);
    auto adjointSourceView = array::make_view<double, 2>(adjointSource);
    adjointSourceView.assign(0.);

    interp.execute_adjoint(adjointSource, adjointTarget);

    constexpr auto tinyNum = 1e-13;
    const auto targetDotTarget = dotProduct(fieldTargetView, fieldTargetView);
    const auto sourceDotAdjointSource = dotProduct(fieldSourceView, adjointSourceView);

    const mpi::Comm& comm_w = mpi::comm();
    std::size_t mpi_rank = comm_w.rank();

    if (mpi_rank == 4) {
        const auto dotProdRatio = targetDotTarget / sourceDotAdjointSource;
        EXPECT_APPROX_EQ(dotProdRatio, 1., tinyNum);
    }
}


}  // namespace test
}  // namespace atlas


int main(int argc, char** argv) {
    return atlas::test::run(argc, argv);
}

