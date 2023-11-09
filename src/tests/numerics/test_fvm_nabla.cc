/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <cmath>
#include <iomanip>
#include <iostream>

#include "eckit/config/Resource.h"

#include "atlas/array/MakeView.h"
#include "atlas/field/Field.h"
#include "atlas/field/FieldSet.h"
#include "atlas/grid/Distribution.h"
#include "atlas/grid/Grid.h"
#include "atlas/grid/Partitioner.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/mesh/Nodes.h"
#include "atlas/meshgenerator.h"
#include "atlas/numerics/Nabla.h"
#include "atlas/numerics/fvm/Method.h"
#include "atlas/option.h"
#include "atlas/output/Gmsh.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/util/Config.h"
#include "atlas/util/Constants.h"
#include "atlas/util/CoordinateEnums.h"
#include "atlas/util/Earth.h"

#include "tests/AtlasTestEnvironment.h"

using namespace eckit;
using namespace atlas::numerics;
using namespace atlas::meshgenerator;
using namespace atlas::grid;

namespace atlas {
namespace test {


int test_levels() {
    static int levels = []() {
        int levels = eckit::Resource<int>("--levels", 0);
        Log::info() << "levels = " << levels << std::endl;
        return levels;
    }();
    return levels;
}

static std::string griduid() {
    return "Slat20";
}

template <typename Value>
array::LocalView<Value, 3> make_vectorview(Field& field) {
    using array::Range;
    return field.levels() ? array::make_view<Value, 3>(field).slice(Range::all(), Range::all(), Range::all())
                          : array::make_view<Value, 2>(field).slice(Range::all(), Range::dummy(), Range::all());
}

template <typename Value>
array::LocalView<Value, 2> make_scalarview(Field& field) {
    using array::Range;
    return field.levels() ? array::make_view<Value, 2>(field).slice(Range::all(), Range::all())
                          : array::make_view<Value, 1>(field).slice(Range::all(), Range::dummy());
}


#define print_min_max_mean(name)                                                \
    do {                                                                        \
        Log::info() << #name << std::endl;                                      \
        Log::info() << std::setprecision(18) << "  min  " << min << std::endl;  \
        Log::info() << std::setprecision(18) << "  max  " << max << std::endl;  \
        Log::info() << std::setprecision(18) << "  mean " << mean << std::endl; \
    } while (0)


//-----------------------------------------------------------------------------

double dual_volume(const Mesh& mesh) {
    const mesh::Nodes& nodes = mesh.nodes();
    int nb_nodes             = nodes.size();
    auto dual_volumes        = array::make_view<double, 1>(nodes.field("dual_volumes"));
    auto is_ghost            = array::make_view<int, 1>(nodes.ghost());
    double area              = 0;
    for (int node = 0; node < nb_nodes; ++node) {
        if (!is_ghost(node)) {
            area += dual_volumes(node);
        }
    }

    mpi::comm().allReduceInPlace(area, eckit::mpi::sum());

    return area;
}

/// @brief Compute magnitude of flow with rotation-angle beta
/// (beta=0 --> zonal, beta=pi/2 --> meridional)
template <typename Value>
void rotated_flow(const fvm::Method& fvm, Field& field, const double& beta) {
    const double radius  = fvm.radius();
    const double USCAL   = 20.;
    const double pvel    = USCAL / radius;
    const double deg2rad = M_PI / 180.;

    auto lonlat_deg = array::make_view<double, 2>(fvm.mesh().nodes().lonlat());
    auto var        = make_vectorview<Value>(field);

    const idx_t nnodes = fvm.mesh().nodes().size();
    const idx_t nlev   = var.shape(1);
    for (idx_t jnode = 0; jnode < nnodes; ++jnode) {
        double x  = lonlat_deg(jnode, LON) * deg2rad;
        double y  = lonlat_deg(jnode, LAT) * deg2rad;
        double Ux = pvel * (std::cos(beta) + std::tan(y) * std::cos(x) * std::sin(beta)) * radius * std::cos(y);
        double Uy = -pvel * std::sin(x) * std::sin(beta) * radius;
        for (idx_t jlev = 0; jlev < nlev; ++jlev) {
            var(jnode, jlev, LON) = Ux;
            var(jnode, jlev, LAT) = Uy;
        }
    }
}

/// @brief Compute magnitude of flow with rotation-angle beta
/// (beta=0 --> zonal, beta=pi/2 --> meridional)
template <typename Value>
void rotated_flow_magnitude(const fvm::Method& fvm, Field& field, const double& beta) {
    const double radius  = fvm.radius();
    const double USCAL   = 20.;
    const double pvel    = USCAL / radius;
    const double deg2rad = M_PI / 180.;

    auto lonlat_deg = array::make_view<double, 2>(fvm.mesh().nodes().lonlat());
    auto var        = make_scalarview<Value>(field);

    const idx_t nnodes = fvm.mesh().nodes().size();
    const idx_t nlev   = var.shape(1);
    for (idx_t jnode = 0; jnode < nnodes; ++jnode) {
        double x  = lonlat_deg(jnode, LON) * deg2rad;
        double y  = lonlat_deg(jnode, LAT) * deg2rad;
        double Ux = pvel * (std::cos(beta) + std::tan(y) * std::cos(x) * std::sin(beta)) * radius * std::cos(y);
        double Uy = -pvel * std::sin(x) * std::sin(beta) * radius;
        for (idx_t jlev = 0; jlev < nlev; ++jlev) {
            var(jnode, jlev) = std::sqrt(Ux * Ux + Uy * Uy);
        }
    }
}

//-----------------------------------------------------------------------------

CASE("test_factory") {
    EXPECT(NablaFactory::has("fvm"));
}

CASE("test_build") {
    Log::info() << "test_build" << std::endl;
    MeshGenerator meshgenerator("structured");
    Mesh mesh      = meshgenerator.generate(Grid("O16"));
    const double R = util::Earth::radius();
    fvm::Method fvm(mesh, util::Config("radius", R));
    Nabla nabla(fvm);

    double spherical_area = 360. * 180.;
    EXPECT(eckit::types::is_approximately_equal(dual_volume(mesh), spherical_area, 5.0));
}

CASE("test_grad") {
    auto radius = option::radius("Earth");
    Grid grid(griduid());
    MeshGenerator meshgenerator("structured");
    Mesh mesh = meshgenerator.generate(grid, Distribution(grid, Partitioner("equal_regions")));
    fvm::Method fvm(mesh, radius | option::levels(test_levels()));
    Nabla nabla(fvm);

    idx_t nnodes = mesh.nodes().size();
    idx_t nlev   = std::max(1, test_levels());

    auto do_test = [&](auto value, double tolerance) {
        using Value = std::decay_t<decltype(value)>;
        FieldSet fields;
        fields.add(fvm.node_columns().createField<Value>(option::name("scalar")));
        fields.add(fvm.node_columns().createField<Value>(option::name("rscalar")));
        fields.add(fvm.node_columns().createField<Value>(option::name("grad") | option::variables(2)));
        fields.add(fvm.node_columns().createField<Value>(option::name("rgrad") | option::variables(2)));
        fields.add(fvm.node_columns().createField<Value>(option::name("xder")));
        fields.add(fvm.node_columns().createField<Value>(option::name("yder")));
        fields.add(fvm.node_columns().createField<Value>(option::name("rxder")));
        fields.add(fvm.node_columns().createField<Value>(option::name("ryder")));

        EXPECT(fields["scalar"].rank() == (test_levels() ? 2 : 1));
        EXPECT(fields["grad"].rank() == fields["scalar"].rank() + 1);
        EXPECT(fields["scalar"].levels() == test_levels());
        EXPECT(fields["grad"].levels() == test_levels());

        rotated_flow_magnitude<Value>(fvm, fields["scalar"], 0.);
        rotated_flow_magnitude<Value>(fvm, fields["rscalar"], M_PI_2 * 0.75);

        nabla.gradient(fields["scalar"], fields["grad"]);
        nabla.gradient(fields["rscalar"], fields["rgrad"]);
        auto xder        = make_scalarview<Value>(fields["xder"]);
        auto yder        = make_scalarview<Value>(fields["yder"]);
        auto rxder       = make_scalarview<Value>(fields["rxder"]);
        auto ryder       = make_scalarview<Value>(fields["ryder"]);
        const auto grad  = make_vectorview<Value>(fields["grad"]);
        const auto rgrad = make_vectorview<Value>(fields["rgrad"]);
        for (idx_t jnode = 0; jnode < nnodes; ++jnode) {
            for (idx_t jlev = 0; jlev < nlev; ++jlev) {
                xder(jnode, jlev)  = grad(jnode, jlev, LON);
                yder(jnode, jlev)  = grad(jnode, jlev, LAT);
                rxder(jnode, jlev) = rgrad(jnode, jlev, LON);
                ryder(jnode, jlev) = rgrad(jnode, jlev, LAT);
            }
        }

        // output to gmsh
        {
            fvm.node_columns().haloExchange(fields);
            output::Gmsh(grid.name() + ".msh").write(mesh);
            output::Gmsh gmsh_fields(grid.name() + "_fields.msh");
            gmsh_fields.write(fields["scalar"]);
            gmsh_fields.write(fields["xder"]);
            gmsh_fields.write(fields["yder"]);
            gmsh_fields.write(fields["rscalar"]);
            gmsh_fields.write(fields["rxder"]);
            gmsh_fields.write(fields["ryder"]);
        }

        double min, max, mean;
        idx_t N;

        fvm.node_columns().minimum(fields["xder"], min);
        fvm.node_columns().maximum(fields["xder"], max);
        fvm.node_columns().mean(fields["xder"], mean, N);
        print_min_max_mean("xder");
        EXPECT_APPROX_EQ(min, 0., tolerance);
        EXPECT_APPROX_EQ(max, 0., tolerance);
        EXPECT_APPROX_EQ(mean, 0., tolerance);

        fvm.node_columns().minimum(fields["yder"], min);
        fvm.node_columns().maximum(fields["yder"], max);
        fvm.node_columns().mean(fields["yder"], mean, N);
        print_min_max_mean("yder");
        EXPECT_APPROX_EQ(min, -3.1141489788326316614e-06, tolerance);
        EXPECT_APPROX_EQ(max, 3.1141489788326316614e-06, tolerance);
        EXPECT_APPROX_EQ(mean, 0., tolerance);


        fvm.node_columns().minimum(fields["rxder"], min);
        fvm.node_columns().maximum(fields["rxder"], max);
        fvm.node_columns().mean(fields["rxder"], mean, N);
        print_min_max_mean("rxder");
        EXPECT_APPROX_EQ(min, -3.02863817262107e-06, tolerance);
        EXPECT_APPROX_EQ(max, +3.02863817262107e-06, tolerance);
        EXPECT_APPROX_EQ(mean, 0., tolerance);

        fvm.node_columns().minimum(fields["ryder"], min);
        fvm.node_columns().maximum(fields["ryder"], max);
        fvm.node_columns().mean(fields["ryder"], mean, N);
        print_min_max_mean("ryder");
        EXPECT_APPROX_EQ(min, -3.114148978832633e-06, tolerance);
        EXPECT_APPROX_EQ(max, +3.114148978832633e-06, tolerance);
        EXPECT_APPROX_EQ(mean, 0., tolerance);
    };
    SECTION("double precision") {
        do_test(double{}, 1.e-20);
    }
   SECTION("single precision") {
        do_test(float{}, 1.e-10);
    }
}


CASE("test_div") {
    const double radius = util::Earth::radius();
    //  const double radius = 1.;
    Grid grid(griduid());
    MeshGenerator meshgenerator("structured");
    Mesh mesh = meshgenerator.generate(grid, Distribution(grid, Partitioner("equal_regions")));
    fvm::Method fvm(mesh, util::Config("radius", radius) | option::levels(test_levels()));
    Nabla nabla(fvm);

    auto do_test = [&](auto value, double tolerance) {
        using Value = std::decay_t<decltype(value)>;

        FieldSet fields;
        fields.add(fvm.node_columns().createField<Value>(option::name("wind") | option::variables(2)));
        fields.add(fvm.node_columns().createField<Value>(option::name("div")));

        rotated_flow<Value>(fvm, fields["wind"], M_PI_2 * 0.75);

        nabla.divergence(fields["wind"], fields["div"]);

        // output to gmsh
        {
            fvm.node_columns().haloExchange(fields);
            output::Gmsh gmsh(grid.name() + "_fields.msh", "a");
            gmsh.write(fields["wind"]);
            gmsh.write(fields["div"]);
        }

        double min, max, mean;
        idx_t N;
        fvm.node_columns().minimum(fields["div"], min);
        fvm.node_columns().maximum(fields["div"], max);
        fvm.node_columns().mean(fields["div"], mean, N);

        // Divergence free flow!
        print_min_max_mean("div");
        EXPECT_APPROX_EQ(min, 0.,tolerance);
        EXPECT_APPROX_EQ(max, 0., tolerance);
        EXPECT_APPROX_EQ(mean, 0., tolerance);
    };
    SECTION("double precision") {
        do_test(double{}, 1.e-18);
    }
    SECTION("single precision") {
        do_test(float{}, 1.e-10);
    }
    
}

CASE("test_curl") {
    Log::info() << "test_curl" << std::endl;
    const double radius = util::Earth::radius();
    //  const double radius = 1.;
    Grid grid(griduid());
    MeshGenerator meshgenerator("structured");
    Mesh mesh = meshgenerator.generate(grid, Distribution(grid, Partitioner("equal_regions")));
    fvm::Method fvm(mesh, util::Config("radius", radius) | option::levels(test_levels()));
    Nabla nabla(fvm);

    auto do_test = [&](auto value, double tolerance) {
        using Value = std::decay_t<decltype(value)>;

        FieldSet fields;
        fields.add(fvm.node_columns().createField<Value>(option::name("wind") | option::variables(2)));
        fields.add(fvm.node_columns().createField<Value>(option::name("vor")));

        rotated_flow<Value>(fvm, fields["wind"], M_PI_2 * 0.75);

        nabla.curl(fields["wind"], fields["vor"]);

        fields.add(fvm.node_columns().createField<Value>(option::name("windgrad") | option::variables(2 * 2)));
        nabla.gradient(fields["wind"], fields["windgrad"]);

        fields.add(fvm.node_columns().createField<Value>(option::name("windX") | option::levels(false)));
        fields.add(fvm.node_columns().createField<Value>(option::name("windY") | option::levels(false)));
        fields.add(fvm.node_columns().createField<Value>(option::name("windXgradX")));
        fields.add(fvm.node_columns().createField<Value>(option::name("windXgradY")));
        fields.add(fvm.node_columns().createField<Value>(option::name("windYgradX")));
        fields.add(fvm.node_columns().createField<Value>(option::name("windYgradY")));
        auto wind     = make_vectorview<Value>(fields["wind"]);
        auto windgrad = make_vectorview<Value>(fields["windgrad"]);

        auto windX      = array::make_view<Value, 1>(fields["windX"]);
        auto windY      = array::make_view<Value, 1>(fields["windY"]);
        auto windXgradX = make_scalarview<Value>(fields["windXgradX"]);
        auto windXgradY = make_scalarview<Value>(fields["windXgradY"]);
        auto windYgradX = make_scalarview<Value>(fields["windYgradX"]);
        auto windYgradY = make_scalarview<Value>(fields["windYgradY"]);
        for (idx_t j = 0; j < windX.size(); ++j) {
            static const idx_t lev0 = 0;
            static const idx_t XdX  = XX * 2 + XX;
            static const idx_t XdY  = XX * 2 + YY;
            static const idx_t YdX  = YY * 2 + XX;
            static const idx_t YdY  = YY * 2 + YY;
            windX(j)                = wind(j, lev0, XX);
            windY(j)                = wind(j, lev0, YY);
            windXgradX(j, lev0)     = windgrad(j, lev0, XdX);
            windXgradY(j, lev0)     = windgrad(j, lev0, XdY);
            windYgradX(j, lev0)     = windgrad(j, lev0, YdX);
            windYgradY(j, lev0)     = windgrad(j, lev0, YdY);
        }

        // output to gmsh
        {
            fvm.node_columns().haloExchange(fields);
            output::Gmsh gmsh(grid.name() + "_fields.msh", "a");
            gmsh.write(fields["vor"]);
            gmsh.write(fields["windX"]);
            gmsh.write(fields["windXgradX"]);
            gmsh.write(fields["windXgradY"]);
            gmsh.write(fields["windY"]);
            gmsh.write(fields["windYgradX"]);
            gmsh.write(fields["windYgradY"]);
            gmsh.write(fields["windgrad"]);
        }

        double min, max, mean;
        idx_t N;

        // Vorticity!
        fvm.node_columns().minimum(fields["vor"], min);
        fvm.node_columns().maximum(fields["vor"], max);
        fvm.node_columns().mean(fields["vor"], mean, N);
        print_min_max_mean("vor");
        EXPECT_APPROX_EQ(min, -6.257451225821150e-06, tolerance);
        EXPECT_APPROX_EQ(max, 6.257451225821150e-06, tolerance);
        EXPECT_APPROX_EQ(mean, 0., tolerance);
    };
    SECTION("double precision") {
        do_test(double{}, 1.e-18);
    }
    SECTION("single precision") {
        do_test(float{}, 1.e-10);
    }
}

CASE("test_lapl") {
    Log::info() << "test_lapl" << std::endl;
    const double radius = util::Earth::radius();
    //  const double radius = 1.;
    Grid grid(griduid());
    MeshGenerator meshgenerator("structured");
    Mesh mesh = meshgenerator.generate(grid, Distribution(grid, Partitioner("equal_regions")));
    fvm::Method fvm(mesh, util::Config("radius", radius) | option::levels(test_levels()));
    Nabla nabla(fvm);


    auto do_test = [&](auto value, double tolerance) {
        using Value = std::decay_t<decltype(value)>;

        FieldSet fields;
        fields.add(fvm.node_columns().createField<Value>(option::name("scal")));
        fields.add(fvm.node_columns().createField<Value>(option::name("lapl")));

        rotated_flow_magnitude<Value>(fvm, fields["scal"], M_PI_2 * 0.75);

        nabla.laplacian(fields["scal"], fields["lapl"]);

        // output to gmsh
        {
            fvm.node_columns().haloExchange(fields);
            output::Gmsh gmsh(grid.name() + "_fields.msh", "a");
            gmsh.write(fields["lapl"]);
        }

        double min, max, mean;
        idx_t N;
        fvm.node_columns().minimum(fields["lapl"], min);
        fvm.node_columns().maximum(fields["lapl"], max);
        fvm.node_columns().mean(fields["lapl"], mean, N);
        print_min_max_mean("lapl");
        EXPECT_APPROX_EQ(min, -6.4088005677811607095e-13, tolerance);
        EXPECT_APPROX_EQ(max, 9.8984499569639476135e-12, tolerance);
        EXPECT_APPROX_EQ(mean, -1.03409e-13, tolerance);
    };
    SECTION("double precision") {
        do_test(double{}, 1.e-16);
    }
    SECTION("single precision") {
        do_test(float{}, 1.e-10);
    }

}

//-----------------------------------------------------------------------------

}  // namespace test
}  // namespace atlas

int main(int argc, char** argv) {
    return atlas::test::run(argc, argv);
}
