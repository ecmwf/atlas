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
#include "atlas/util/function/SolidBodyRotation.h"

#include "tests/AtlasTestEnvironment.h"

using namespace eckit;
using namespace atlas::numerics;
using namespace atlas::meshgenerator;
using namespace atlas::grid;

namespace atlas {
namespace test {


// This test relates to JIRA issues METV-2657 , MIR-459

static std::string griduid() {
    //return "O80";
    //return "Slat80";
    return "Slat720x360";
}
static double radius() {
    return util::Earth::radius();
}
static double beta_in_degrees() {
    return 90.;
}

static bool output_gmsh() {
    return false;
}

static bool mask_polar_values() {
    return false;
}

static int metric_approach() {
    // Experimental!!!
    // approach = 0  ORIGINAL, DEFAULT
    //   metric_term cos(y) is multiplied with each wind component
    //   --> cell-interface:  avg = 0.5*( cos(y1)*u1 + cos(y2)*u2 )
    // approach = 1:
    //   metric_term cos(0.5*(y1+y2)) is used at cell interface
    //   --> cell-interface: avg = 0.5*(u1+u2)*cos(0.5*(y1+y2))
    // Results seem to indicate that approach=0 is overall better, although approach=1
    // seems to handle pole slightly better (error factor 2 to 4 times lower)
    return 0;
}

FieldSet analytical_fields(const fvm::Method& fvm) {
    const double radius      = fvm.radius();

    auto lonlat = array::make_view<double, 2>(fvm.mesh().nodes().lonlat());

    FieldSet fields;
    auto add_scalar_field = [&](const std::string& name) {
        return array::make_view<double, 1>(fields.add(fvm.node_columns().createField<double>(option::name(name))));
    };
    auto add_vector_field = [&](const std::string& name) {
        return array::make_view<double, 2>(fields.add(fvm.node_columns().createField<double>(
            option::name(name) | option::type("vector") | option::variables(2))));
    };
    auto f      = add_scalar_field("f");
    auto uv     = add_vector_field("uv");
    auto u      = add_scalar_field("u");
    auto v      = add_scalar_field("v");
    auto grad_f = add_vector_field("ref_grad_f");
    auto dfdx   = add_scalar_field("ref_dfdx");
    auto dfdy   = add_scalar_field("ref_dfdy");
    auto div    = add_scalar_field("ref_div");
    auto vor    = add_scalar_field("ref_vor");

    auto flow          = atlas::util::function::SolidBodyRotation{beta_in_degrees(), radius};
    auto is_ghost      = array::make_view<int, 1>(fvm.mesh().nodes().ghost());
    const idx_t nnodes = fvm.mesh().nodes().size();
    for (idx_t jnode = 0; jnode < nnodes; ++jnode) {
        if (is_ghost(jnode)) {
            continue;
        }
        double x = lonlat(jnode, LON);
        double y = lonlat(jnode, LAT);

        flow.wind(x, y, u(jnode), v(jnode));
        flow.vordiv(x, y, vor(jnode), div(jnode));
        f(jnode) = flow.windMagnitudeSquared(x, y);
        flow.windMagnitudeSquaredGradient(x, y, dfdx(jnode), dfdy(jnode));

        uv(jnode, XX)     = u(jnode);
        uv(jnode, YY)     = v(jnode);
        grad_f(jnode, XX) = dfdx(jnode);
        grad_f(jnode, YY) = dfdy(jnode);
    }
    fields.set_dirty();
    fields.haloExchange();

    return fields;
}

//-----------------------------------------------------------------------------

CASE("test_analytical") {
    Grid grid(griduid(), GlobalDomain(-180.));

    Mesh mesh = MeshGenerator{"structured"}.generate(grid);
    fvm::Method fvm(mesh, option::radius(radius()));
    Nabla nabla(fvm, util::Config("metric_approach", metric_approach()));
    FieldSet fields = analytical_fields(fvm);

    Field div    = fields.add(fvm.node_columns().createField<double>(option::name("div")));
    Field vor    = fields.add(fvm.node_columns().createField<double>(option::name("vor")));
    Field grad_f = fields.add(fvm.node_columns().createField<double>(option::name("grad_f") | option::variables(2)));
    Field dfdx   = fields.add(fvm.node_columns().createField<double>(option::name("dfdx")));
    Field dfdy   = fields.add(fvm.node_columns().createField<double>(option::name("dfdy")));

    auto split = [](const Field& vector, Field& component_x, Field& component_y) {
        auto v = array::make_view<double, 2>(vector);
        auto x = array::make_view<double, 1>(component_x);
        auto y = array::make_view<double, 1>(component_y);
        for (idx_t j = 0; j < v.shape(0); ++j) {
            x(j) = v(j, XX);
            y(j) = v(j, YY);
        }
    };


    ATLAS_TRACE_SCOPE("gradient") nabla.gradient(fields["f"], grad_f);
    split(grad_f, dfdx, dfdy);
    ATLAS_TRACE_SCOPE("divergence") nabla.divergence(fields["uv"], div);
    ATLAS_TRACE_SCOPE("vorticity") nabla.curl(fields["uv"], vor);

    auto do_mask_polar_values = [&](Field& field, double mask) {
        using Topology  = atlas::mesh::Nodes::Topology;
        using Range     = atlas::array::Range;
        auto node_flags = array::make_view<int, 1>(fvm.node_columns().nodes().flags());
        auto is_polar   = [&](idx_t j) {
            return Topology::check(node_flags(j), Topology::BC | Topology::NORTH) ||
                   Topology::check(node_flags(j), Topology::BC | Topology::SOUTH);
        };
        auto apply = [&](array::LocalView<double, 2>&& view) {
            for (idx_t j = 0; j < view.shape(0); ++j) {
                if (is_polar(j)) {
                    for (idx_t v = 0; v < view.shape(1); ++v) {
                        view(j, v) = mask;
                    }
                }
            }
        };
        if (field.rank() == 1) {
            apply(array::make_view<double, 1>(field).slice(Range::all(), Range::dummy()));
        }
        else if (field.rank() == 2) {
            apply(array::make_view<double, 2>(field).slice(Range::all(), Range::all()));
        }
    };

    for (auto fieldname : std::vector<std::string>{"dfdx", "dfdy", "div", "vor"}) {
        auto err_field  = fields.add(fvm.node_columns().createField<double>(option::name("err_" + fieldname)));
        auto err2_field = fields.add(fvm.node_columns().createField<double>(option::name("err2_" + fieldname)));
        auto fld        = array::make_view<double, 1>(fields[fieldname]);
        auto ref        = array::make_view<double, 1>(fields["ref_" + fieldname]);
        auto err        = array::make_view<double, 1>(fields["err_" + fieldname]);
        auto err2       = array::make_view<double, 1>(fields["err2_" + fieldname]);
        for (idx_t j = 0; j < fld.shape(0); ++j) {
            err(j)  = fld(j) - ref(j);
            err2(j) = err(j) * err(j);
        }

        if (mask_polar_values()) {
            do_mask_polar_values(fields["err_" + fieldname], 0.);
            do_mask_polar_values(fields["err2_" + fieldname], 0.);
        }
    }

    fields.haloExchange();


    // output to gmsh
    if (output_gmsh()) {
        output::Gmsh{"mesh_2d.msh", util::Config("coordinates", "lonlat")}.write(mesh);
        output::Gmsh{"mesh_3d.msh", util::Config("coordinates", "xyz")}.write(mesh);
        output::Gmsh{"fields.msh"}.write(fields);
    }

    auto minmax_within_error = [&](const std::string& name, double error) {
        Field field = fields["err_" + name];
        error       = std::abs(error);
        double min, max;
        fvm.node_columns().minimum(field, min);
        fvm.node_columns().maximum(field, max);
        bool success = true;
        if (min < -error) {
            Log::warning() << "minumum " << min << " smaller than error " << -error << std::endl;
            success = false;
        }
        if (max > error) {
            Log::warning() << "maximum " << max << " greater than error " << error << std::endl;
            success = false;
        }
        Log::info() << name << "\t: minmax error between { " << min << " , " << max << " }" << std::endl;
        return success;
    };
    EXPECT(minmax_within_error("dfdx", 1.e-11));
    EXPECT(minmax_within_error("dfdy", 1.e-11));
    EXPECT(minmax_within_error("div", 1.e-16));
    EXPECT(minmax_within_error("vor", 1.5e-9));

    auto rms_within_error = [&](const std::string& name, double error) {
        Field field = fields["err2_" + name];
        double mean;
        idx_t N;
        fvm.node_columns().mean(field, mean, N);
        double rms   = std::sqrt(mean / double(N));
        bool success = true;
        if (rms > error) {
            Log::warning() << "rms " << rms << " greater than error " << error << std::endl;
            success = false;
        }
        Log::info() << name << "\t: rms error = " << rms << std::endl;
        return success;
    };
    EXPECT(rms_within_error("dfdx", 1.e-14));
    EXPECT(rms_within_error("dfdy", 1.e-14));
    EXPECT(rms_within_error("div", 5.e-20));
    EXPECT(rms_within_error("vor", 5.e-13));

    // error for vorticity seems too high ?
}

//-----------------------------------------------------------------------------

}  // namespace test
}  // namespace atlas

int main(int argc, char** argv) {
    return atlas::test::run(argc, argv);
}
