/*
 * (C) Copyright 1996- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <cmath>
#include <fstream>
#include <map>
#include <unordered_map>

#include "eckit/geometry/Sphere.h"
#include "eckit/log/Bytes.h"
#include "eckit/log/JSON.h"
#include "eckit/types/FloatCompare.h"

#include "atlas/array.h"
#include "atlas/array/MakeView.h"
#include "atlas/field.h"
#include "atlas/grid.h"
#include "atlas/grid/Spacing.h"
#include "atlas/grid/detail/spacing/CustomSpacing.h"
#include "atlas/interpolation/Interpolation.h"
#include "atlas/interpolation/method/unstructured/ConservativeSphericalPolygonInterpolation.h"
#include "atlas/mesh.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/mesh/actions/Build2DCellCentres.h"
#include "atlas/meshgenerator.h"
#include "atlas/option.h"
#include "atlas/output/Gmsh.h"
#include "atlas/runtime/AtlasTool.h"
#include "atlas/util/Config.h"
#include "atlas/util/function/MDPI_functions.h"
#include "atlas/util/function/SolidBodyRotation.h"
#include "atlas/util/function/SphericalHarmonic.h"
#include "atlas/util/function/VortexRollup.h"

#include "tests/AtlasTestEnvironment.h"


namespace atlas {


class AtlasParallelInterpolation : public AtlasTool {
    int execute(const AtlasTool::Args& args) override;
    std::string briefDescription() override { return "Demonstration of parallel conservative interpolation"; }
    std::string usage() override {
        return name() +
               "[OPTIONS]... [--help]";
    }

    int numberOfPositionalArguments() override { return -1; }
    int minimumPositionalArguments() override { return 0; }

public:
    AtlasParallelInterpolation(int argc, char* argv[]): AtlasTool(argc, argv) {
        // Grid options
        add_option(new eckit::option::Separator("Grid options"));
        add_option(new SimpleOption<std::string>("source.grid", "source gridname"));
        add_option(new SimpleOption<std::string>("source.partitioner", "source partitioner name (spherical-polygon, lonlat-polygon, brute-force)"));
        add_option(new SimpleOption<std::string>("target.grid", "target gridname"));
        add_option(new SimpleOption<std::string>("target.partitioner", "target partitioner name (equal_regions, regular_bands, equal_bands)"));
        add_option(new SimpleOption<std::string>("source.functionspace",
                                                 "source functionspace, to override source grid default"));
        add_option(new SimpleOption<std::string>("target.functionspace",
                                                 "target functionspace, to override target grid default"));
        add_option(new SimpleOption<long>("source.halo", "default=2"));
        add_option(new SimpleOption<long>("target.halo", "default=0 for CellColumns and 1 for NodeColumns"));

        // Interpolation options
        add_option(new eckit::option::Separator("Interpolation options"));
        add_option(new SimpleOption<long>("order", "Interpolation order. Supported: 1, 2 (default=1)"));
        add_option(new SimpleOption<bool>("normalise_intersections",
                                          "Normalize polygon intersections so that interpolation weights sum to 1."));
        add_option(new SimpleOption<bool>("validate",
                                          "Enable extra validations at cost of performance. For debugging purpose."));
        add_option(new SimpleOption<bool>("matrix_free", "Do not store matrix for consecutive interpolations"));

        // Interpolation statistics options
        add_option(new eckit::option::Separator("Interpolation statistics options"));
        add_option(new SimpleOption<bool>("statistics.all", "Enable all statistics"));
        add_option(
            new SimpleOption<bool>("statistics.timings", "Enable statistics on interpolation setup and execution timings"));
        add_option(
            new SimpleOption<bool>("statistics.intersection", "Enable statistics on polygon intersections"));
        add_option(new SimpleOption<bool>("statistics.accuracy",
                                          "Enable statistics w.r.t. an analytical solution"));
        add_option(
            new SimpleOption<bool>("statistics.conservation", "Enable statistics on mass conservation"));

        // Output options
        add_option(new eckit::option::Separator("Output options"));
        add_option(new SimpleOption<bool>(
            "output-gmsh", "Output gmsh files src_mesh.msh, tgt_mesh.msh, src_field.msh, tgt_field.msh"));
        add_option(new SimpleOption<std::string>("gmsh.coordinates", "Mesh coordinates [xy,lonlat,xyz]"));
        add_option(new SimpleOption<bool>("gmsh.ghost", "output of ghost"));

        add_option(new SimpleOption<bool>("output-json", "Output json file with run information"));
        add_option(new SimpleOption<std::string>("json.file", "File path for json output"));

        // Initial condition options
        add_option(new eckit::option::Separator("Initial condition options"));
        add_option(new SimpleOption<std::string>(
            "init", "Setup initial source field [ constant, spherical_harmonic, vortex_rollup (default), solid_body_rotation_wind_magnitude ]"));
        add_option(new SimpleOption<double>("solid_body_rotation.angle", "Angle of solid body rotation (default = 0.)"));
        add_option(new SimpleOption<double>("vortex_rollup.t", "Value that controls vortex rollup (default = 0.5)"));
        add_option(new SimpleOption<double>("constant.value", "Value that is assigned in case init==constant)"));
        add_option(new SimpleOption<long>("spherical_harmonic.n", "total wave number 'n' of a spherical harmonic"));
        add_option(new SimpleOption<long>("spherical_harmonic.m", "zonal wave number 'm' of a spherical harmonic"));
    }

    struct Timers {
        using StopWatch = atlas::runtime::trace::StopWatch;
        StopWatch target_setup;
        StopWatch source_setup;
        StopWatch initial_condition;
        StopWatch interpolation_setup;
        StopWatch interpolation_execute;
    } timers;
};

std::function<double(const PointLonLat&)> get_init(const eckit::LocalConfiguration& args) {
    std::string init;
    args.get("init", init = "vortex_rollup");
    if (init == "vortex_rollup") {
        double t;
        args.get("vortex_rollup.t", t = 1.);
        return [t](const PointLonLat& p) { return util::function::vortex_rollup(p.lon(), p.lat(), t); };
    }
    else if (init == "spherical_harmonic") {
        int n = 2;
        int m = 2;
        args.get("spherical_harmonic.n", n);
        args.get("spherical_harmonic.m", m);

        bool caching = true;  // true -> warning not thread-safe
        util::function::SphericalHarmonic Y(n, m, caching);
        return [Y](const PointLonLat& p) { return Y(p.lon(), p.lat()); };
    }
    else if (init == "constant") {
        double value;
        args.get("constant.value", value = 1.);
        return [value](const PointLonLat&) { return value; };
    }
    else if (init == "solid_body_rotation_wind_magnitude") {
        double beta;
        args.get("solid_body_rotation.angle", beta = 0.);
        util::function::SolidBodyRotation sbr(beta);
        return [sbr](const PointLonLat& p) { return sbr.windMagnitude(p.lon(), p.lat()); };
    }
    else if (init == "MDPI_sinusoid") {
        auto sbr = util::function::MDPI_sinusoid;
        return [sbr](const PointLonLat& p) { return sbr(p.lon(), p.lat()); };
    }
    else if (init == "MDPI_harmonic") {
        auto sbr = util::function::MDPI_harmonic;
        return [sbr](const PointLonLat& p) { return sbr(p.lon(), p.lat()); };
    }
    else if (init == "MDPI_vortex") {
        auto sbr = util::function::MDPI_vortex;
        return [sbr](const PointLonLat& p) { return sbr(p.lon(), p.lat()); };
    }
    else if (init == "MDPI_gulfstream") {
        auto sbr = util::function::MDPI_gulfstream;
        return [sbr](const PointLonLat& p) { return sbr(p.lon(), p.lat()); };
    }
    else {
        if (args.has("init")) {
            Log::error() << "Bad value for \"init\": \"" << init << "\" not recognised." << std::endl;
            ATLAS_NOTIMPLEMENTED;
        }
    }
    ATLAS_THROW_EXCEPTION("Should not be here");
}

int AtlasParallelInterpolation::execute(const AtlasTool::Args& args) {\
    eckit::LocalConfiguration config(args);

    auto get_grid = [](std::string grid_name) {
        int grid_number = std::atoi( grid_name.substr(1, grid_name.size()).c_str() );
        if (grid_name.at(0) == 'P') {
            Log::info() << "P-grid number: " << grid_number << std::endl;
            ATLAS_ASSERT(grid_number > 3);
            std::vector<double> y = {90, 89.9999, 0, -90};
            auto xspace = StructuredGrid::XSpace( grid::LinearSpacing(0, 360, grid_number, false) );
            auto yspace = StructuredGrid::YSpace( new grid::spacing::CustomSpacing( y.size(), y.data() ) );
            return StructuredGrid(xspace, yspace);
        }
        else {
            return StructuredGrid{grid_name};
        }
    };
    auto src_grid = get_grid(config.getString("source.grid", "H32"));
    auto tgt_grid = get_grid(config.getString("target.grid", "H32"));

    auto create_functionspace = [&](Mesh& mesh, int halo, std::string type) -> FunctionSpace {
        if (type.empty()) {
            type = "NodeColumns";
            if (mesh.grid().type() == "healpix" || mesh.grid().type() == "cubedsphere") {
                type = "CellColumns";
            }
        }
        if (type == "CellColumns") {
            auto fspace = functionspace::CellColumns(mesh, option::halo(halo));
            if (! mesh.cells().has_field("lonlat")) {
                mesh::actions::Build2DCellCentres{"lonlat"}(mesh);
            }
            return fspace;
        }
        else if (type == "NodeColumns") {
            return functionspace::NodeColumns(mesh, option::halo(std::max(1,halo)));
        }
        ATLAS_THROW_EXCEPTION("FunctionSpace " << type << " is not recognized.");
    };

    timers.target_setup.start();
    auto tgt_mesh = Mesh{tgt_grid, grid::Partitioner(config.getString("target.partitioner", "regular_bands"))};
    auto tgt_functionspace =
        create_functionspace(tgt_mesh, config.getLong("target.halo", 0), config.getString("target.functionspace", ""));
    auto tgt_field = tgt_functionspace.createField<double>();
    timers.target_setup.stop();

    timers.source_setup.start();
    auto src_meshgenerator =
        MeshGenerator{src_grid.meshgenerator() | option::halo(2) | util::Config("pole_elements", "")};
    auto src_partitioner = grid::MatchingPartitioner{tgt_mesh, util::Config("partitioner", config.getString("source.partitioner", "spherical-polygon"))};
    auto src_mesh        = src_meshgenerator.generate(src_grid, src_partitioner);
    auto src_functionspace =
        create_functionspace(src_mesh, config.getLong("source.halo", 2), config.getString("source.functionspace", ""));
    auto src_field = src_functionspace.createField<double>();
    timers.source_setup.stop();

    {
        ATLAS_TRACE("Initial condition");
        timers.initial_condition.start();
        const auto lonlat = array::make_view<double, 2>(src_functionspace.lonlat());
        auto src_view     = array::make_view<double, 1>(src_field);
        auto f            = get_init(config);
        for (idx_t n = 0; n < lonlat.shape(0); ++n) {
            src_view(n) = f(PointLonLat{lonlat(n, LON), lonlat(n, LAT)});
        }
        src_field.set_dirty(true);
        timers.initial_condition.stop();
    }

    timers.interpolation_setup.start();
    auto interpolation =
        Interpolation(option::type("conservative-spherical-polygon") | config, src_functionspace, tgt_functionspace);
    Log::info() << interpolation << std::endl;
    timers.interpolation_setup.stop();

    timers.interpolation_execute.start();
    auto metadata = interpolation.execute(src_field, tgt_field);
    tgt_field.haloExchange();
    timers.interpolation_execute.stop();

    Field src_conservation_field;
    {
        using Statistics = interpolation::method::ConservativeSphericalPolygonInterpolation::Statistics;
        Statistics stats(metadata);
        if (config.getBool("statistics.accuracy", false) || config.getBool("statistics.all", false)) {
            stats.compute_accuracy(interpolation, tgt_field, get_init(config), &metadata);
        }
    }

    Log::info() << "interpolation metadata: \n";
    {
        eckit::JSON json(Log::info(), eckit::JSON::Formatting::indent(2));
        json << metadata;
    }
    Log::info() << std::endl;

    if (config.getBool("output-gmsh", false)) {
        if (config.getBool("gmsh.ghost", false)) {
            ATLAS_TRACE("halo exchange target");
            tgt_field.haloExchange();
        }
        util::Config gmsh_config(config.getSubConfiguration("gmsh"));
        output::Gmsh{"src_mesh.msh", gmsh_config}.write(src_mesh);
        output::Gmsh{"src_field.msh", gmsh_config}.write(src_field);
        output::Gmsh{"tgt_mesh.msh", gmsh_config}.write(tgt_mesh);
        output::Gmsh{"tgt_field.msh", gmsh_config}.write(tgt_field);
        if (src_conservation_field) {
            output::Gmsh{"src_conservation_field.msh", config}.write(src_conservation_field);
        }
    }

    if (config.getBool("output-json", false)) {
        util::Config output;
        output.set("setup.source.grid", config.getString("source.grid"));
        output.set("setup.target.grid", config.getString("target.grid"));
        output.set("setup.source.functionspace", src_functionspace.type());
        output.set("setup.target.functionspace", tgt_functionspace.type());
        output.set("setup.source.halo", config.getLong("source.halo", 2));
        output.set("setup.target.halo", config.getLong("target.halo", 0));
        output.set("setup.interpolation.order", config.getInt("order", 1));
        output.set("setup.interpolation.normalise_intersections", config.getBool("normalise_intersections", false));
        output.set("setup.interpolation.validate", config.getBool("validate", false));
        output.set("setup.interpolation.matrix_free", config.getBool("matrix-free", false));
        output.set("setup.init", config.getString("init", "vortex_rollup"));

        output.set("runtime.mpi", mpi::size());
        output.set("runtime.omp", atlas_omp_get_max_threads());
        output.set("atlas.build_type", ATLAS_BUILD_TYPE);

        output.set("timings.target.setup", timers.target_setup.elapsed());
        output.set("timings.source.setup", timers.source_setup.elapsed());
        output.set("timings.initial_condition", timers.initial_condition.elapsed());
        output.set("timings.interpolation.setup", timers.interpolation_setup.elapsed());
        output.set("timings.interpolation.execute", timers.interpolation_execute.elapsed());

        output.set("interpolation", metadata);

        eckit::PathName json_filepath(config.getString("json.file", "out.json"));
        std::ostringstream ss;
        eckit::JSON json(ss, eckit::JSON::Formatting::indent(4));
        json << output;

        eckit::FileStream file(json_filepath, "w");
        std::string str = ss.str();
        file.write(str.data(), str.size());
        file.close();
    }


    return success();
}

}  // namespace atlas


int main(int argc, char* argv[]) {
    atlas::AtlasParallelInterpolation tool(argc, argv);
    return tool.start();
}
