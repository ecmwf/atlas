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

using CSPInterpolation = interpolation::method::ConservativeSphericalPolygonInterpolation;

class AtlasEOAComputation : public AtlasTool {
    int execute(const AtlasTool::Args& args) override;
    std::string briefDescription() override { return "Computation of experimental order of convergence for Atlas' interpolation methods"; }
    std::string usage() override {
        return name() +
               " [--source.grid=gridname]"
               " [--target.grid=gridname] [Options]... [--help]";
    }

    int numberOfPositionalArguments() override { return -1; }
    int minimumPositionalArguments() override { return 0; }

public:
    AtlasEOAComputation(int argc, char* argv[]): AtlasTool(argc, argv) {
        add_option(new Separator("Options: PingPong"));
        add_option(new SimpleOption<long>("pingpong.nremaps","Number of remaps src->tgt->src"));

        add_option(new Separator("Options: Interpolation"));
        add_option(new SimpleOption<std::string>("type","Type of interpolation: bilinear, cubic, ..."));
        add_option(new SimpleOption<std::string>("order","Order of interpolation, when applicable"));
        add_option(new SimpleOption<std::string>("interpolation.structured","Is the interpolation for structured grids"));
        add_option(new SimpleOption<long>("k-nearest-neighbours", "The number of neighbours in k-nearest-neighbours"));
        
        add_option(new Separator("Options: Initial data"));
        add_option(new SimpleOption<bool>("init_via_highres", "Get initial data by remapping a highres grid data"));
        add_option(new SimpleOption<std::string>(
            "init", "Setup initial source field [ constant, spherical_harmonic, vortex_rollup (default), solid_body_rotation_wind_magnitude, MDPI_sinusoid, MDPI_harmonic, MDPI_vortex, MDPI_gulfstream ]"));
        add_option(new SimpleOption<double>("solid_body_rotation.angle", "Angle of solid body rotation (default = 0.)"));
        add_option(new SimpleOption<double>("vortex_rollup.t", "Value that controls vortex rollup (default = 0.5)"));
        add_option(new SimpleOption<double>("constant.value", "Value that is assigned in case init==constant)"));
        add_option(new SimpleOption<long>("spherical_harmonic.n", "total wave number 'n' of a spherical harmonic"));
        add_option(new SimpleOption<long>("spherical_harmonic.m", "zonal wave number 'm' of a spherical harmonic"));
    
        add_option(new Separator("Options: Input"));
        //add_option(new SimpleOption<std::string>("input.meshgenerator.pole_elements", "default = pentagons"));
        //add_option(new SimpleOption<bool>("input.meshed",  "Use function spaces based on mesh, also required for gmsh output"));
        //add_option(new SimpleOption<eckit::PathName>("input.gmsh.file", "Input gmsh file. If not provided, no gmsh output will be performed"));
        //add_option(new SimpleOption<std::string>("input.gmsh.coordinates", "Mesh coordinates: [xy, lonlat, xyz]"));

        add_option(new Separator("Options: Output"));
        add_option(new SimpleOption<bool>("output-gmsh", "Output gmsh file"));
        //add_option(new SimpleOption<bool>("output.cell_centred", "Overwrite defaults"));
        //add_option(new SimpleOption<std::string>("output.meshgenerator.pole_elements", "default = pentagons"));
        //add_option(new SimpleOption<bool>("output.meshed", "Use function spaces based on mesh, also required for gmsh output"));
        add_option(new SimpleOption<eckit::PathName>("output.gmsh.file", "Output gmsh file. If not provided, no gmsh output will be performed"));
        add_option(new SimpleOption<std::string>("output.gmsh.coordinates", "Mesh coordinates: [xy, lonlat, xyz]"));
       
        //add_option(new Separator("Optinos: Advanced / Debugging"));
        //add_option(new SimpleOption<bool>("double_remap", "default: true for target data on nodes"));
    }

    struct Timers {
        using StopWatch = atlas::runtime::trace::StopWatch;
        StopWatch target_setup;
        StopWatch source_setup;
        StopWatch initial_condition;
        StopWatch interpolation_fw_setup;
        StopWatch interpolation_bw_setup;
        StopWatch interpolation_fw_execute;
        StopWatch interpolation_bw_execute;
    } timers;
};

std::function<double(const PointLonLat&)> get_init(const AtlasTool::Args& args) {
    std::string init;
    args.get("init", init = "MDPI_vortex");
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
        util::function::SolidBodyRotation func(beta);
        return [func](const PointLonLat& p) { return func.windMagnitude(p.lon(), p.lat()); };
    }
    else if (init == "MDPI_sinusoid") {
        auto func = util::function::MDPI_sinusoid;
        return [func](const PointLonLat& p) { return func(p.lon(), p.lat()); };
    }
    else if (init == "MDPI_harmonic") {
        auto func = util::function::MDPI_harmonic;
        return [func](const PointLonLat& p) { return func(p.lon(), p.lat()); };
    }
    else if (init == "MDPI_vortex") {
        auto func = util::function::MDPI_vortex;
        return [func](const PointLonLat& p) { return func(p.lon(), p.lat()); };
    }
    else if (init == "MDPI_gulfstream") {
        auto func = util::function::MDPI_gulfstream;
        return [func](const PointLonLat& p) { return func(p.lon(), p.lat()); };
    }
    else {
        if (args.has("init")) {
            Log::error() << "Bad value for \"init\": \"" << init << "\" not recognised." << std::endl;
            ATLAS_NOTIMPLEMENTED;
        }
    }
    ATLAS_THROW_EXCEPTION("Should not be here");
}

FunctionSpace create_functionspace(Mesh& mesh, int halo, std::string type, bool struct_interpolation) {
    if (type.empty()) {
        type = "NodeColumns";
        if (mesh.grid().type() == "healpix" || mesh.grid().type() == "cubedsphere") {
            type = "CellColumns";
        }
    }
    if (struct_interpolation) {
        type = "StructuredColumns";
    }
    if (type == "CellColumns") {
        if (!mesh.cells().has_field("lonlat")) {
            mesh::actions::Build2DCellCentres{"lonlat"}(mesh);
        }
        return functionspace::CellColumns(mesh, option::halo(halo));
    }
    else if (type == "NodeColumns") {
        return functionspace::NodeColumns(mesh, option::halo(std::max(1,halo)));
    }
    else {
        return functionspace::StructuredColumns(mesh.grid(), option::halo(halo));
    }
    ATLAS_THROW_EXCEPTION("FunctionSpace " << type << " is not recognized.");
}


int AtlasEOAComputation::execute(const AtlasTool::Args& args) {
    ATLAS_ASSERT(atlas::mpi::size() == 1);

    std::stringstream sstream;
    int nremaps = args.getInt("pingpong.nremaps", 1);

    Grid src_grid = StructuredGrid( args.getString("source.grid", "O16") );
    Grid tgt_grid = StructuredGrid( args.getString("source.grid", "O32") );

    // setup interpolators
    //
    timers.source_setup.start();
    auto src_meshgenerator =
        MeshGenerator{src_grid.meshgenerator() | option::halo(2) | util::Config("pole_elements", "")};
    auto src_mesh        = src_meshgenerator.generate(src_grid);
    auto src_fs =
        create_functionspace(src_mesh, 4, args.getString("source.functionspace", ""), args.getBool("interpolation.structured", false));
    auto src_field = src_fs.createField<double>();
    timers.source_setup.stop();

    timers.initial_condition.start();
    const auto lonlat = array::make_view<double, 2>(src_fs.lonlat());
    auto src_view     = array::make_view<double, 1>(src_field);
    auto f            = get_init(args);
    for (idx_t n = 0; n < lonlat.shape(0); ++n) {
        src_view(n) = f(PointLonLat{lonlat(n, LON), lonlat(n, LAT)});
    }
    src_field.set_dirty(true);
    src_field.haloExchange();
    timers.initial_condition.stop();

    timers.target_setup.start();
    auto tgt_meshgenerator =
        MeshGenerator{tgt_grid.meshgenerator() | option::halo(2) | util::Config("pole_elements", "")};
    auto tgt_mesh        = tgt_meshgenerator.generate(tgt_grid);
    auto tgt_fs =
        create_functionspace(tgt_mesh, 1, args.getString("target.functionspace", ""), args.getBool("interpolation.structured", false));
    auto tgt_field = tgt_fs.createField<double>();
    timers.target_setup.stop();

    timers.interpolation_fw_setup.start();
    auto interpolation_fw =
        Interpolation(args, src_fs, tgt_fs);
    timers.interpolation_fw_setup.stop();
    timers.interpolation_bw_setup.start();
    auto interpolation_bw =
        Interpolation(args, tgt_fs, src_fs);
    timers.interpolation_bw_setup.stop();

    if (args.getBool("output-gmsh", false)) {
        if (args.getBool("gmsh.ghost", false)) {
            ATLAS_TRACE("halo exchange target");
            src_field.haloExchange();
        }
        util::Config config(args.getSubConfiguration("gmsh"));
        output::Gmsh{src_grid.name()+"_src.msh", config}.write(src_mesh);
        output::Gmsh{src_grid.name()+"_srcfield_init.msh", config}.write(src_field);
    }

    for (int irepeat = 0; irepeat < nremaps; irepeat++) {
        std::cout << "== remap step " << irepeat + 1 << " / " << nremaps << std::endl;
        timers.interpolation_fw_execute.start();
        interpolation_fw.execute(src_field, tgt_field);
        timers.interpolation_fw_execute.stop();
        timers.interpolation_bw_execute.start();
        interpolation_bw.execute(tgt_field, src_field);
        timers.interpolation_bw_execute.stop();
    }

    if (args.getBool("output-gmsh", false)) {
        if (args.getBool("gmsh.ghost", false)) {
            ATLAS_TRACE("halo exchange target");
            tgt_field.haloExchange();
        }
        util::Config config(args.getSubConfiguration("gmsh"));
        output::Gmsh{tgt_grid.name()+"_tgt.msh", config}.write(tgt_mesh);
        output::Gmsh{tgt_grid.name()+"_tgtfield_from_"+src_grid.name()+"_after_"+std::to_string(nremaps)+"remaps.msh", config}.write(tgt_field);
    }

    if (args.getBool("output-json", true)) {
        util::Config output;
        output.set("setup.source.grid", src_grid.name());
        output.set("setup.target.grid", tgt_grid.name());
        output.set("setup.source.functionspace", src_fs.type());
        output.set("setup.target.functionspace", tgt_fs.type());
        output.set("setup.source.halo", args.getLong("source.halo", 2));
        output.set("setup.target.halo", args.getLong("target.halo", 0));
        output.set("setup.interpolation.type", args.getInt("order", 1));
        output.set("setup.interpolation.order", args.getInt("order", 1));
        output.set("setup.interpolation.normalise_intersections", args.getBool("normalise-intersections", false));
        output.set("setup.interpolation.validate", args.getBool("validate", false));
        output.set("setup.interpolation.matrix_free", args.getBool("matrix-free", false));
        output.set("setup.init", args.getString("init", "vortex_rollup"));

        output.set("runtime.mpi", mpi::size());
        output.set("runtime.omp", atlas_omp_get_max_threads());
        output.set("atlas.build_type", ATLAS_BUILD_TYPE);

        output.set("timings.target.setup", timers.target_setup.elapsed());
        output.set("timings.source.setup", timers.source_setup.elapsed());
        output.set("timings.initial_condition", timers.initial_condition.elapsed());
        output.set("timings.interpolation_fw.setup", timers.interpolation_fw_setup.elapsed());
        output.set("timings.interpolation_bw.setup", timers.interpolation_bw_setup.elapsed());
        output.set("timings.interpolation_fw.execute", timers.interpolation_fw_execute.elapsed());
        output.set("timings.interpolation_bw.execute", timers.interpolation_bw_execute.elapsed());

        sstream.str("");
        sstream << "json_" << src_grid.name() << "_" << tgt_grid.name();
        std::cout << "output to: " <<  sstream.str() << std::endl;
        eckit::PathName json_filepath(sstream.str());
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
    atlas::AtlasEOAComputation tool(argc, argv);
    return tool.start();
}
