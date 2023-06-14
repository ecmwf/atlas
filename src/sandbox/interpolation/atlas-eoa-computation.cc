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
    std::string briefDescription() override { return "Computation of experimental order of accuracy for Atlas' interpolation methods"; }
    std::string usage() override {
        return name() +
               " [--source.grid=gridname] "
               "[--target.grid=gridname] [OPTION]... [--help]";
    }

    int numberOfPositionalArguments() override { return -1; }
    int minimumPositionalArguments() override { return 0; }

public:
    AtlasEOAComputation(int argc, char* argv[]): AtlasTool(argc, argv) {
        add_option(new Separator("EOC options"));
        add_option(new SimpleOption<long>("eoc.grid.cycles","Simulation cycles for determining EOC"));
        add_option(new SimpleOption<long>("eoc.grid.maxres","Maximal grid resolution in computing EOC"));
        add_option(new SimpleOption<std::string>("source.grid_type", "Source grid type: O, N, L, ..."));
        add_option(new SimpleOption<std::string>("target.grid_type", "Target grid type: O, N, L, ..."));
        add_option(new SimpleOption<long>("eoc.refine-source","Refine source to compute EOC, otherwise refine target"));

        add_option(new Separator("Interpolation options"));
        add_option(new SimpleOption<std::string>("type","Type of interpolation: bilinear, cubic, ..."));
        add_option(new SimpleOption<std::string>("order","Order of interpolation, when applicable"));
        add_option(new SimpleOption<std::string>("interpolation.structured","Is the interpolation for structured grids"));
        
        add_option(new Separator("Initial data"));
        add_option(new SimpleOption<std::string>(
            "init", "Setup initial source field [ constant, spherical_harmonic, vortex_rollup (default), solid_body_rotation_wind_magnitude, MDPI_sinusoid, MDPI_harmonic, MDPI_vortex, MDPI_gulfstream ]"));
        add_option(new SimpleOption<double>("solid_body_rotation.angle", "Angle of solid body rotation (default = 0.)"));
        add_option(new SimpleOption<double>("vortex_rollup.t", "Value that controls vortex rollup (default = 0.5)"));
        add_option(new SimpleOption<double>("constant.value", "Value that is assigned in case init==constant)"));
        add_option(new SimpleOption<long>("spherical_harmonic.n", "total wave number 'n' of a spherical harmonic"));
        add_option(new SimpleOption<long>("spherical_harmonic.m", "zonal wave number 'm' of a spherical harmonic"));
    
        add_option(new Separator("Input options"));
        //add_option(new SimpleOption<std::string>("input.meshgenerator.pole_elements", "default = pentagons"));
        //add_option(new SimpleOption<bool>("input.meshed",  "Use function spaces based on mesh, also required for gmsh output"));
        //add_option(new SimpleOption<eckit::PathName>("input.gmsh.file", "Input gmsh file. If not provided, no gmsh output will be performed"));
        //add_option(new SimpleOption<std::string>("input.gmsh.coordinates", "Mesh coordinates: [xy, lonlat, xyz]"));

        add_option(new Separator("Output options"));
        add_option(new SimpleOption<bool>("output-gmsh", "Output gmsh file"));
        //add_option(new SimpleOption<bool>("output.cell_centred", "Overwrite defaults"));
        //add_option(new SimpleOption<std::string>("output.meshgenerator.pole_elements", "default = pentagons"));
        //add_option(new SimpleOption<bool>("output.meshed", "Use function spaces based on mesh, also required for gmsh output"));
        //add_option(new SimpleOption<eckit::PathName>("output.gmsh.file", "Output gmsh file. If not provided, no gmsh output will be performed"));
        //add_option(new SimpleOption<std::string>("output.gmsh.coordinates", "Mesh coordinates: [xy, lonlat, xyz]"));
       
        add_option(new Separator("Advanced / Debugging"));
        //add_option(new SimpleOption<bool>("double_remap", "default: true for target data on nodes"));
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

std::function<double(const PointLonLat&)> get_init(const AtlasTool::Args& args) {
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


void compute_errors(const Field source, const Field target,
        std::function<double(const PointLonLat&)> func,
        Mesh src_mesh, Mesh tgt_mesh) {
    auto src_vals             = array::make_view<double, 1>(source);
    auto tgt_vals             = array::make_view<double, 1>(target);
    const auto src_node_ghost = array::make_view<int, 1>(src_mesh.nodes().ghost());
    const auto src_node_halo  = array::make_view<int, 1>(src_mesh.nodes().halo());
    const auto tgt_node_ghost = array::make_view<int, 1>(tgt_mesh.nodes().ghost());
    const auto tgt_node_halo  = array::make_view<int, 1>(tgt_mesh.nodes().halo());

    // get tgt_points and tgt_areas
    // TODO: no need for polygon intersections here, we just need src_points and src_areas
    auto src_fs = create_functionspace(src_mesh, 2, "NodeColumns", 0);
    auto tgt_fs = create_functionspace(tgt_mesh, 2, "NodeColumns", 0);
    auto interpolation = CSPInterpolation();
    interpolation.do_setup(src_fs, tgt_fs);

    std::cout << "source field size             : " << source.size() << std::endl;
    std::cout << "target field size             : " << target.size() << std::endl;
    std::cout << "src_mesh nodes                : " << src_mesh.nodes().size() << std::endl;
    std::cout << "tgt_mesh nodes                : " << tgt_mesh.nodes().size() << std::endl;
    std::cout << "CSP-interpolation.src_npoints : " << interpolation.src_npoints() << std::endl;
    std::cout << "CSP-interpolation.tgt_npoints : " << interpolation.tgt_npoints() << std::endl;

    // compute error to the analytical solution on the target
    double tgt_mass_pos = 0.;
    double tgt_mass_neg = 0.;
    double terr_remap_l2 = 0;
    double terr_remap_linf = -1.;
    int cc = 0;
    int ncc = 0;
    for (idx_t tpt = 0; tpt < interpolation.tgt_npoints(); ++tpt) {
        if (tgt_node_ghost(tpt) or tgt_node_halo(tpt)) {
            cc++;
            continue;
        }
        ncc++;
        if (tgt_vals(tpt) > 0.) {
            tgt_mass_pos += tgt_vals(tpt) * interpolation.tgt_areas(tpt);
        }
        else {
            tgt_mass_neg -= tgt_vals(tpt) * interpolation.tgt_areas(tpt);
        }
        auto p = interpolation.tgt_points(tpt);
        PointLonLat pll;
        eckit::geometry::Sphere::convertCartesianToSpherical(1., p, pll);
        double err_l = std::abs(tgt_vals(tpt) - func(pll));
        terr_remap_l2 += err_l * err_l * interpolation.tgt_areas(tpt);
        terr_remap_linf = std::max(terr_remap_linf, err_l);
    }
    std::cout << "tgt points (omitted), (considered) : " << cc << ", " << ncc << std::endl;

    cc = 0;
    ncc = 0;
    // compute the conservation error and the projection errors serr_remap_*
    double src_mass_pos = 0.;
    double src_mass_neg = 0.;
    double serr_remap_l2 = 0;
    double serr_remap_linf = -1.;
    for (idx_t spt = 0; spt < interpolation.src_npoints(); ++spt) {
        if (src_node_ghost(spt) or src_node_halo(spt)) {
            cc++;
            continue;
        }
        ncc++;
        if (src_vals(spt) > 0.) {
            src_mass_pos += src_vals(spt) * interpolation.src_areas(spt);
        }
        else {
            src_mass_neg -= src_vals(spt) * interpolation.src_areas(spt);
        }
        auto p = interpolation.src_points(spt);
        PointLonLat pll;
        eckit::geometry::Sphere::convertCartesianToSpherical(1., p, pll);
        double err_l = std::abs(src_vals(spt) - func(pll));
        serr_remap_l2 += err_l * err_l * interpolation.src_areas(spt);
        serr_remap_linf = std::max(serr_remap_linf, err_l);
    }
    std::cout << "src points (omitted), (considered) : " << cc << ", " << ncc << std::endl;

    serr_remap_l2 = std::sqrt(serr_remap_l2);
    terr_remap_l2 = std::sqrt(terr_remap_l2);
    double err_cons_pos = tgt_mass_pos - src_mass_pos;
    double err_cons_neg = tgt_mass_neg - src_mass_neg;
    double err_cons = err_cons_pos - err_cons_neg;
    double src_mass = src_mass_pos + src_mass_neg;
    double tgt_mass = tgt_mass_pos + tgt_mass_neg;
    std::cout << "src l2 error            : " << serr_remap_l2 << std::endl;
    std::cout << "src l_inf error         : " << serr_remap_linf << std::endl;
    std::cout << "tgt l2 error            : " << terr_remap_l2 << std::endl;
    std::cout << "tgt l_inf error         : " << terr_remap_linf << std::endl;
    std::cout << "conservation error      : " << err_cons << std::endl;
    std::cout << "total tgt mas           : " << tgt_mass << std::endl;
    std::cout << "total src mas           : " << src_mass << std::endl;
    std::cout << "rel. cons. error on src : " << 100. * err_cons/src_mass << " %" << std::endl;
    std::cout << "rel. cons. error on tgt : " << 100. * err_cons/tgt_mass << " %" << std::endl;
}


int AtlasEOAComputation::execute(const AtlasTool::Args& args) {
    ATLAS_ASSERT(atlas::mpi::size() == 1);

    std::stringstream sstream;
    int eoc_cycles = args.getInt("eoc.grid.cycles", 3); 
    int eoc_maxres = args.getInt("eoc.grid.maxres", 128); 
    std::string sgrid_type = args.getString("source.grid_type", "O");
    std::string tgrid_type = args.getString("target.grid_type", "O");
    bool refine_source = args.getBool("eoc.refine-source", true);

    Grid src_grid;
    Grid tgt_grid;
    if (refine_source) {
        sstream << tgrid_type << eoc_maxres;
        tgt_grid = StructuredGrid(sstream.str());
    }
    else {
        int gres = eoc_maxres;
        int cycles = eoc_cycles;
        for (; cycles--; gres/=2);
        sstream << sgrid_type << gres;
        src_grid = StructuredGrid(sstream.str());
    }

    for (int gres = eoc_maxres; eoc_cycles--; gres/=2) {
        sstream.str("");
        if (refine_source) {
            sstream << sgrid_type << gres;
            src_grid = StructuredGrid(sstream.str());
        }
        else {
            sstream << tgrid_type << gres;
            tgt_grid = StructuredGrid(sstream.str());
        }
        std::cout << src_grid.name() << " --> " << tgt_grid.name() << std::endl;

        timers.target_setup.start();
        auto tgt_mesh = Mesh{tgt_grid, grid::Partitioner(args.getString("target.partitioner", "serial"))};
        auto tgt_functionspace =
            create_functionspace(tgt_mesh, 2, args.getString("target.functionspace", ""), args.getBool("interpolation.structured", false));
        auto tgt_field = tgt_functionspace.createField<double>();
        timers.target_setup.stop();

        timers.source_setup.start();
        auto src_meshgenerator =
            MeshGenerator{src_grid.meshgenerator() | option::halo(2) | util::Config("pole_elements", "")};
        auto src_partitioner = grid::MatchingPartitioner{tgt_mesh, util::Config("partitioner",args.getString("source.partitioner", "spherical-polygon"))};
        auto src_mesh        = src_meshgenerator.generate(src_grid, src_partitioner);
        auto src_functionspace =
            create_functionspace(src_mesh, 2, args.getString("source.functionspace", ""), args.getBool("interpolation.structured", false));
        auto src_field = src_functionspace.createField<double>();
        timers.source_setup.stop();

        {
            ATLAS_TRACE("Initial condition");
            timers.initial_condition.start();
            const auto lonlat = array::make_view<double, 2>(src_functionspace.lonlat());
            auto src_view     = array::make_view<double, 1>(src_field);
            auto f            = get_init(args);
            for (idx_t n = 0; n < lonlat.shape(0); ++n) {
                src_view(n) = f(PointLonLat{lonlat(n, LON), lonlat(n, LAT)});
            }
            src_field.set_dirty(true);
            timers.initial_condition.start();
        }

        timers.interpolation_setup.start();
        auto interpolation =
            Interpolation(args, src_functionspace, tgt_functionspace);
        //Log::info() << interpolation << std::endl;
        timers.interpolation_setup.stop();


        timers.interpolation_execute.start();
        auto metadata = interpolation.execute(src_field, tgt_field);
        timers.interpolation_execute.stop();

        compute_errors(src_field, tgt_field, get_init(args), src_mesh, tgt_mesh);

        // API not yet acceptable
        Field src_conservation_field;
        {
            //using Statistics = interpolation::method::ConservativeSphericalPolygonInterpolation::Statistics;
            //Statistics stats(metadata);
            //stats.accuracy(interpolation, tgt_field, get_init(args), &metadata);
            //src_conservation_field = stats.diff(interpolation, src_field, tgt_field);
            //src_conservation_field.set_dirty(true);
            //src_conservation_field.haloExchange();
        }

        // skip metadata, most interpolation methods do not provide them anyways
        /*
        Log::info() << "interpolation metadata: \n";
        {
            eckit::JSON json(Log::info(), eckit::JSON::Formatting::indent(2));
            json << metadata;
        }
        Log::info() << std::endl;
        */

        if (args.getBool("output-gmsh", false)) {
            if (args.getBool("gmsh.ghost", false)) {
                ATLAS_TRACE("halo exchange target");
                tgt_field.haloExchange();
            }
            util::Config config(args.getSubConfiguration("gmsh"));
            sstream.str(src_grid.name());
            output::Gmsh{sstream.str()+".msh", config}.write(src_mesh);
            output::Gmsh{sstream.str()+"_field.msh", config}.write(src_field);
            sstream.str(tgt_grid.name());
            output::Gmsh{sstream.str()+".msh", config}.write(tgt_mesh);
            output::Gmsh{sstream.str()+"_field.msh", config}.write(tgt_field);
            if (src_conservation_field) {
                output::Gmsh{"src_conservation_field.msh", config}.write(src_conservation_field);
            }
        }

        if (args.getBool("output-json", true)) {
            util::Config output;
            output.set("setup.source.grid", src_grid.name());
            output.set("setup.target.grid", tgt_grid.name());
            output.set("setup.source.functionspace", src_functionspace.type());
            output.set("setup.target.functionspace", tgt_functionspace.type());
            output.set("setup.source.halo", args.getLong("source.halo", 2));
            output.set("setup.target.halo", args.getLong("target.halo", 0));
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
            output.set("timings.interpolation.setup", timers.interpolation_setup.elapsed());
            output.set("timings.interpolation.execute", timers.interpolation_execute.elapsed());

            output.set("interpolation", metadata);

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
    }

    return success();
}

}  // namespace atlas


int main(int argc, char* argv[]) {
    atlas::AtlasEOAComputation tool(argc, argv);
    return tool.start();
}
