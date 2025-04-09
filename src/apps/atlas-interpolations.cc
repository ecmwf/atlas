/*
 * (C) Copyright 2025 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <istream>

#include "atlas/array.h"
#include "atlas/field/Field.h"
#include "atlas/field/FieldSet.h"
#include "atlas/functionspace/CellColumns.h"
#include "atlas/functionspace/NodeColumns.h"
#include "atlas/functionspace/PointCloud.h"
#include "atlas/functionspace/StructuredColumns.h"
#include "atlas/grid.h"
#include "atlas/interpolation.h"
#include "atlas/interpolation/AssembleGlobalMatrix.h"
#include "atlas/interpolation/method/MethodFactory.h"
#include "atlas/linalg/sparse.h"
#include "atlas/linalg/sparse/MakeEckitSparseMatrix.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/meshgenerator.h"
#include "atlas/output/Gmsh.h"
#include "atlas/output/ScripIO.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/runtime/AtlasTool.h"
#include "atlas/util/CoordinateEnums.h"
#include "atlas/util/function/VortexRollup.h"

#include "tests/AtlasTestEnvironment.h"

using atlas::functionspace::PointCloud;
using atlas::functionspace::NodeColumns;
using atlas::functionspace::StructuredColumns;
using atlas::util::Config;

using Matrix      = atlas::linalg::SparseMatrixStorage;


namespace atlas {


// This app can:
//
// 1) list all available interpolators
//      atlas-interpolations --list
// 2) list a limited set of available interpolators for a given source and a target grids
//      atlas-interpolations --sgrid <sgrid> --tgrid <tgrid> --list
// 3) compute and store the global interpolation matrix to disk:
//      atlas-interpolations --sgrid <sgrid> --tgrid <tgrid> --interpolation <intrp> --format <fmt>
// 4) read in a stored global interpolation matrix from the disc:
//      atlas-interpolations --sgrid <sgrid> --tgrid <tgrid> --format <fmt> --read --matrix <matrix>
//
// Formats <fmt> can be 'eckit' (binary) or 'scrip' (netcdf).
//
// This app hides configurations of meshes, function spaces, halo size, etc. from the user.
//

class AtlasInterpolations : public AtlasTool {
    int execute(const AtlasTool::Args& args) override;
    std::string briefDescription() override {
        return "Atlas intepolation app to perform MPI-parallel interpolations for given source "
                "and target grids and/or generate remapping matrices.\n"; }
    std::string usage() override {
        return name() + " [OPTION] ... [--help]";
    }
    std::string longDescription() override {
        return "The 'sgrid' and 'tgrid' arguments are source and target grids name, "
               "for instance: N80, F40, O24, L64x33, CS-ED-12, etc. See './bin/atlas-grids --list' "
               "for the complete list of available grids. See './bin/atlas-interpolation --list' "
               "for a list of named interpolations. See './bin/atlas-interpolation --sgrid <grid1> "
               "--tgrid <grid2> --list' for the list of available interpolations for that grid pair.\n"
               "\n";
    }

    int numberOfPositionalArguments() override { return -1; }
    int minimumPositionalArguments() override { return 0; }
    void interpolation_config(Config&);
    void create_fspaces(Config&, const Grid& input_grid, const Grid& output_grid,
        FunctionSpace& fs_in, FunctionSpace& fs_out, grid::Distribution& dist_in, grid::Distribution& dist_out);
    Matrix read_matrix(std::string matrix_name, std::string matrix_format);
    void write_matrix(const Matrix& mat, std::string matrix_name, std::string matrix_format);
    std::string get_matrix_name(std::string& sgrid, std::string& tgrid, std::string& interp, std::string format = std::string());
    bool init_interpolation_database();
    bool available_combination(Grid sgrid, Grid tgrid, std::string interpolation);
    bool grid_type_same(const char* c1, const char* c2);
    bool get_available_interpolations(Grid sgrid, Grid tgrid, std::vector<std::string>& serial_interpolations,
        std::vector<std::string>& parallel_interpolations);

    // there is a file database of functioning interpolation combinations which 
    // is read into the member idatabase_ via InterpolationStatus class
    struct InterpolationStatus {
        std::string mpi_tasks;
        std::string interpolation;
        std::string sgrid;
        std::string tgrid;
        std::string status;

        friend std::istream& operator>>(std::istream& in, InterpolationStatus& istat) {
            std::array<std::string, 7> line;
            in >> line[0] >> line[1] >> line[2] >> line[3] >> line[4] >> line[5] >> line[6];
            istat.mpi_tasks = line[0];
            istat.interpolation = line[2];
            istat.sgrid = line[3];
            istat.tgrid = line[4];
            istat.status = line[6];
            return in;
        }
    };
    std::vector<InterpolationStatus> idatabase_;

    struct Timers {
        using StopWatch = atlas::runtime::trace::StopWatch;
        StopWatch functionspace_setup;
        StopWatch interpolation_setup;
        StopWatch interpolation_exe;
        StopWatch global_matrix_read;
        StopWatch global_matrix_setup;
        StopWatch global_matrix_exe;
    } timers;

public:
    AtlasInterpolations(int argc, char* argv[]): AtlasTool(argc, argv) {
        add_option(new SimpleOption<bool>("list", "list available interpolations"));
        add_option(new SimpleOption<bool>("force", "do not check if the interpolator is available"));
        add_option(new SimpleOption<std::string>("sgrid", "source grid"));
        add_option(new SimpleOption<bool>("scell", "switch to cell-centred data on the source grid (only for interpolation = conservative-spherical-polygon)"));
        add_option(new SimpleOption<std::string>("tgrid", "target grid"));
        add_option(new SimpleOption<bool>("tcell", "switch to cell-centred data on the target grid (only for interpolation = conservative-spherical-polygon)"));
        add_option(new SimpleOption<std::string>("interpolation", "interpolation methods"));
        add_option(new SimpleOption<bool>("read", "use a provided remapping matrix"));
        add_option(new SimpleOption<std::string>("matrix", "name of the remapping matrix"));
        add_option(new SimpleOption<std::string>("format", "format of the remapping matrix: eckit, SCRIP"));
        add_option(new SimpleOption<bool>("test", "write the interpolated result from a test source field"));
        add_option(new SimpleOption<bool>("partest", "write the interpolated result from a test source field with the parallel interpolator"));
        add_option(new SimpleOption<bool>("3d", "write the interpolated result for visualisation on the sphere"));
    }
};


int AtlasInterpolations::execute(const AtlasTool::Args& args) {
    ATLAS_TRACE("main");
    if (args.has("force")) {
        Log::info() << "\nWARNING Running an untested combination of interpolators and grids.\n\n";
    }
    std::string sgrid_name = "O8";
    args.get("sgrid", sgrid_name);
    std::string tgrid_name = "O32";
    args.get("tgrid", tgrid_name);
    auto sgrid = Grid{sgrid_name};
    auto tgrid = Grid{tgrid_name};
    std::string sdata_name = "node data";
    if (args.has("scell")) {
        sdata_name = "cell data";
    }
    std::string tdata_name = "node data";
    if (args.has("tcell")) {
        tdata_name = "cell data";
    }

    if ((args.has("sgrid") || args.has("tgrid")) && args.has("list")) {
        std::vector<std::string> serial_interpolations;
        std::vector<std::string> parallel_interpolations;
        if (not get_available_interpolations(sgrid, tgrid, serial_interpolations, parallel_interpolations)) {
            Log::info() << "Please generate an Atlas interpolation database." <<std::endl;
            return failed();
        };

        auto list_tabbed = [](std::vector<std::string>& l) {
            if (l.size()) {
                for (auto i : l) {
                    Log::info() << "\t" << i << std::endl;
                    if (i == "conservative-spherical-polygon") {
                        Log::info() << "\t" << "conservative-spherical-polygon-2" << '\n';
                    }
                }
            }
            else {
                Log::info() << "\t no interpolations available" << std::endl;
            }
        };

        Log::info() << "Available serial interpolations for sgrid = " << sgrid_name << " and " << "tgrid = " << tgrid_name << ":" << std::endl;
        list_tabbed(serial_interpolations);
        Log::info() << "Available parallel interpolations for sgrid = " << sgrid_name << " and " << "tgrid = " << tgrid_name << ":" << std::endl;
        list_tabbed(parallel_interpolations);
        return success();
    }
    else if (args.has("list")) {
        Log::info() << "Usage: " << usage() << std::endl;
        Log::info() << "Available interpolators:" << std::endl;
        std::set<std::string> interpolation_types;
        for (const auto& key : interpolation::MethodFactory::keys()) {
            interpolation_types.insert(key);
        }
        for (const auto& b : interpolation_types) {
            Log::info() << "\t" << b << '\n';
        }
        Log::info() << "\t" << "conservative-spherical-polygon-2" << '\n';
        return success();
    }
    std::string interpolation_method = "finite-element";
    args.get("interpolation", interpolation_method);
    Log::info() << "\tsource grid\t:\t" << sgrid_name << " (" << sdata_name <<")\n\ttarget grid\t:\t" << tgrid_name
        << " (" << tdata_name <<")\n\tinterpolation\t:\t" << interpolation_method << std::endl;
    Config config;
    config.set("src_cell_data", args.has("scell"));
    config.set("tgt_cell_data", args.has("tcell"));
    config.set("interpolation_method", interpolation_method);
    std::string matrix_format = "eckit";
    args.get("format", matrix_format);

    if (args.has("read")) {
        ATLAS_TRACE("Setup");
        Matrix gmatrix;
        std::string matrix_name("");
        args.get("matrix", matrix_name);
        if (matrix_name == "") {
            matrix_name = get_matrix_name(sgrid_name, tgrid_name, interpolation_method);
        }
        FunctionSpace src_fs;
        FunctionSpace tgt_fs;
        grid::Distribution src_dist, tgt_dist;
        Interpolation interpolator;

        ATLAS_TRACE_SCOPE("Setup source and target") {
            timers.functionspace_setup.start();
            if (interpolation_method.find("nearest-neighbour") != std::string::npos) {
                interpolation_method = "finite-element";
            }
            create_fspaces(config, sgrid, tgrid, src_fs, tgt_fs, src_dist, tgt_dist);
            timers.functionspace_setup.stop();
        }
        ATLAS_TRACE_SCOPE("Setup interpolator") {
            ATLAS_TRACE_SCOPE("Read matrix on rank 0")
            if (mpi::comm().rank() == 0) {
                timers.global_matrix_read.start();
                gmatrix = read_matrix(matrix_name, matrix_format);
                timers.global_matrix_read.stop();
                Log::info() << "Global matrix read timer\t: " << timers.global_matrix_read.elapsed() * 1000. << " [ms]"  << std::endl;
                auto eckit_gmatrix = atlas::linalg::make_non_owning_eckit_sparse_matrix(gmatrix);
                eckit_gmatrix.print(Log::info());
                Log::info() << std::endl;
            }
            mpi::comm().barrier();

            timers.global_matrix_setup.start();
            auto matrix = interpolation::distribute_global_matrix(tgt_dist, src_fs, tgt_fs, gmatrix);
            timers.global_matrix_setup.stop();

            timers.interpolation_setup.start();
            interpolation::MatrixCache cache(std::move(matrix));
            interpolator = Interpolation(config, src_fs, tgt_fs, cache);
            timers.interpolation_setup.stop();
        }

        double timer_functionspace_setup = timers.functionspace_setup.elapsed() / mpi::comm().size();
        double timer_interpolation_setup = timers.interpolation_setup.elapsed() / mpi::comm().size();
        double timer_global_matrix_setup = timers.global_matrix_setup.elapsed() / mpi::comm().size();
        mpi::comm().allReduceInPlace(&timer_functionspace_setup, 1, eckit::mpi::sum());
        mpi::comm().allReduceInPlace(&timer_interpolation_setup, 1, eckit::mpi::sum());
        mpi::comm().allReduceInPlace(&timer_global_matrix_setup, 1, eckit::mpi::sum());
        Log::info() << "Grid + FunctionSpace timer\t: " << timer_functionspace_setup * 1000. << " [ms]"  << std::endl;
        Log::info() << "Interpolation setup timer\t: " << timer_interpolation_setup * 1000. << " [ms]"  << std::endl;
        Log::info() << "Global matrix setup timer\t: " << timer_global_matrix_setup * 1000. << " [ms]"  << std::endl;
    
        if (args.has("test")) {
            ATLAS_TRACE("test");
            // Allocate and initialise own memory here to show possibilities
            // Note: reading a field from disc is an extra feature
            auto src_field = interpolator.source().createField<double>();
            auto tgt_field = interpolator.target().createField<double>();
            auto src_lonlat = array::make_view<double, 2>(interpolator.source().lonlat());
            ATLAS_TRACE_SCOPE("initialize source") {
                auto src_field_v = array::make_view<double, 1>(src_field);
                for (idx_t i = 0; i < src_fs.size(); ++i) {
                    src_field_v[i] = util::function::vortex_rollup(src_lonlat(i, 0), src_lonlat(i, 1), 1.);
                }
            }

            timers.interpolation_exe.start();
            interpolator.execute(src_field, tgt_field);
            timers.interpolation_exe.stop();

            double timer_interpolation_exe = timers.interpolation_exe.elapsed() / mpi::comm().size();
            mpi::comm().allReduceInPlace(&timer_interpolation_exe, 1, eckit::mpi::sum());
            Log::info() << "Interpolation execute timer     : " << timer_interpolation_exe * 1000. << " [ms]"  << std::endl;

            ATLAS_TRACE_SCOPE("output from the read-in matrix") {
                std::string tgt_name = "tfield_" + matrix_name;
                std::string coords = (args.has("3d") ? "3d" : "lonlat");
                output::Gmsh gmsh(tgt_name + ".msh", Config("coordinates", coords) | Config("ghost", "true"));
                if( functionspace::NodeColumns(tgt_field.functionspace())) {
                    Log::info() << "storing distributed remapped field '" << tgt_name << ".msh'." << std::endl;
                    gmsh.write(functionspace::NodeColumns(tgt_field.functionspace()).mesh());
                    tgt_field.haloExchange();
                    gmsh.write(tgt_field);
                }
                else if( functionspace::CellColumns(tgt_field.functionspace())) {
                    Log::info() << "storing distributed remapped field '" << tgt_name << ".msh'." << std::endl;
                    gmsh.write(functionspace::CellColumns(tgt_field.functionspace()).mesh());
                    tgt_field.haloExchange();
                    gmsh.write(tgt_field);
                }
                else {
                    Log::info() << "Target function space is NOT NodeColumns, hence no remapped test field will be stored." << std::endl;
                }
            }
        }
    }
    else {
        FunctionSpace src_fs;
        FunctionSpace tgt_fs;
        Interpolation interpolator;
        ATLAS_TRACE_SCOPE("Setup") {
        ATLAS_TRACE_SCOPE("Setup source and target") {
        timers.functionspace_setup.start();
        if (not args.has("force") && not available_combination(sgrid, tgrid, interpolation_method)) {
            Log::info() << "This interpolator is not available for these grids." << std::endl;
            Log::info() << "Run\n\t./bin/atlas-interpolations --sgrid " << sgrid.name() << " --tgrid " << tgrid.name()
                << "--list\n to see the list of available interpolators for these grid." << std::endl;
            return failed();
        };
        grid::Distribution dist_in, dist_out;
        create_fspaces(config, sgrid, tgrid, src_fs, tgt_fs, dist_in, dist_out);
        timers.functionspace_setup.stop();
        }

        ATLAS_TRACE_SCOPE("Setup interpolator") {
        timers.interpolation_setup.start();
        if (interpolation_method == "grid-box-average") {
            interpolator = Interpolation{ config, sgrid, tgrid };
        }
        else {
            interpolator = Interpolation{ config, src_fs, tgt_fs };
        }
        timers.interpolation_setup.stop();
        }
        }

        double timer_functionspace_setup = timers.functionspace_setup.elapsed() / mpi::comm().size();
        double timer_interpolation_setup = timers.interpolation_setup.elapsed() / mpi::comm().size();
        mpi::comm().allReduceInPlace(&timer_functionspace_setup, 1, eckit::mpi::sum());
        mpi::comm().allReduceInPlace(&timer_interpolation_setup, 1, eckit::mpi::sum());
        Log::info() << "Grid + FunctionSpace timer\t: " << timer_functionspace_setup * 1000. << " [ms]"  << std::endl;
        Log::info() << "Interpolation setup timer\t: " << timer_interpolation_setup * 1000. << " [ms]"  << std::endl;

        if (args.has("partest")) {
            ATLAS_TRACE_SCOPE("partest") {
                auto tgt_field = interpolator.target().createField<double>();
                auto field_in = interpolator.source().createField<double>();
                auto lonlat_in = array::make_view<double,2>(interpolator.source().lonlat());
                auto view_in = array::make_view<double,1>(field_in);
                for(idx_t j = 0; j < field_in.size(); ++j) {
                    view_in(j) = util::function::vortex_rollup(lonlat_in(j,0), lonlat_in(j,1), 1.);
                }
                field_in.set_dirty();
                field_in.haloExchange();
                ATLAS_TRACE_SCOPE("Interpolation") {
                    field_in.set_dirty();
                    field_in.haloExchange();
                    interpolator.execute(field_in, tgt_field);
                }
                std::string tgt_name = "par-tfield_" + get_matrix_name(sgrid_name, tgrid_name, interpolation_method);
                output::Gmsh gmsh(tgt_name + ".msh", Config("coordinates", "lonlat") | Config("ghost", "true"));
                if( functionspace::NodeColumns(tgt_field.functionspace())) {
                    Log::info() << "storing distributed remapped field '" << tgt_name << ".msh'." << std::endl;
                    gmsh.write(functionspace::NodeColumns(tgt_field.functionspace()).mesh());
                    tgt_field.haloExchange();
                    gmsh.write(tgt_field);
                }
                else {
                    Log::info() << "Target function space is NOT NodeColumns, hence no remapped test field will be stored." << std::endl;
                }
            }
        }

        auto matrix = interpolation::assemble_global_matrix(interpolator);

        if (mpi::comm().rank() == 0) {

            ATLAS_TRACE_SCOPE("store matrix") {
                std::string matrix_name;
                args.get("matrix", matrix_name);
                if (matrix_name == "") {
                    matrix_name = get_matrix_name(sgrid_name, tgrid_name, interpolation_method, matrix_format);
                }
                write_matrix(matrix, matrix_name, matrix_format);
            }

            if (args.has("test")) {
                ATLAS_TRACE("test");
                // Allocate and initialise own memory here to show possibilities
                std::vector<double> src_data(sgrid.size());
                std::vector<double> tgt_data(tgrid.size());

                ATLAS_TRACE_SCOPE("initialize source") {
                    idx_t n{0};
                    for (auto p : sgrid.lonlat()) {
                        src_data[n++] = util::function::vortex_rollup(p.lon(), p.lat(), 1.);
                    }
                }
                auto src = eckit::linalg::Vector(src_data.data(), src_data.size());
                auto tgt = eckit::linalg::Vector(tgt_data.data(), tgt_data.size());
                auto eckit_matrix = atlas::linalg::make_non_owning_eckit_sparse_matrix(matrix);
                timers.global_matrix_exe.start();
                eckit::linalg::LinearAlgebraSparse::backend().spmv(eckit_matrix, src, tgt);
                timers.global_matrix_exe.stop();
                Log::info() << "Global matrix-multiply timer  \t: " << 1000.*timers.global_matrix_exe.elapsed()  << " [ms]" << std::endl;
                Log::info() << "Global matrix non-zero entries\t: " << matrix.nnz() << std::endl;
                Log::info() << "Global matrix footprint       \t: " << eckit_matrix.footprint() << " B" << std::endl;

                ATLAS_TRACE_SCOPE("output from proc 0") {
                    mpi::Scope mpi_scope("self"); 
                    std::string tgt_name = "tfield_cache_" + get_matrix_name(sgrid_name, tgrid_name, interpolation_method);
                    auto tgt_field = Field(tgt_name, tgt_data.data(), array::make_shape(tgt_data.size())); // wrap
                    output::Gmsh gmsh(tgt_name + ".msh", Config("coordinates", "lonlat") | Config("ghost", "true"));
                    gmsh.write(MeshGenerator("structured", util::Config("three_dimensional",true)).generate(tgrid));
                    Log::info() << "storing remapped field '" << tgt_field.name() << ".msh'." << std::endl;
                    gmsh.write(tgt_field, StructuredColumns(tgrid));
                }
            }
        }
    }

    return success();
}


//-----------------------------------------------------------------------------


void AtlasInterpolations::interpolation_config(Config& config) {
    std::string interpolation_method;
    config.get("interpolation_method", interpolation_method);
    config.set("type", interpolation_method);
    config.set("halo", 1);
    if (interpolation_method.find("cubic") != std::string::npos) {
        config.set("halo", 2);
    }
    if (interpolation_method == "k-nearest-neighbours") {
        config.set("k-nearest-neighbours", 4);
        config.set("halo", 2);
    }
    if (interpolation_method == "conservative-spherical-polygon") {
        config.set("order", 1);
        config.set("halo", 1);
    }
    if (interpolation_method == "conservative-spherical-polygon-2") {
        config.set("type", "conservative-spherical-polygon");
        config.set("order", 2);
        config.set("halo", 2);
    }
    config.set("name", interpolation_method);
}


std::string AtlasInterpolations::get_matrix_name(std::string& sgrid, std::string& tgrid,
    std::string& interpolation_name, std::string matrix_format) {
    if (matrix_format == "eckit") {
        return "remap_" + sgrid + "_" + tgrid + "_" + interpolation_name + ".eckit";
    }
    if (matrix_format == "scrip") {
        return "remap_" + sgrid + "_" + tgrid + "_" + interpolation_name + ".nc";
    }
    if (matrix_format.empty()) {
        return "remap_" + sgrid + "_" + tgrid + "_" + interpolation_name;
    }
    ATLAS_NOTIMPLEMENTED;
}


void AtlasInterpolations::create_fspaces(Config& config, const Grid& input_grid, const Grid& output_grid,
        FunctionSpace& fs_in, FunctionSpace& fs_out, grid::Distribution& dist_in, grid::Distribution& dist_out) {
    interpolation_config(config);
    auto interpolation_str = config.getString("type");
    if (interpolation_str == "finite-element" || interpolation_str == "unstructured-bilinear-lonlat") {
        Mesh inmesh;
        ATLAS_TRACE_SCOPE("source functionspace") {
            if (input_grid.type() == "ORCA") {
                if (mpi::size() > 1 ) {
                    ATLAS_NOTIMPLEMENTED; // Cannot use orca in parallel as a source!
                }
                dist_in = grid::Distribution(input_grid, grid::Partitioner{"serial"});
                inmesh = Mesh(input_grid, dist_in);
                fs_in = functionspace::NodeColumns(inmesh);
            }
            else {
                dist_in = grid::Distribution(input_grid, input_grid.partitioner());
                inmesh = Mesh(input_grid, dist_in);
                fs_in = functionspace::NodeColumns(inmesh);
                // fs_in = functionspace::NodeColumns(inmesh, option::halo(1));
            }
        }
        ATLAS_TRACE_SCOPE("target functionspace") {
            auto partitioner = mpi::size() == 1 ? grid::Partitioner("serial") : grid::MatchingPartitioner(inmesh);
            dist_out = grid::Distribution(output_grid, partitioner);
            if( output_grid.type() == "ORCA") {
                // fs_out = functionspace::NodeColumns(Mesh(output_grid, dist_out));
                fs_out = functionspace::PointCloud(output_grid, dist_out);
            }
            else {
                fs_out = functionspace::PointCloud(output_grid, dist_out);
            }
        }
    }
    else if (interpolation_str == "conservative-spherical-polygon") {
        bool src_cell_data = config.getBool("src_cell_data");
        bool tgt_cell_data = config.getBool("tgt_cell_data");
        auto tgt_mesh_config = output_grid.meshgenerator() | option::halo(0);
        int halo_size;
        config.get("halo", halo_size);
        Mesh tgt_mesh;
        ATLAS_TRACE_SCOPE("target functionspace") {
            bool orca = (output_grid.type() == "ORCA");
            if (orca) {
                if (mpi::size() > 1 ) {
                    ATLAS_NOTIMPLEMENTED; // the source must match target but cannot be for orca grids
                }
                dist_out = grid::Distribution(output_grid, grid::Partitioner{"serial"});
                tgt_mesh = MeshGenerator(tgt_mesh_config).generate(output_grid);
            }
            else {
                dist_out = grid::Distribution(output_grid, grid::Partitioner{output_grid.partitioner()});
                tgt_mesh = MeshGenerator(tgt_mesh_config).generate(output_grid, dist_out);
            }
            if (tgt_cell_data) {
                fs_out = functionspace::CellColumns(tgt_mesh, option::halo(0));
            }
            else {
                fs_out = functionspace::NodeColumns(tgt_mesh, option::halo(orca ? 0 : 1));
            }
        }
        ATLAS_TRACE_SCOPE("source functionspace") {
            bool orca = (input_grid.type() == "ORCA");
            halo_size = (orca ? 0 : halo_size);
            Config src_mesh_config;
            src_mesh_config = input_grid.meshgenerator() | option::halo(halo_size);
            Mesh src_mesh;
            auto partitioner = mpi::size() == 1 ? grid::Partitioner("serial") : grid::MatchingPartitioner(tgt_mesh);
            dist_in = grid::Distribution(input_grid, partitioner);
            src_mesh = MeshGenerator(src_mesh_config).generate(input_grid, dist_in);
            if (src_cell_data) {
                fs_in = functionspace::CellColumns(src_mesh, option::halo(halo_size));
            }
            else {
                fs_in = functionspace::NodeColumns(src_mesh, option::halo(halo_size + (orca ? 0 : 1)));
            }
        }
    }
    else if (interpolation_str == "nearest-neighbour" || interpolation_str == "k-nearest-neighbours"
        || interpolation_str == "grid-box-average") {
        Mesh inmesh;
        ATLAS_TRACE_SCOPE("source functionspace") {
            dist_in = grid::Distribution(input_grid, grid::Partitioner{input_grid.partitioner()});
            fs_in = PointCloud(input_grid, dist_in);
            inmesh = Mesh(input_grid);
        }
        ATLAS_TRACE_SCOPE("target functionspace") {
            auto partitioner = mpi::size() == 1 ? grid::Partitioner("serial") : grid::MatchingPartitioner(inmesh);
            dist_out = grid::Distribution(output_grid, partitioner);
            fs_out = functionspace::PointCloud(output_grid, dist_out);
        }
    }
    else {
        ATLAS_TRACE_SCOPE("source functionspace") {
            dist_in = grid::Distribution(input_grid, grid::Partitioner{input_grid.partitioner()});
            fs_in = functionspace::StructuredColumns(input_grid, dist_in, config);
        }
        ATLAS_TRACE_SCOPE("target functionspace") {
            auto partitioner = mpi::size() == 1 ? grid::Partitioner("serial") : grid::MatchingPartitioner(fs_in);
            dist_out = grid::Distribution(output_grid, partitioner);
            fs_out = functionspace::PointCloud(output_grid, dist_out, config);
        }
    }
}


Matrix AtlasInterpolations::read_matrix(std::string matrix_name, std::string format) {
    if (format == "eckit") {
        eckit::linalg::SparseMatrix eckit_matrix;
        Log::info() << "reading matrix '" << matrix_name << std::endl;
        eckit_matrix.load(matrix_name + ".eckit");
        return atlas::linalg::make_sparse_matrix_storage(std::move(eckit_matrix));
    }
    else if (format == "scrip") {
        Log::info() << "reading matrix '" << matrix_name << std::endl;
        return ScripIO::read(matrix_name + ".nc");
    }
    else {
        ATLAS_NOTIMPLEMENTED;
    }
    return Matrix{};
}


void AtlasInterpolations::write_matrix(const Matrix& matrix, std::string matrix_name, std::string format) {
    Log::info() << "storing matrix '" << matrix_name << std::endl;
    if (format == "eckit") {
        auto eckit_matrix = linalg::make_non_owning_eckit_sparse_matrix(matrix);
        eckit_matrix.save(matrix_name);
    }
    else if (format == "scrip") {
        ScripIO::write(matrix, matrix_name);
    }
    else {
        ATLAS_NOTIMPLEMENTED;
    }
}


bool AtlasInterpolations::grid_type_same(const char* c1, const char* c2) {
    auto gridtype = [](const char* c) {
        if ((c[0] == 'e') || (c[0] == 'O' && c[1] == 'R')) {
            return std::string("ORCA");
        }
        if (c[0] == 'H') {
            return std::string("healpix");
        }
        if ((c[0] == 'f') || (c[0] == 'D') || (c[0] == 't') || ((c[0] == 'C') && (c[1] == 'O')) || ((c[0] == 'N') && (c[1] == 'G'))) {
            return std::string("fesom");
        }
        if (c[0] == 'N' || c[0] == 'O' || c[0] == 'F' || c[0] == 'L' || c[0] == 'S') {
            return std::string("structured");
        }
        if ((c[0] == 'C') && (c[1] == 'S')) {
            return std::string("cubedsphere");
        }
        return std::string("unsupported");
    };
    auto gt1 = gridtype(c1);
    return (gt1 == gridtype(c2) && gt1 != "unsupported");
};


bool AtlasInterpolations::init_interpolation_database() {
    std::ifstream ifile("db_atlas_interpolations");
    if (!ifile) {
        std::cerr << "ERROR Cannot open ./db_atlas_interpolations.\n";
        return false;
    }
    // skip comment lines
    while (ifile.peek() == '#') {
        ifile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    }
    InterpolationStatus istat;
    while (ifile >> istat) {
        idatabase_.push_back(istat);
    }
    ifile.close();
    ATLAS_ASSERT(idatabase_.size() > 0);
    return true;
}


bool AtlasInterpolations::get_available_interpolations(Grid sgrid, Grid tgrid, std::vector<std::string>& serial_interpolations,
    std::vector<std::string>& parallel_interpolations) {
    if (idatabase_.size() == 0) {
        if (not init_interpolation_database()) {
            return false;
        };
    }
    auto filter_interpolations = [this, &sgrid, &tgrid](const char mpi_tasks, std::vector<std::string>& interp) { 
        for (int i = 0; i < idatabase_.size(); ++i) {
            if (std::find(interp.begin(), interp.end(), idatabase_[i].interpolation) != interp.end()) {
                continue;
            }
            if (idatabase_[i].status == "[OK]" && idatabase_[i].mpi_tasks[0] == mpi_tasks) {
                if (grid_type_same(idatabase_[i].sgrid.c_str(), sgrid.name().c_str())) {
                    for (int j = 0; j < idatabase_.size(); ++j) {
                        if (idatabase_[j].status == "[OK]" && idatabase_[i].interpolation == idatabase_[j].interpolation && idatabase_[j].mpi_tasks[0] == mpi_tasks) {
                            if (grid_type_same(idatabase_[j].tgrid.c_str(), tgrid.name().c_str())) {
                                interp.emplace_back(idatabase_[i].interpolation);
                                break;
                            }
                        }
                    }
                }
            }
        }
    };
    filter_interpolations('1', serial_interpolations);
    filter_interpolations('4', parallel_interpolations);
    std::sort(serial_interpolations.begin(), serial_interpolations.end());
    std::sort(parallel_interpolations.begin(), parallel_interpolations.end());
    return true;
}


bool AtlasInterpolations::available_combination(Grid sgrid, Grid tgrid, std::string interpolation_name) {
    if (idatabase_.size() == 0) {
        if (not init_interpolation_database()) {
            return false;
        };
    }
    std::vector<std::string> serial_interpolations;
    std::vector<std::string> parallel_interpolations;
    if (not get_available_interpolations(sgrid, tgrid, serial_interpolations, parallel_interpolations)) {
        return false;
    }

    if (mpi::comm().size() == 1) {
        return (std::find(serial_interpolations.begin(), serial_interpolations.end(), interpolation_name) != serial_interpolations.end());
    }
    else {
        return (std::find(parallel_interpolations.begin(), parallel_interpolations.end(), interpolation_name) != parallel_interpolations.end());
    }
}

} // namespace


int main(int argc, char** argv) {
    atlas::AtlasInterpolations tool(argc, argv);
    return tool.start();
}