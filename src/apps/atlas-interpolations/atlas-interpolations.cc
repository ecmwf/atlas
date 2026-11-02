/*
 * (C) Copyright 2025 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <unordered_set>

#include "eckit/filesystem/PathName.h"

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
#include "atlas/io/atlas-io.h"
#include "atlas/linalg/sparse.h"
#include "atlas/linalg/sparse/SparseMatrixTriplet.h"

#include "atlas/mesh/Mesh.h"
#include "atlas/meshgenerator.h"
#include "atlas/output/Gmsh.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/runtime/AtlasTool.h"
#include "atlas/util/Checksum.h"
#include "atlas/util/CoordinateEnums.h"
#include "atlas/util/KDTree.h"
#include "atlas/util/function/VortexRollup.h"


#include "atlas/util/function/MDPI_functions.h"
#include "atlas/util/function/SolidBodyRotation.h"
#include "atlas/util/function/SphericalHarmonic.h"
#include "atlas/util/function/VortexRollup.h"

#include "MaskIO.h"
#include "ScripIO.h"

#if ATLAS_HAVE_GRIB
#include "Grib.h"
#include "GribFileReader.h"
#endif

using atlas::functionspace::CellColumns;
using atlas::functionspace::NodeColumns;
using atlas::functionspace::PointCloud;
using atlas::functionspace::StructuredColumns;
using atlas::util::Config;

using Matrix    = atlas::linalg::SparseMatrixStorage;
using StopWatch = atlas::runtime::trace::StopWatch;

namespace atlas {

std::vector<double> tgt_ref;


// This app can:
//
// 1) list all available interpolators
//      atlas-interpolations --list
// 2) compute and store the global interpolation matrix to disk:
//      atlas-interpolations --s.grid=<sgrid> --t.grid=<tgrid> --i.type=<interp> --output-matrix --matrix.format=<fmt>
// 3) read in a stored global interpolation matrix from the disk:
//      atlas-interpolations --s.grid=<sgrid> --t.grid=<tgrid> --i.type=<interp> --read-matrix --matrix.format=<fmt>
// 4) extract the interpolation stencil for a few target points into a lightweight ASCII file for
//    visualisation/study (see plot_stencil.py), without writing the entire interpolated field:
//      atlas-interpolations --s.grid=<sgrid> --t.grid=<tgrid> --i.type=<interp> --output-stencil \
//                           --stencil.target-lonlat=<lon,lat[;lon,lat...]> [--s.mask=<mask> --t.mask=<mask>]
//
// Formats <fmt> can be 'eckit' (binary) or 'scrip' (netcdf).
//
// This app hides configurations of meshes, function spaces, halo size, etc. from the user.
//

class AtlasInterpolations : public AtlasTool {
    int execute(const AtlasTool::Args& args) override;
    std::string briefDescription() override {
        return "Atlas intepolation app to perform MPI-parallel interpolations for given source "
               "and target grids and/or generate remapping matrices.\n";
    }
    std::string usage() override { return name() + " [OPTION] ... [--help]"; }
    std::string longDescription() override {
        return "The 's.grid' and 't.grid' arguments are source and target grids name, "
               "for instance: N80, F40, O24, L64x33, CS-ED-12, etc. See './bin/atlas-grids --list' "
               "for the complete list of available grids. See './bin/atlas-interpolation --list' "
               "for a list of named interpolations. See './bin/atlas-interpolation --s.grid <grid1> "
               "--t.grid <grid2> --list' for the list of available interpolations for that grid pair.\n"
               "\n";
    }

    int numberOfPositionalArguments() override { return -1; }
    int minimumPositionalArguments() override { return 0; }
    Matrix read_matrix(std::string matrix_name, std::string matrix_format);
    void write_matrix(const Matrix& mat, std::string matrix_name, std::string matrix_format);
    void print_matrix(const Matrix&, std::ostream&);

    struct Timers {
        StopWatch functionspace_setup;
        StopWatch interpolation_setup;
        StopWatch interpolation_exe;
        StopWatch global_matrix_read;
        StopWatch global_matrix_setup;
        StopWatch global_matrix_exe;
    } timers;

public:
    AtlasInterpolations(int argc, char* argv[]): AtlasTool(argc, argv) {
        // Commands
        add_option(new eckit::option::Separator("Commands"));
        add_option(new SimpleOption<bool>("list", "List available interpolation types"));
        add_option(new SimpleOption<bool>("interpolate", "Test an interpolation type"));
        add_option(new SimpleOption<bool>(
            "output-stencil", "Extract interpolation stencil(s) for a few target points into a lightweight ASCII "
                              "file for visualisation/study (source+target points, weights, masks, and source mesh "
                              "connectivity when available), without writing the full interpolated field. "
                              "Configure with the 'stencil.*' options. Best used in serial (1 MPI task)."));

        // Configuration
        add_option(new eckit::option::Separator("Interpolation options"));
        add_option(new SimpleOption<std::string>(
            "i.knn", "number of nearest neighbours in case i.type=knn or k-nearest-neighbours. (default=1)"));
        add_option(new SimpleOption<bool>("i.normalise", "Normalise weights"));
        add_option(new SimpleOption<std::string>("s.grid", "source grid"));
        add_option(new SimpleOption<std::string>("s.grib", "file path to input"));
        add_option(new SimpleOption<std::string>("s.mask", "source mask"));
        add_option(new SimpleOption<std::string>("t.grid", "target grid"));
        add_option(new SimpleOption<std::string>("t.mask", "target mask"));
        add_option(new SimpleOption<std::string>("i.type", "interpolation type"));

        add_option(new eckit::option::Separator("Advanced configuration"));
        add_option(new SimpleOption<std::string>("s.fs", "source function space [structured, nodes, cells, points]"));
        add_option(new SimpleOption<std::string>("t.fs", "target function space [structured, nodes, cells, points]"));
        add_option(new SimpleOption<long>("s.halo", "override default source halo"));
        add_option(new SimpleOption<long>("t.halo", "override default target halo"));
        add_option(new SimpleOption<std::string>("p.type", "partitioner type"));

        add_option(new eckit::option::Separator("Advanced commands"));
        add_option(new SimpleOption<bool>("read-matrix", "Read matrix to speed up interpolator"));
        add_option(new SimpleOption<bool>("test-matrix", "Test matrix in serial"));
        add_option(new SimpleOption<bool>("test-interpolator", "Test interpolator in parallel"));

        add_option(new eckit::option::Separator("Output options"));
        add_option(new SimpleOption<bool>("output-matrix", "Write interpolation matrix"));
        add_option(new SimpleOption<std::string>(
            "matrix.name", "Name of the remapping matrix. If not provided a default unique name will be chosen"));
        add_option(new SimpleOption<std::string>("matrix.format", "Format of the remapping matrix: eckit, scrip"));
        add_option(new SimpleOption<bool>("output-gmsh", "Set to enable gmsh output for tests"));
        add_option(new SimpleOption<bool>("output-gmsh-source", "Include source field in Gmsh output"));
        add_option(
            new SimpleOption<std::string>("gmsh.coordinates", "Choose the coordinates in gmsh output: {lonlat, xyz}"));
        add_option(new SimpleOption<std::string>("gmsh.missing_value.policy",
                                                 "Missing-value policy in Gmsh output: {skip, fill, nan, preserve}"));
        add_option(new SimpleOption<double>("gmsh.missing_value.fill",
                                            "Fill value used by the Gmsh missing-value policy 'fill'"));
        add_option(new SimpleOption<std::string>("gmsh.masked_value.policy",
                                                 "Masked-value policy in Gmsh output: {skip, fill, nan, preserve}"));
        add_option(new SimpleOption<double>("gmsh.masked_value.fill",
                                            "Fill value used by the Gmsh masked-value policy 'fill'"));

        add_option(new SimpleOption<bool>("output-checksum", "Compute checksums"));
        add_option(new SimpleOption<std::string>("checksum.file", "Path of files for checksums [default='checksum']"));

        // Stencil extraction options
        add_option(new eckit::option::Separator("Stencil extraction options (used with --stencil)"));
        add_option(new SimpleOption<std::string>(
            "stencil.target-lonlat", "Target location(s) to extract; the nearest target grid point is used. Format: "
                                     "'lon,lat' or several separated by ';', e.g. --stencil.target-lonlat=10,45;-30,60"));
        add_option(new SimpleOption<std::string>(
            "stencil.target-index",
            "Target grid point(s) to extract, given as 1-based global index/indices, comma-separated, "
            "e.g. --stencil.target-index=1024,2048"));
        add_option(new SimpleOption<long>(
            "stencil.source-context", "Number of nearest source grid points added around each stencil point to give "
                                      "context (e.g. to see the source mask around the stencil). (default=8)"));
        add_option(new SimpleOption<long>(
            "stencil.target-context", "Number of nearest target grid points added around the target point to give "
                                      "context. (default: same as --stencil.source-context)"));
        add_option(new SimpleOption<std::string>("stencil.file", "Output ASCII file for stencils [default='stencil.txt']"));

        // Initial condition options
        add_option(new eckit::option::Separator("Initial condition options"));
        add_option(new SimpleOption<std::string>("init",
                                                 "Setup initial source field [ constant, spherical_harmonic, "
                                                 "vortex_rollup (default), solid_body_rotation_wind_magnitude, "
                                                 "atlas-field-magnitude (EXPERIMENTAL) ]"));
        add_option(new SimpleOption<std::string>(
            "init-file", "AtlasIO file used by experimental init=atlas-field-magnitude"));
        add_option(
            new SimpleOption<double>("solid_body_rotation.angle", "Angle of solid body rotation (default = 0.)"));
        add_option(new SimpleOption<double>("vortex_rollup.t", "Value that controls vortex rollup (default = 0.5)"));
        add_option(new SimpleOption<double>("constant.value", "Value that is assigned in case init==constant)"));
        add_option(new SimpleOption<long>("spherical_harmonic.n", "total wave number 'n' of a spherical harmonic"));
        add_option(new SimpleOption<long>("spherical_harmonic.m", "zonal wave number 'm' of a spherical harmonic"));

        add_option(new SimpleOption<double>("missing_value", "Value used with --test-matrix to mark missing values"));
    }
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

        util::function::SphericalHarmonic Y(n, m);
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

class SourceInitialCondition {
public:
    SourceInitialCondition(const AtlasTool::Args& args, const Grid& grid) {
        const std::string init = args.getString("init", "vortex_rollup");
        if (init != "atlas-field-magnitude") {
            function_ = get_init(args);
            return;
        }

        Log::warning() << "\n"
                       << "**************************************************************************\n"
                       << "* WARNING: --init=atlas-field-magnitude IS EXPERIMENTAL                  *\n"
                       << "* Its interface and initialisation semantics may change.                 *\n"
                       << "* It expects two fields in the input file:                               *\n"
                       << "*   the first and second components of the vector field.                 *\n"
                       << "* Verify the input fields and missing-value handling before using it.    *\n"
                       << "**************************************************************************"
                       << std::endl;

        if (!args.has("init-file")) {
            ATLAS_THROW_EXCEPTION("--init=atlas-field-magnitude requires --init-file");
        }
        const std::string path = args.getString("init-file");
        std::string record_grid_name;
        std::vector<double> first;
        std::vector<double> second;
        io::RecordReader reader(path);
        reader.read("grid.name", record_grid_name);
        reader.read("fields[0].array", first);
        reader.read("fields[1].array", second);
        reader.wait();

        if (record_grid_name != grid.name()) {
            ATLAS_THROW_EXCEPTION("Initial-condition file '" << path << "' is for grid '" << record_grid_name
                                                               << "', expected '" << grid.name() << "'");
        }
        if (first.size() != grid.size() || second.size() != grid.size()) {
            ATLAS_THROW_EXCEPTION("Initial-condition file '" << path << "' contains field sizes " << first.size()
                                                               << " and " << second.size() << ", expected "
                                                               << grid.size());
        }

        values_ = std::move(first);
        const double missing_value = args.getDouble("missing_value", std::numeric_limits<double>::lowest());
        for (std::size_t i = 0; i < values_.size(); ++i) {
            if (values_[i] == missing_value || second[i] == missing_value) {
                values_[i] = missing_value;
            }
            else {
                values_[i] = std::hypot(values_[i], second[i]);
            }
        }
    }

    double operator()(gidx_t index, const PointLonLat& point) const {
        return values_.empty() ? function_(point) : values_.at(static_cast<std::size_t>(index));
    }

private:
    std::function<double(const PointLonLat&)> function_;
    std::vector<double> values_;
};

double elapsed_ms(const StopWatch& timer, bool single_rank = false) {
    double timer_elapsed = timer.elapsed() * 1000.;
    if (single_rank) {
        return timer_elapsed;
    }
    mpi::comm().allReduceInPlace(&timer_elapsed, 1, eckit::mpi::sum());
    return timer_elapsed / mpi::size();
}

std::string get_interpolation_method(const AtlasTool::Args& args) {
    std::string method   = args.getString("i.type", "finite-element");
    auto method_contains = [method](std::string_view substr) { return method.find(substr) != std::string_view::npos; };
    if (method_contains("cons")) {
        return method_contains("2") ? "cons2" : "cons";
    }
    return method;
}

std::string get_extension(const std::string& str) {
    std::size_t i = str.rfind('.', str.length());
    if (i != std::string::npos) {
        return (str.substr(i + 1, str.length() - i));
    }
    return ("");
}

std::string get_basename(const std::string& str) {
    std::string basename = str;
    std::string ext      = get_extension(str);
    if (ext == "eckit" || ext == "nc") {
        std::size_t lastindex = str.find_last_of(".");
        basename              = str.substr(0, lastindex);
    }
    return basename;
}

bool string_contains(std::string_view str, std::string_view substr) {
    return str.find(substr) != std::string_view::npos;
};

bool string_starts_with(std::string_view str, std::string_view substr) {
    return str.rfind(substr, 0) == 0;
};

Grid get_sgrid(const AtlasTool::Args& args) {
    ATLAS_ASSERT(args.has("s.grib") or args.has("s.grid"),
                 "Could not detect grid because s.grib or s.grid is not defined");
    if (args.has("s.grid")) {
        return Grid(args.getString("s.grid"));
    }
    else if (args.has("s.grib")) {
#if ATLAS_HAVE_GRIB
        // Read from GRIB
        GribFileReader grib_reader(args.getString("s.grib"));
        ATLAS_ASSERT(grib_reader.count());
        return grib_reader.grib().get_grid();
#else
        throw_Exception("GRIB support is not available", Here());
#endif
    }
    throw_Exception("Could not detect grid because s.grid or s.grib is not defined", Here());
}

Grid get_tgrid(const AtlasTool::Args& args) {
    return Grid(args.getString("t.grid", "O32"));
}

bool gridpoints_are_cells(const Grid& g) {
    if (g.type() == "healpix") {
        return true;
    }
    return false;
}

std::string get_matrix_format(const AtlasTool::Args& args) {
    if (args.has("matrix.format")) {
        return args.getString("matrix.format");
    }
    if (args.has("matrix.name")) {
        auto ext = get_extension(args.getString("matrix.name"));
        if (ext == "eckit") {
            return "eckit";
        }
        if (ext == "nc") {
            return "scrip";
        }
    }
    return "eckit";
}

std::string get_matrix_name(const AtlasTool::Args& args) {
    if (args.has("matrix.name")) {
        auto matrix = args.getString("matrix.name");
        return get_basename(matrix);
    }
    auto sgrid              = get_sgrid(args).name();
    auto tgrid              = args.getString("t.grid");
    auto interpolation_name = get_interpolation_method(args);
    return "remap_" + sgrid + "_" + tgrid + "_" + interpolation_name;
}

std::string get_mask_format(const std::string& mask) {
    auto ext = get_extension(mask);
    if (ext == "nc") {
        return "scrip";
    }
    else if (ext == "atlas") {
        return "atlas";
    }
    ATLAS_NOTIMPLEMENTED;
}

void read_mask(const std::string& mask_name, array::mdspan<int, array::dims<1>> mask) {
    if (get_mask_format(mask_name) == "scrip") {
        ScripIO::read_mask(mask_name, mask);
    }
    else if (get_mask_format(mask_name) == "atlas") {
        MaskIO::read(mask_name, mask);
    }
    else {
        ATLAS_NOTIMPLEMENTED;
    }
}

void set_grid_mask(const std::string& mask_name, FunctionSpace& functionspace) {
    Field global_indices = functionspace.createField(functionspace.global_index(), option::global());
    array::make_view<gidx_t, 1>(global_indices).assign(0);
    functionspace.gather(functionspace.global_index(), global_indices);

    functionspace.setMask([&](int* mask, size_t size) {
        std::vector<int> grid_mask(functionspace.grid().size());
        read_mask(mask_name, array::mdspan<int, array::dims<1>>(grid_mask.data(), grid_mask.size()));

        auto global_index = array::make_view<gidx_t, 1>(global_indices);
        ATLAS_ASSERT(global_index.size() == size);
        for (idx_t i = 0; i < global_index.size(); ++i) {
            if (global_index(i) == 0) {
                mask[i] = 0;
                continue;
            }
            const gidx_t grid_index = global_index(i) - 1;
            ATLAS_ASSERT(grid_index >= 0 && grid_index < grid_mask.size());
            mask[i] = grid_mask[grid_index];
        }
    });
}

Config get_interpolation_config(const Grid& sgrid, const Grid& tgrid, const AtlasTool::Args& args) {
    Config c;
    std::string type = get_interpolation_method(args);
    if (type == "cons") {
        c.set("type", "conservative-spherical-polygon");
        c.set("order", 1);
        c.set("src_cell_data", gridpoints_are_cells(sgrid));
        c.set("tgt_cell_data", gridpoints_are_cells(tgrid));
        bool normalise = true;
        args.get("i.normalise", normalise);
        c.set("normalise", normalise);
    }
    else if (type == "cons2") {
        c.set("type", "conservative-spherical-polygon");
        c.set("order", 2);
        c.set("src_cell_data", gridpoints_are_cells(sgrid));
        c.set("tgt_cell_data", gridpoints_are_cells(tgrid));
        bool normalise = true;
        args.get("i.normalise", normalise);
        c.set("normalise", normalise);
    }
    else {
        c.set("type", type);
    }
    if (type == "k-nearest-neighbours" || type == "knn") {
        c.set("type", "k-nearest-neighbours");
        c.set("k-nearest-neighbours", 1);
        if (args.has("i.knn")) {
            c.set("k-nearest-neighbours", args.getLong("i.knn"));
        }
        if (args.has("i.k-nearest-neighbours")) {
            c.set("k-nearest-neighbours", args.getLong("i.k-nearest-neighbours"));
        }
    }
    return c;
}

Mesh create_mesh_with_gridpoints_as_nodes(const Grid& grid, const grid::Distribution& dist, int halo) {
    if (gridpoints_are_cells(grid)) {
        if (StructuredGrid(grid)) {
            return MeshGenerator("structured", option::halo(halo)).generate(grid, dist);
        }
        else {
            ATLAS_THROW_EXCEPTION("Could not detect the mesh generator for grid \""
                                  << grid.name() << "\" that treats grid points as nodes");
        }
    }
    return Mesh(grid, dist, option::halo(halo));
}

std::pair<FunctionSpace, FunctionSpace> get_fs(const Grid& sgrid, const Grid& tgrid, const AtlasTool::Args& args) {
    enum class fs_type
    {
        unset,
        structured,
        cells,
        nodes,
        points
    };
    auto to_fs_type = [](const std::string& str) -> fs_type {
        if (str == "structured") {
            return fs_type::structured;
        }
        if (str == "cells") {
            return fs_type::cells;
        }
        if (str == "nodes") {
            return fs_type::nodes;
        }
        if (str == "points") {
            return fs_type::points;
        }
        return fs_type::unset;
    };
    FunctionSpace sfs, tfs;
    std::string interpolation_method = get_interpolation_method(args);
    fs_type sfs_type                 = fs_type::unset;
    fs_type tfs_type                 = fs_type::unset;
    bool target_follows_source       = true;
    int shalo                        = 0;
    int thalo                        = 0;
    if (string_starts_with(interpolation_method, "cons")) {
        int order             = string_starts_with(interpolation_method, "cons2") ? 2 : 1;
        target_follows_source = false;
        sfs_type              = gridpoints_are_cells(sgrid) ? fs_type::cells : fs_type::nodes;
        tfs_type              = gridpoints_are_cells(tgrid) ? fs_type::cells : fs_type::nodes;
        if (tfs_type == fs_type::nodes) {
            thalo = 1;
        }
        if (mpi::size() == 1 && (order == 2 || sfs_type == fs_type::nodes)) {
            shalo = 1;
        }
        if (mpi::size() == 1 && tgrid.type() == "ORCA") {
            thalo = 0;
        }
        if (mpi::size() == 1 && sgrid.type() == "ORCA") {
            shalo = 0;
        }
    }
    else if (string_starts_with(interpolation_method, "structured-")) {
        if (not StructuredGrid(sgrid)) {
            ATLAS_THROW_EXCEPTION("Interpolation method " << interpolation_method
                                                          << " requires that the source grid is a StructuredGrid");
        }
        sfs_type = fs_type::structured;
        tfs_type = fs_type::points;
        if (args.has("t.fs")) {
            tfs_type = to_fs_type(args.getString("t.fs"));
        }
        shalo = 1;
        if (string_contains(interpolation_method, "cubic")) {
            shalo = 2;
        }
    }
    else if (string_starts_with(interpolation_method, "grid-box")) {
        if (not StructuredGrid(sgrid)) {
            ATLAS_THROW_EXCEPTION("Interpolation method " << interpolation_method
                                                          << " requires that the source grid is a StructuredGrid");
        }
        if (not StructuredGrid(tgrid)) {
            ATLAS_THROW_EXCEPTION("Interpolation method " << interpolation_method
                                                          << " requires that the target grid is a StructuredGrid");
        }
        sfs_type = fs_type::structured;
        tfs_type = fs_type::structured;
    }
    else if (string_contains(interpolation_method, "nearest")) {
        if (StructuredGrid(sgrid)) {
            sfs_type = fs_type::structured;
        }
        else {
            sfs_type = fs_type::points;
        }
        tfs_type = fs_type::points;
        if (args.has("t.fs")) {
            tfs_type = to_fs_type(args.getString("t.fs"));
        }
        if (args.has("s.fs")) {
            sfs_type = to_fs_type(args.getString("s.fs"));
        }
    }
    else {
        sfs_type = fs_type::nodes;
        tfs_type = fs_type::points;
        if (args.has("t.fs")) {
            tfs_type = to_fs_type(args.getString("t.fs"));
        }
    }
    if (target_follows_source && sfs_type == fs_type::points && mpi::size() > 1) {
        if (StructuredGrid(sgrid)) {
            sfs_type = fs_type::structured;
        }
        else {
            Log::warning()
                << "WARNING: The source will be meshed for reason of creating the distribution of the target grid"
                << std::endl;
            sfs_type = fs_type::nodes;
        }
    }
    args.get("s.halo", shalo);
    args.get("t.halo", thalo);

    if (target_follows_source && sgrid.type() == "ORCA" && mpi::size() > 1) {
        ATLAS_THROW_EXCEPTION("ORCA grid is not yet supported as source in parallel for interpolation method \""
                              << interpolation_method << "\"");
    }
    if (not target_follows_source && tgrid.type() == "ORCA" && mpi::size() > 1) {
        ATLAS_THROW_EXCEPTION("ORCA grid is not yet supported as target in parallel for interpolation method \""
                              << interpolation_method << "\"");
    }
    if (string_starts_with(interpolation_method, "grid-box") && shalo > 0) {
        ATLAS_THROW_EXCEPTION("Interpolation method \"" << interpolation_method
                                                        << "\" only supported with --s.halo=0 (default)");
    }

    auto create_grid_dist = [&](const Grid& grid) {
        grid::Distribution dist;
        if (mpi::size() == 1) {
            dist = grid::Distribution(grid, grid::Partitioner{"serial"});
        }
        else {
            if (args.has("p")) {
                dist = grid::Distribution(grid, grid::Partitioner{args.getString("p.type")});
            }
            else {
                dist = grid::Distribution(grid, grid.partitioner());
            }
        }
        return dist;
    };

    auto create_matching_grid_dist = [&](const Grid& grid, const FunctionSpace& fs_to_match) {
        grid::Distribution dist;
        if (mpi::size() == 1) {
            dist = grid::Distribution(grid, grid::Partitioner{"serial"});
        }
        else {
            if (auto fs = functionspace::NodeColumns(fs_to_match)) {
                dist = grid::Distribution{grid, grid::MatchingPartitioner{fs.mesh()}};
            }
            else if (auto fs = functionspace::CellColumns(fs_to_match)) {
                dist = grid::Distribution{grid, grid::MatchingPartitioner{fs.mesh()}};
            }
            else {
                dist = grid::Distribution{grid, grid::MatchingPartitioner{fs_to_match}};
            }
        }
        return dist;
    };

    auto create_fs = [&](fs_type type, const Grid& grid, const grid::Distribution& dist, int halo) {
        FunctionSpace fs;
        if (type == fs_type::cells) {
            auto mesh = Mesh(grid, dist, option::halo(halo));
            fs        = CellColumns(mesh, option::halo(halo));
        }
        else if (type == fs_type::nodes) {
            auto mesh = create_mesh_with_gridpoints_as_nodes(grid, dist, halo);
            fs        = NodeColumns(mesh, option::halo(halo));
        }
        else if (type == fs_type::structured) {
            fs = StructuredColumns(grid, dist, option::halo(halo));
        }
        else if (type == fs_type::points) {
            fs = PointCloud(grid, dist, option::halo(halo));
        }
        else {
            ATLAS_NOTIMPLEMENTED;
        }
        return fs;
    };

    if (target_follows_source) {
        // This is the usual case for most interpolation methods
        grid::Distribution sdist, tdist;
        ATLAS_TRACE_SCOPE("Partition source grid") {
            sdist = create_grid_dist(sgrid);
        }
        ATLAS_TRACE_SCOPE("Create source function space") {
            sfs = create_fs(sfs_type, sgrid, sdist, shalo);
        }
        ATLAS_TRACE_SCOPE("Partition target grid") {
            tdist = create_matching_grid_dist(tgrid, sfs);
        }
        ATLAS_TRACE_SCOPE("Create target function space") {
            tfs = create_fs(tfs_type, tgrid, tdist, thalo);
        }
    }
    else {
        // This is only needed for conservative interpolation
        ATLAS_ASSERT(string_starts_with(interpolation_method, "cons"));
        grid::Distribution sdist, tdist;
        ATLAS_TRACE_SCOPE("Partition target grid") {
            tdist = create_grid_dist(tgrid);
        }
        ATLAS_TRACE_SCOPE("Create target function space") {
            tfs = create_fs(tfs_type, tgrid, tdist, thalo);
        }
        ATLAS_TRACE_SCOPE("Partition source grid") {
            sdist = create_matching_grid_dist(sgrid, tfs);
        }
        ATLAS_TRACE_SCOPE("Create source function space") {
            sfs = create_fs(sfs_type, sgrid, sdist, shalo);
        }
    }
    return std::make_pair(sfs, tfs);
}


Config gmsh_field_config(const AtlasTool::Args& args, const util::Config& default_config = util::Config() ) {
    Config config(default_config);
    if (args.has("gmsh.missing_value.policy")) {
        config.set("missing_value.policy", args.getString("gmsh.missing_value.policy"));
    }
    if (args.has("gmsh.masked_value.policy")) {
        config.set("masked_value.policy", args.getString("gmsh.masked_value.policy"));
    }
    if (args.has("gmsh.missing_value.fill")) {
        config.set("missing_value.fill", args.getDouble("gmsh.missing_value.fill"));
    }
    if (args.has("gmsh.masked_value.fill")) {
        config.set("masked_value.fill", args.getDouble("gmsh.masked_value.fill"));
    }
    std::string default_policy = "skip";
    if (not config.has("missing_value.policy")) {
        config.set("missing_value.policy", default_policy);
    }
    if (not config.has("masked_value.policy")) {
        config.set("masked_value.policy", default_policy);
    }
    return config;
}

void set_missing_value_metadata(const AtlasTool::Args& args, Field& field) {
    double missing_value = std::numeric_limits<double>::lowest();
    args.get("missing_value", missing_value);
    field.metadata().set("missing_value", missing_value);
    field.metadata().set("missing_value_type", "equals");
    // field.metadata().set("missing_value_type", "approximately-equals");
    // field.metadata().set("missing_value_epsilon", 1.);
}

template <typename View>
void gmsh_output(const std::string& name, const Grid& grid, const View& field, const std::string& mask_option,
                 const AtlasTool::Args& args) {
    mpi::Scope mpi_scope("self");
    Log::info() << "Storing field '" << name << ".msh'" << std::endl;
    ATLAS_ASSERT(field.size() <= grid.size());
    std::string coords = args.getString("gmsh.coordinates", "lonlat");
    output::Gmsh gmsh(name + ".msh", Config("coordinates", coords) | Config("ghost", "true"));
    auto mesh = Mesh(grid);
    gmsh.write(mesh);
    auto fs = gridpoints_are_cells(grid) ? FunctionSpace(CellColumns(mesh)) : FunctionSpace(NodeColumns(mesh));

    if (args.has(mask_option)) {
        const std::string mask_name = args.getString(mask_option);
        set_grid_mask(mask_name, fs);
    }
    auto f = fs.createField<double>(option::name(name));
    set_missing_value_metadata(args, f);

    auto v          = array::make_view<double, 1>(f);
    auto g          = array::make_view<gidx_t, 1>(fs.global_index());
    size_t min_size = std::min<size_t>(f.size(), field.size());
    for (size_t i = 0; i < min_size; ++i) {
        v(i) = field[g(i) - 1];
    }
    fs.haloExchange(f);

    // By default, preserve missing values when a mask is present, so we can observe when the
    // mask is not covering the present missing values
    util::Config gmsh_field_default_config;
    if (fs.hasMask()) {
        gmsh_field_default_config.set("missing_value.policy", "preserve");
    }

    gmsh.write(f, gmsh_field_config(args, gmsh_field_default_config));
}

std::vector<std::string> split_string(const std::string& str, char delim) {
    std::vector<std::string> out;
    std::stringstream ss(str);
    std::string item;
    while (std::getline(ss, item, delim)) {
        // trim surrounding whitespace
        auto b = item.find_first_not_of(" \t");
        auto e = item.find_last_not_of(" \t");
        if (b != std::string::npos) {
            out.push_back(item.substr(b, e - b + 1));
        }
    }
    return out;
}

std::vector<int> load_grid_mask(const AtlasTool::Args& args, const std::string& mask_option, const Grid& grid) {
    // Default: everything valid (mask==1). Masks are stored/read in global grid-point order.
    std::vector<int> mask(grid.size(), 1);
    if (args.has(mask_option)) {
        read_mask(args.getString(mask_option), array::mdspan<int, array::dims<1>>(mask.data(), mask.size()));
    }
    return mask;
}

util::IndexKDTree build_lonlat_kdtree(const std::vector<PointLonLat>& lonlat) {
    util::IndexKDTree search;
    search.reserve(lonlat.size());
    for (idx_t i = 0; i < static_cast<idx_t>(lonlat.size()); ++i) {
        search.insert(lonlat[i], i);
    }
    search.build();
    return search;
}

// Extract the interpolation stencil for a small selection of target points and write it to a
// lightweight ASCII file, suitable for studying/visualising a single interpolation with a Python tool.
//
// The stencil of a target point is the row of the (global) interpolation matrix: the set of source
// grid points and their weights that contribute to that target point. Around the stencil we also add
// some 'context' points (nearest neighbours) so that e.g. the source and target masks can be inspected
// in the neighbourhood of the stencil: context source points around the stencil sources, and context
// target points around the target point. When a function space is mesh-based, the mesh cells touching
// the stencil/context are written as well, to visualise the connectivity.
void extract_stencils(const Grid& sgrid, const Grid& tgrid, const Matrix& matrix, const FunctionSpace& src_fs,
                      const FunctionSpace& tgt_fs, const AtlasTool::Args& args) {
    ATLAS_TRACE();

    if (matrix.rows() == 0 || matrix.cols() == 0) {
        Log::warning() << "Stencil extraction skipped: interpolation matrix is empty." << std::endl;
        return;
    }
    ATLAS_ASSERT(static_cast<gidx_t>(matrix.rows()) == static_cast<gidx_t>(tgrid.size()),
                 "Stencil extraction requires the global matrix (rows == target grid size)");
    ATLAS_ASSERT(static_cast<gidx_t>(matrix.cols()) == static_cast<gidx_t>(sgrid.size()),
                 "Stencil extraction requires the global matrix (cols == source grid size)");

    // --- Coordinates (global grid-point order) ---
    std::vector<PointLonLat> src_lonlat(sgrid.size());
    std::vector<PointLonLat> tgt_lonlat(tgrid.size());
    ATLAS_TRACE_SCOPE("Load coordinates") {
        idx_t n = 0;
        for (auto p : sgrid.lonlat()) {
            src_lonlat[n++] = p;
        }
        n = 0;
        for (auto p : tgrid.lonlat()) {
            tgt_lonlat[n++] = p;
        }
    }

    // --- Masks (global grid-point order; default all valid) ---
    const bool has_src_mask = args.has("s.mask");
    const bool has_tgt_mask = args.has("t.mask");
    std::vector<int> src_mask = load_grid_mask(args, "s.mask", sgrid);
    std::vector<int> tgt_mask = load_grid_mask(args, "t.mask", tgrid);

    // --- Context configuration and KDTrees ---
    const long source_context_k = std::max<long>(0, args.getLong("stencil.source-context", 8));
    // Target context defaults to the same amount as the source context when not given explicitly.
    const long target_context_k = std::max<long>(0, args.getLong("stencil.target-context", source_context_k));

    // The target KDTree is needed to select targets by lonlat and to gather target context points.
    util::IndexKDTree tsearch;
    if (target_context_k > 0 || args.has("stencil.target-lonlat")) {
        ATLAS_TRACE_SCOPE("Build target KDTree") {
            tsearch = build_lonlat_kdtree(tgt_lonlat);
        }
    }

    // --- Select target points (result: 0-based row indices) ---
    std::vector<gidx_t> targets;
    if (args.has("stencil.target-index")) {
        for (const auto& s : split_string(args.getString("stencil.target-index"), ',')) {
            gidx_t idx1 = static_cast<gidx_t>(std::stol(s));
            if (idx1 < 1 || idx1 > static_cast<gidx_t>(tgrid.size())) {
                ATLAS_THROW_EXCEPTION("--stencil.target-index=" << idx1 << " is out of range [1," << tgrid.size()
                                                                << "]");
            }
            targets.push_back(idx1 - 1);
        }
    }
    if (args.has("stencil.target-lonlat")) {
        for (const auto& pt : split_string(args.getString("stencil.target-lonlat"), ';')) {
            auto ll = split_string(pt, ',');
            if (ll.size() != 2) {
                ATLAS_THROW_EXCEPTION("--stencil.target-lonlat expects 'lon,lat' pairs separated by ';', got '" << pt
                                                                                                             << "'");
            }
            PointLonLat p{std::stod(ll[0]), std::stod(ll[1])};
            targets.push_back(tsearch.closestPoint(p).payload());
        }
    }
    if (targets.empty()) {
        targets.push_back(tgrid.size() / 2);
        Log::warning() << "Stencil: neither --stencil.target-lonlat nor --stencil.target-index given; defaulting to "
                          "target global index "
                       << (targets.front() + 1) << std::endl;
    }
    // De-duplicate while preserving order
    {
        std::unordered_set<gidx_t> seen;
        std::vector<gidx_t> unique_targets;
        for (auto t : targets) {
            if (seen.insert(t).second) {
                unique_targets.push_back(t);
            }
        }
        targets.swap(unique_targets);
    }

    // --- KDTree over source grid points for context lookup ---
    util::IndexKDTree ssearch;
    if (source_context_k > 0) {
        ATLAS_TRACE_SCOPE("Build source KDTree") {
            ssearch = build_lonlat_kdtree(src_lonlat);
        }
    }

    // --- Optional source/target mesh connectivity ---
    auto mesh_of = [](const FunctionSpace& fs) -> Mesh {
        if (auto nc = NodeColumns(fs)) {
            return nc.mesh();
        }
        if (auto cc = CellColumns(fs)) {
            return cc.mesh();
        }
        return Mesh();
    };
    Mesh src_mesh = mesh_of(src_fs);
    Mesh tgt_mesh = mesh_of(tgt_fs);

    auto matrix_view = linalg::make_host_view<double, int>(matrix);

    const std::string filename = args.getString("stencil.file", "stencil.txt");
    std::ofstream out(filename);
    if (!out) {
        ATLAS_THROW_EXCEPTION("Could not open stencil output file '" << filename << "'");
    }
    out << std::setprecision(10);

    // --- Header ---
    out << "# tool=atlas-interpolation-stencil\n";
    out << "# format=1\n";
    out << "# source_grid=" << sgrid.name() << "\n";
    out << "# target_grid=" << tgrid.name() << "\n";
    out << "# interpolation=" << get_interpolation_method(args) << "\n";
    out << "# source_mask=" << (has_src_mask ? "present" : "none") << "\n";
    out << "# target_mask=" << (has_tgt_mask ? "present" : "none") << "\n";
    out << "# source_context=" << source_context_k << "\n";
    out << "# target_context=" << target_context_k << "\n";
    out << "# source_mesh=" << (src_mesh ? "present" : "none") << "\n";
    out << "# target_mesh=" << (tgt_mesh ? "present" : "none") << "\n";
    out << "# n_stencils=" << targets.size() << "\n";
    out << "#\n";
    out << "# Record types (whitespace-separated, first token identifies the record):\n";
    out << "#   T gidx lon lat mask                    : target point\n";
    out << "#   S gidx lon lat weight mask             : source point in the stencil (has a weight)\n";
    out << "#   X gidx lon lat mask                    : context source point (no weight, for context only)\n";
    out << "#   Y gidx lon lat mask                    : context target point (target neighbour, for context)\n";
    out << "#   C n gidx0 gidx1 ... gidx{n-1}          : source mesh cell (source node global indices)\n";
    out << "#   D n gidx0 gidx1 ... gidx{n-1}          : target mesh cell (target node global indices)\n";
    out << "# Global indices (gidx) are 1-based, within their own (source or target) grid.\n";
    out << "# Each stencil is delimited by 'BEGIN <i>' / 'END <i>'.\n";
    out << "#\n";

    // Emit the mesh cells that touch any of the 'involved' node global indices (1-based) as records
    // with the given tag ('C' for source, 'D' for target).
    auto emit_mesh_cells = [&out](char tag, const Mesh& mesh, const std::unordered_set<gidx_t>& involved) {
        if (!mesh || involved.empty()) {
            return;
        }
        const auto& conn = mesh.cells().node_connectivity();
        auto node_gidx   = array::make_view<gidx_t, 1>(mesh.nodes().global_index());
        for (idx_t e = 0; e < mesh.cells().size(); ++e) {
            const idx_t nb_cell_nodes = conn.cols(e);
            bool touches              = false;
            for (idx_t ni = 0; ni < nb_cell_nodes; ++ni) {
                if (involved.count(node_gidx(conn(e, ni)))) {
                    touches = true;
                    break;
                }
            }
            if (!touches) {
                continue;
            }
            out << tag << " " << nb_cell_nodes;
            for (idx_t ni = 0; ni < nb_cell_nodes; ++ni) {
                out << " " << node_gidx(conn(e, ni));
            }
            out << "\n";
        }
    };

    for (std::size_t is = 0; is < targets.size(); ++is) {
        const gidx_t row = targets[is];  // 0-based target grid index

        // Stencil: source columns + weights for this target row
        std::vector<gidx_t> stencil_cols;  // 0-based source grid indices
        std::vector<double> stencil_wgts;
        linalg::sparse_matrix_for_each_row(static_cast<std::size_t>(row), matrix_view,
                                           [&](int /*r*/, int col, double val) {
                                               stencil_cols.push_back(col);
                                               stencil_wgts.push_back(val);
                                           });

        // Members of stencil (1-based gidx) for quick lookup
        std::unordered_set<gidx_t> stencil_set;
        for (auto c : stencil_cols) {
            stencil_set.insert(c + 1);
        }

        // Context: union of nearest source neighbours around each stencil point, excluding stencil points
        std::vector<gidx_t> context_cols;  // 0-based
        std::unordered_set<gidx_t> context_set;
        if (source_context_k > 0) {
            for (auto c : stencil_cols) {
                auto neighbours = ssearch.closestPoints(src_lonlat[c], source_context_k).payloads();
                for (auto nb : neighbours) {
                    const gidx_t g = nb + 1;
                    if (stencil_set.count(g) == 0 && context_set.insert(g).second) {
                        context_cols.push_back(nb);
                    }
                }
            }
        }

        // All source global indices (1-based) involved in this stencil (for cell filtering)
        std::unordered_set<gidx_t> involved = stencil_set;
        involved.insert(context_set.begin(), context_set.end());

        // Context: nearest target neighbours around the target point (excluding the target itself)
        std::vector<gidx_t> tgt_context;  // 0-based target grid indices
        if (target_context_k > 0) {
            auto neighbours = tsearch.closestPoints(tgt_lonlat[row], target_context_k + 1).payloads();
            for (auto nb : neighbours) {
                if (nb != row) {
                    tgt_context.push_back(nb);
                }
            }
        }
        // Target global indices (1-based) involved (target + its context) for target-cell filtering
        std::unordered_set<gidx_t> tgt_involved;
        tgt_involved.insert(row + 1);
        for (auto t : tgt_context) {
            tgt_involved.insert(t + 1);
        }

        out << "BEGIN " << is << "\n";

        // Target
        out << "T " << (row + 1) << " " << tgt_lonlat[row].lon() << " " << tgt_lonlat[row].lat() << " "
            << tgt_mask[row] << "\n";

        // Context target points (no weight)
        for (auto t : tgt_context) {
            out << "Y " << (t + 1) << " " << tgt_lonlat[t].lon() << " " << tgt_lonlat[t].lat() << " " << tgt_mask[t]
                << "\n";
        }

        if (stencil_cols.empty()) {
            out << "# (empty stencil: target has no source contributions / missing value)\n";
        }

        // Stencil source points (with weights)
        for (std::size_t k = 0; k < stencil_cols.size(); ++k) {
            const gidx_t c = stencil_cols[k];
            out << "S " << (c + 1) << " " << src_lonlat[c].lon() << " " << src_lonlat[c].lat() << " " << stencil_wgts[k]
                << " " << src_mask[c] << "\n";
        }

        // Context source points (no weight)
        for (auto c : context_cols) {
            out << "X " << (c + 1) << " " << src_lonlat[c].lon() << " " << src_lonlat[c].lat() << " " << src_mask[c]
                << "\n";
        }

        // Source mesh cells touching any involved source point
        emit_mesh_cells('C', src_mesh, involved);

        // Target mesh cells touching the target point or its context
        emit_mesh_cells('D', tgt_mesh, tgt_involved);

        out << "END " << is << "\n";
    }

    Log::info() << "Stencil extraction: wrote " << targets.size() << " stencil(s) to '" << filename << "'" << std::endl;
}

void test_matrix(const Grid& sgrid, const Grid& tgrid, const Matrix& matrix, const AtlasTool::Args& args,
                 const SourceInitialCondition& initial_condition) {
    ATLAS_TRACE();

    // Allocate and initialise own memory here to show possibilities
    std::vector<double> sdata(sgrid.size());
    std::vector<double> tdata(tgrid.size());

    ATLAS_TRACE_SCOPE("initialize source") {
        if (args.has("init") || not args.has("s.grib")) {
            ATLAS_DEBUG("using init");
            idx_t n{0};
            for (auto p : sgrid.lonlat()) {
                sdata[n] = initial_condition(n, p);
                ++n;
            }
        }
        else {
            ATLAS_DEBUG("using grib");

            ATLAS_ASSERT(args.has("s.grib"));
#if ATLAS_HAVE_GRIB
            GribFileReader grib_reader(args.getString("s.grib"));
            grib_reader.grib().get_values(sdata.data(), sdata.size());
#endif
        }
    }

    ATLAS_TRACE_SCOPE("sparse-matrix-multiply") {
        StopWatch timer_serial_sparse_matrix_multiply;
        timer_serial_sparse_matrix_multiply.start();
        atlas::linalg::sparse_matrix_multiply(atlas::linalg::make_host_view<double>(matrix),
                                              atlas::array::make_view<double, 1>(sdata.data(), matrix.cols()),
                                              atlas::array::make_view<double, 1>(tdata.data(), matrix.rows()));
        timer_serial_sparse_matrix_multiply.stop();
        Log::info() << "Serial sparse-matrix-multiply timer  \t: "
                    << elapsed_ms(timer_serial_sparse_matrix_multiply.elapsed(), true) << " [ms]" << std::endl;
        Log::info() << "Serial sparse-matrix non-zero entries\t: " << matrix.nnz() << std::endl;
    }

    ATLAS_TRACE_SCOPE("Apply missing values to target") {
        double missing_value = std::numeric_limits<double>::lowest();
        args.get("missing_value", missing_value);
        auto m = atlas::linalg::make_host_view<double>(matrix);
        auto is_matrix_row_empty = [&m](idx_t r) -> bool {
            return m.outer()[r + 1] - m.outer()[r] == 0;
        };
        for (idx_t r = 0; r < m.rows(); ++r) {
            if (is_matrix_row_empty(r)) {
                tdata[r] = missing_value;
            }
        }
    }

    if (args.getBool("output-checksum", false)) {
        ATLAS_TRACE("checksum target-field");

        // The matrix version is possibly not bit-identical with the interpolator version
        // This is because the matrix version uses grid points without halo, wrapping around the globe
        size_t min_size = std::min(tgt_ref.size(), tdata.size());
        for (size_t i = 0; i < min_size; ++i) {
            if (std::abs(tgt_ref[i] - tdata[i]) < 1.e-14) {
                tdata[i] = tgt_ref[i];
            }
        }

        std::ofstream checksum_file(args.getString("checksum.file", "checksum"), std::ios::app);
        auto target_checksum = atlas::util::checksum(tdata.data(), matrix.rows());
        checksum_file << std::setw(8) << target_checksum << "    [test-matrix]   checksum of target field" << std::endl;
    }

#if 0
    {
        auto it = tgrid.lonlat().begin();
        for (size_t i=0; i<tdata.size(); ++i) {
            if (tdata[i] < -2.0) {
                std::vector<double> src_val;
                std::vector<double> weights;
                linalg::sparse_matrix_for_each_row(i, atlas::linalg::make_host_view<double>(matrix), [&](auto row, auto col, auto val) {
                    weights.emplace_back(val);
                    src_val.emplace_back(sdata[col]);
                });

                Log::info() << i << "\t" << tdata[i] << "\t" << *it << std::endl;
                Log::info() << "\t weights " << weights << std::endl; 
                Log::info() << "\t src_val " << src_val << std::endl; 
            }
            ++it;
        }
    }
#endif

    if (args.getBool("output-gmsh", false)) {
        ATLAS_TRACE_SCOPE("Gmsh serial output") {
            if (args.getBool("output-gmsh-source", false)) {
                std::string sname = "test_matrix_source_" + get_matrix_name(args);
                auto sview        = atlas::array::make_view<double, 1>(sdata.data(), matrix.cols());
                gmsh_output(sname, sgrid, sview, "s.mask", args);
            }

            std::string tname = "test_matrix_target_" + get_matrix_name(args);
            auto tview        = atlas::array::make_view<double, 1>(tdata.data(), matrix.rows());
            gmsh_output(tname, tgrid, tview, "t.mask", args);
        }
    }
}

int AtlasInterpolations::execute(const AtlasTool::Args& args) {
    ATLAS_TRACE("main");

    constexpr int preserve_existing_value = 0; // when 0 it preserves existing values
    if (::setenv("ATLAS_INTERPOLATION_FALLBACK", "nearest-neighbour", preserve_existing_value) != 0) {
        ATLAS_THROW_EXCEPTION("Could not set ATLAS_INTERPOLATION_FALLBACK");
    }

    bool invalid_option = false;
    for (const std::string& arg : args) {
        if (!arg.empty() && arg.front() == '-') {
            Log::error() << "Invalid option '" << arg << "'. Options must start with '--'." << std::endl;
            invalid_option = true;
        }
    }
    if (invalid_option) {
        return failed();
    }

    if (args.has("list")) {
        Log::info() << "Usage: " << usage() << std::endl;
        Log::info() << "Available interpolators:" << std::endl;
        std::set<std::string> interpolation_types;
        for (const auto& key : interpolation::MethodFactory::keys()) {
            interpolation_types.insert(key);
        }
        for (const auto& b : interpolation_types) {
            if (b == "conservative-spherical-polygon") {
                Log::info() << "\t" << "cons" << '\n';
                Log::info() << "\t" << "cons2" << '\n';
            }
            else {
                Log::info() << "\t" << b << '\n';
            }
        }
        return success();
    }

    bool output_checksum           = args.getBool("output-checksum", false);
    std::string checksum_file_path = args.getString("checksum.file", "checksum");
    if (output_checksum) {
        eckit::PathName checksum_file(checksum_file_path);
        if (checksum_file.exists()) {
            checksum_file.unlink();
        }
    }

    auto sgrid             = get_sgrid(args);
    auto tgrid             = get_tgrid(args);
    const SourceInitialCondition initial_condition(args, sgrid);
    std::string sdata_name = gridpoints_are_cells(sgrid) ? "cell data" : "node data";
    std::string tdata_name = gridpoints_are_cells(tgrid) ? "cell data" : "node data";

    std::string interpolation_method = get_interpolation_method(args);
    Log::info() << "\tsource grid\t:\t" << sgrid.name() << " (" << sdata_name << ")\n\ttarget grid\t:\t" << tgrid.name()
                << " (" << tdata_name << ")\n\tinterpolation\t:\t" << interpolation_method << std::endl;

    std::string matrix_format = get_matrix_format(args);
    std::string matrix_name   = get_matrix_name(args);
    bool matrix_tested        = false;
    Matrix matrix;
    if (args.getBool("read-matrix", false)) {
        if (mpi::comm().rank() == 0) {
            ATLAS_TRACE_SCOPE("Read matrix on rank 0") {
                timers.global_matrix_read.start();
                matrix = read_matrix(matrix_name, matrix_format);
                timers.global_matrix_read.stop();
                Log::info() << "Global matrix read timer\t: " << elapsed_ms(timers.global_matrix_read, true) << " [ms]"
                            << std::endl;
                print_matrix(matrix, Log::info());
                Log::info() << std::endl;
                if (output_checksum) {
                    ATLAS_TRACE("checksum input-matrix");
                    std::ofstream checksum_file(checksum_file_path, std::ios::app);
                    auto matrix_view    = linalg::make_host_view<double, int>(matrix);
                    auto outer_checksum = util::checksum(matrix_view.inner(), matrix_view.inner_size());
                    auto inner_checksum = util::checksum(matrix_view.outer(), matrix_view.outer_size());
                    auto value_checksum = util::checksum(matrix_view.value(), matrix_view.value_size());
                    checksum_file << std::setw(8) << outer_checksum << "    [read-matrix]   checksum of matrix.outer()"
                                  << std::endl;
                    checksum_file << std::setw(8) << inner_checksum << "    [read-matrix]   checksum of matrix.inner()"
                                  << std::endl;
                    checksum_file << std::setw(8) << value_checksum << "    [read-matrix]   checksum of matrix.value()"
                                  << std::endl;
                }
            }
        }
        if (args.getBool("test-matrix", false)) {
            if (mpi::comm().rank() == 0) {
                test_matrix(sgrid, tgrid, matrix, args, initial_condition);
            }
            matrix_tested = true;
        }
    }

    FunctionSpace src_fs;
    FunctionSpace tgt_fs;
    Interpolation interpolator;

    ATLAS_TRACE_SCOPE("Setup source and target") {
        timers.functionspace_setup.start();
        std::tie(src_fs, tgt_fs) = get_fs(sgrid, tgrid, args);
        timers.functionspace_setup.stop();
    }

    interpolation::Cache cache;
    if (args.getBool("read-matrix", false)) {
        mpi::comm().barrier();
        timers.global_matrix_setup.start();
        auto distributed_matrix = interpolation::distribute_global_matrix(src_fs, tgt_fs, matrix);
        timers.global_matrix_setup.stop();
        cache = interpolation::MatrixCache{std::move(distributed_matrix)};
        Log::info() << "Global matrix setup timer\t: " << elapsed_ms(timers.global_matrix_setup) << " [ms]"
                    << std::endl;
    }

    matrix.clear();

    if (args.getBool("interpolate", false) || args.getBool("output-matrix", false) ||
        args.getBool("test-interpolator", false) || args.getBool("output-stencil", false) ||
        (args.getBool("test-matrix", false) && !matrix_tested)) {

        auto configure_mask = [&](const std::string& mask_option, FunctionSpace& fs) {
            if (args.has(mask_option)) {
                set_grid_mask(args.getString(mask_option), fs);
            }
        };

        configure_mask("s.mask", src_fs);
        configure_mask("t.mask", tgt_fs);

        ATLAS_TRACE_SCOPE("Setup interpolator") {
            timers.interpolation_setup.start();
            auto config = get_interpolation_config(sgrid, tgrid, args);
            if (string_starts_with(interpolation_method, "grid-box")) {
                interpolator = Interpolation(config, sgrid, tgrid);
            }
            else {
                interpolator = Interpolation(config, src_fs, tgt_fs, cache);
            }
            timers.interpolation_setup.stop();
            ATLAS_DEBUG_VAR(interpolator.failedInterpolations().size());
        }

        Log::info() << "Grid + FunctionSpace timer\t: " << elapsed_ms(timers.functionspace_setup) << " [ms]"
                    << std::endl;
        Log::info() << "Interpolation setup timer\t: " << elapsed_ms(timers.interpolation_setup) << " [ms]"
                    << std::endl;

        if (args.getBool("test-interpolator", false)) {
            ATLAS_TRACE("Test test_interpolator");
            auto src_field  = interpolator.source().createField<double>();
            auto tgt_field  = interpolator.target().createField<double>();
            auto src_lonlat = array::make_view<double, 2>(interpolator.source().lonlat());
            auto src_global_index = array::make_view<gidx_t, 1>(interpolator.source().global_index());
            ATLAS_TRACE_SCOPE("initialize source") {
                auto src_field_v = array::make_view<double, 1>(src_field);
                for (idx_t i = 0; i < src_field_v.size(); ++i) {
                    src_field_v[i] = initial_condition(
                        src_global_index(i) - 1, PointLonLat{src_lonlat(i, 0), src_lonlat(i, 1)});
                }
            }
            src_field.set_dirty(false);  // all values are up to date

            timers.interpolation_exe.start();
            interpolator.execute(src_field, tgt_field);
            timers.interpolation_exe.stop();
            Log::info() << "Interpolation execute timer     : " << elapsed_ms(timers.interpolation_exe) << " [ms]"
                        << std::endl;

#if 0
            // Set the target values marked as missing to zero
            auto tgt_field_v = array::make_view<double,1>(tgt_field);
            for (idx_t i = 0; i < tgt_fs.size(); ++i) {
                if (tgt_field_v[i] == 9999.) {
                    tgt_field_v[i] = 0;
                }
            }
#endif

            Field tgt_field_global;
            if (output_checksum) {
                ATLAS_TRACE("checksum target-field");
                tgt_field_global        = tgt_fs.createField(tgt_field, option::global());
                auto tgt_field_global_v = array::make_view<double, 1>(tgt_field_global);
                tgt_fs.gather(tgt_field, tgt_field_global);
#if 0
                for( idx_t i=0; i<tgt_field_global_v.size(); ++i) {
                    if (tgt_field_global_v(i) == 9999.) {
                        tgt_field_global_v(i) = 0.;
                    }
                }
#endif

                if (mpi::rank() == 0) {
                    auto target_checksum = util::checksum(tgt_field_global_v.data(), tgt_field_global_v.size());
                    std::ofstream checksum_file(checksum_file_path, std::ios::app);
                    checksum_file << std::setw(8) << target_checksum << "    [test-interp]   checksum of target field"
                                  << std::endl;

                    // for further comparison with test-matrix
                    tgt_ref.resize(tgt_field_global_v.size());
                    for (idx_t i = 0; i < tgt_field_global_v.size(); ++i) {
                        tgt_ref[i] = tgt_field_global_v[i];
                    }
                }
            }

            if (args.getBool("output-gmsh", false)) {
                Mesh tmesh;
                if (auto fs = NodeColumns(tgt_fs)) {
                    tmesh = fs.mesh();
                }
                if (auto fs = CellColumns(tgt_fs)) {
                    tmesh = fs.mesh();
                }
                if (not tmesh) {
                    Log::warning()
                        << "Target function space is NOT mesh based, hence no Gmsh output will be created for this"
                        << std::endl;
                }
                else {
                    tgt_field.haloExchange();
                    set_missing_value_metadata(args, tgt_field);
                    std::string tgt_name = "test_interpolator_target_" + matrix_name;
                    Log::info() << "Storing field '" << tgt_name << ".msh'." << std::endl;
                    std::string coords = args.getString("gmsh.coordinates", "lonlat");
                    output::Gmsh gmsh(tgt_name + ".msh", Config("coordinates", coords) | Config("ghost", "false"));
                    gmsh.write(tmesh);
                    gmsh.write(tgt_field, gmsh_field_config(args));
                }
            }
        }
    }

    if (args.getBool("output-matrix", false) || args.getBool("output-stencil", false) ||
        (args.getBool("test-matrix", false) && !matrix_tested)) {
        ATLAS_TRACE("Matrix operations");
        ATLAS_TRACE_SCOPE("assemble_global_matrix") {
            matrix = interpolation::assemble_global_matrix(tgrid.size(), sgrid.size(), interpolator);
        }
        if (mpi::comm().rank() == 0) {
            if (args.getBool("output-stencil", false)) {
                extract_stencils(sgrid, tgrid, matrix, src_fs, tgt_fs, args);
            }
            if (args.getBool("output-matrix", false)) {
                write_matrix(matrix, matrix_name, matrix_format);
                print_matrix(matrix, Log::info());
                if (output_checksum) {
                    ATLAS_TRACE("checksum output-matrix");
                    std::ofstream checksum_file(checksum_file_path, std::ios::app);
                    auto matrix_view    = linalg::make_host_view<double, int>(matrix);
                    auto outer_checksum = util::checksum(matrix_view.inner(), matrix_view.inner_size());
                    auto inner_checksum = util::checksum(matrix_view.outer(), matrix_view.outer_size());
                    auto value_checksum = util::checksum(matrix_view.value(), matrix_view.value_size());
                    checksum_file << std::setw(8) << outer_checksum << "    [output-matrix] checksum of matrix.outer()"
                                  << std::endl;
                    checksum_file << std::setw(8) << inner_checksum << "    [output-matrix] checksum of matrix.inner()"
                                  << std::endl;
                    checksum_file << std::setw(8) << value_checksum << "    [output-matrix] checksum of matrix.value()"
                                  << std::endl;
                }
            }
            if (args.getBool("test-matrix", false) && !matrix_tested) {
                test_matrix(sgrid, tgrid, matrix, args, initial_condition);
            }
        }
    }

    mpi::comm().barrier();
    if (mpi::rank() == 0) {
        if (output_checksum) {
            auto checksum_file = eckit::PathName(checksum_file_path);
            if (checksum_file.exists()) {
                std::ifstream f(checksum_file);
                if (f.is_open()) {
                    Log::info() << "\nContent of checksum file " << checksum_file << std::endl;
                    Log::info() << f.rdbuf() << std::endl;
                }
            }
        }
    }

    return success();
}

//-----------------------------------------------------------------------------

Matrix AtlasInterpolations::read_matrix(std::string matrix_name, std::string format) {
    if (format == "eckit") {
        Log::info() << "Reading matrix from file '" << matrix_name << ".eckit'" << std::endl;
        eckit::linalg::SparseMatrix eckit_matrix;
        eckit_matrix.load(matrix_name + ".eckit");
        return atlas::linalg::make_sparse_matrix_storage(std::move(eckit_matrix));
    }
    else if (format == "scrip") {
        Log::info() << "Reading matrix from file '" << matrix_name << ".nc'" << std::endl;
        return ScripIO::read_matrix(matrix_name + ".nc");
    }
    else {
        ATLAS_THROW_EXCEPTION("Matrix format " << format << " is not recognised. Recognized are {eckit,scrip}");
    }
    return Matrix{};
}


void AtlasInterpolations::write_matrix(const Matrix& matrix, std::string matrix_name, std::string format) {
    ATLAS_TRACE();
    if (format == "eckit") {
        Log::info() << "Writing matrix in eckit format to file '" << matrix_name << ".eckit'" << std::endl;
        auto eckit_matrix = linalg::make_non_owning_eckit_sparse_matrix(matrix);
        eckit_matrix.save(matrix_name + ".eckit");
    }
    else if (format == "scrip") {
        Log::info() << "Writing matrix in scrip format to file '" << matrix_name << ".nc'" << std::endl;
        ScripIO::write_matrix(matrix, matrix_name + ".nc");
    }
    else {
        ATLAS_THROW_EXCEPTION("Matrix format " << format << " is not recognised. Recognized are {eckit,scrip}");
    }
}

void AtlasInterpolations::print_matrix(const Matrix& matrix, std::ostream& out) {
    out << "    matrix.rows : " << matrix.rows() << '\n';
    out << "    matrix.cols : " << matrix.cols() << '\n';
    out << "    matrix.nnz  : " << matrix.nnz() << '\n';
}


}  // namespace atlas


int main(int argc, char** argv) {
    atlas::AtlasInterpolations tool(argc, argv);
    return tool.start();
}
