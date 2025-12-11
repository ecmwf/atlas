/*
 * (C) Copyright 2025 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <fstream>
#include <iomanip>

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
#include "atlas/linalg/sparse.h"
#include "atlas/linalg/sparse/MakeEckitSparseMatrix.h"
#include "atlas/linalg/sparse/SparseMatrixToTriplets.h"
#include "atlas/linalg/sparse/SparseMatrixTriplet.h"


#include "atlas/mesh/Mesh.h"
#include "atlas/meshgenerator.h"
#include "atlas/output/Gmsh.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/runtime/AtlasTool.h"
#include "atlas/util/CoordinateEnums.h"
#include "atlas/util/Checksum.h"
#include "atlas/util/function/VortexRollup.h"


#include "atlas/util/function/MDPI_functions.h"
#include "atlas/util/function/SolidBodyRotation.h"
#include "atlas/util/function/SphericalHarmonic.h"
#include "atlas/util/function/VortexRollup.h"

#include "AtlasIO.h"
#include "ScripIO.h"

#if ATLAS_HAVE_GRIB
#include "GribFileReader.h"
#include "Grib.h"
#endif

using atlas::functionspace::PointCloud;
using atlas::functionspace::NodeColumns;
using atlas::functionspace::CellColumns;
using atlas::functionspace::StructuredColumns;
using atlas::util::Config;

using Matrix = atlas::linalg::SparseMatrixStorage;
using StopWatch = atlas::runtime::trace::StopWatch;


namespace atlas::linalg {
    
template <typename Value, typename Index>
class TripletView {
public:
    TripletView(Index& row, Index& column, Value& value):
        row_{row}, column_{column}, value_{value} {}
    Index& row() { return row_; }
    Index& column() { return column_; }
    Value& value() { return value_; }
    const Index& row() const { return row_; }
    const Index& column() const { return column_; }
    const Value& value() const { return value_; }
private:
    Value& value_;
    Index& row_;
    Index& column_;
};
template <typename Value, typename Index>
class TripletsView {
    public:
    class iterator {
        public:
        iterator(std::size_t idx, TripletsView* parent): idx_{idx}, parent_{parent} {}
        iterator& operator++() { ++idx_; return *this; }
        bool operator!=(const iterator& other) const { return idx_ != other.idx_; }
        TripletView<Value, Index> operator*() const {
            return TripletView<Value, Index>{parent_->row(idx_), parent_->column(idx_), parent_->value(idx_)};
        }
        private:
        std::size_t idx_;
        TripletsView* parent_;
    };
    class const_iterator {
    public:
        const_iterator(std::size_t idx, const TripletsView* parent): idx_{idx}, parent_{parent} {}
        const_iterator& operator++() { ++idx_; return *this; }
        bool operator!=(const const_iterator& other) const { return idx_ != other.idx_; }
        TripletView<const Value, const Index> operator*() const {
            return TripletView<const Value, const Index>{parent_->row(idx_), parent_->column(idx_), parent_->value(idx_)};
        }
    private:
        std::size_t idx_;
        const TripletsView* parent_;
    };

    template <typename A1, typename A2, typename A3>
    TripletsView(std::vector<Index, A1>& rows, std::vector<Index, A2>& columns, std::vector<Value,A3>& values) :
        TripletsView(rows.data(), columns.data(), values.data(), values.size()) {
            ATLAS_ASSERT(rows.size() == columns.size());
            ATLAS_ASSERT(values.size() == rows.size());
        }

    TripletsView(Index* rows, Index* columns, Value* values, std::size_t size):
        rows_{rows}, columns_{columns}, values_{values}, size_{size} {}

    iterator begin() { return iterator{0, this}; }  
    iterator end() { return iterator{size_, this}; }
    const_iterator begin() const { return const_iterator{0, this}; }  
    const_iterator end() const { return const_iterator{size_, this}; }
    std::size_t size() const { return size_; }
    Index& row(std::size_t idx) { return rows_[idx]; }
    Index& column(std::size_t idx) { return columns_[idx]; }
    Value& value(std::size_t idx) { return values_[idx]; }
    const Index& row(std::size_t idx) const { return rows_[idx]; }
    const Index& column(std::size_t idx) const { return columns_[idx]; }
    const Value& value(std::size_t idx) const { return values_[idx]; }
private:
    Index* rows_;
    Index* columns_;
    Value* values_;
    std::size_t size_;
};


interpolation::MatrixCache create_matrix_cache(Interpolation interpolation) {
    interpolation::MatrixCache interpolation_cache = interpolation.createCache();
    auto interpolation_matrix_storage = interpolation_cache.matrix();
    auto interpolation_matrix = atlas::linalg::make_host_view<eckit::linalg::Scalar,eckit::linalg::Index>(interpolation_matrix_storage);

    if (interpolation.failedInterpolations().size() == 0) {
        return interpolation_cache;
    }

    Interpolation fallback(option::type("nearest-neighbour"), interpolation.source(), interpolation.target());
    interpolation::MatrixCache fallback_interpolation_cache = interpolation.createCache();
    auto fallback_interpolation_matrix_storage = interpolation_cache.matrix();
    auto fallback_interpolation_matrix = atlas::linalg::make_host_view<eckit::linalg::Scalar,eckit::linalg::Index>(fallback_interpolation_matrix_storage);

    auto missing_rows = array::make_view<int,1>(interpolation.failedInterpolations());
    std::set<idx_t> missing_set;
    for( idx_t i=0; i<missing_rows.size(); ++i ) {
        missing_set.insert( missing_rows(i) );
    }

    std::vector<Triplet<eckit::linalg::Scalar, eckit::linalg::Index>> triplets;
    triplets.reserve(interpolation_matrix.rows());
    for (std::size_t r=0; r<interpolation_matrix.rows(); ++r) {
        if (missing_set.find(r) != missing_set.end()) {
            for (idx_t c = interpolation_matrix.outer()[r]; c < interpolation_matrix.outer()[r + 1]; ++c) {
                triplets.emplace_back(r, interpolation_matrix.inner()[c], interpolation_matrix.value()[c]);
            }
        }
        else {
            for (idx_t c = fallback_interpolation_matrix.outer()[r]; c < fallback_interpolation_matrix.outer()[r + 1]; ++c) {
                triplets.emplace_back(r, fallback_interpolation_matrix.inner()[c], fallback_interpolation_matrix.value()[c]);
            }
        }
    }
    return interpolation::MatrixCache{
        make_sparse_matrix_storage_from_triplets(
            interpolation_matrix.rows(),
            interpolation_matrix.cols(),
            triplets,
            true
        )
    };
}


}


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
        add_option(new SimpleOption<bool>("list",  "List available interpolation types"));
        add_option(new SimpleOption<bool>("interpolate", "Test an interpolation type"));

        // Configuration
        add_option(new eckit::option::Separator("Interpolation options"));
        add_option(new SimpleOption<std::string>("i.knn", "number of nearest neighbours in case i.type=knn or k-nearest-neighbours. (default=1)"));
        add_option(new SimpleOption<bool>("i.normalise", "Normalise weights"));
        add_option(new SimpleOption<std::string>("s.grid", "source grid"));
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
        add_option(new SimpleOption<std::string>("matrix.name", "Name of the remapping matrix. If not provided a default unique name will be chosen"));
        add_option(new SimpleOption<std::string>("matrix.format", "Format of the remapping matrix: eckit, SCRIP"));

        add_option(new SimpleOption<bool>("output-gmsh", "Set to enable gmsh output for tests"));
        add_option(new SimpleOption<std::string>("gmsh.coordinates", "Choose the coordinates in gmsh output: {lonlat, xyz}"));

        add_option(new SimpleOption<bool>("output-checksum", "Compute checksums"));
        add_option(new SimpleOption<std::string>("checksum.file", "Path of files for checksums [default='checksum']"));


        // Initial condition options
        add_option(new eckit::option::Separator("Initial condition options"));
        add_option(new SimpleOption<std::string>(
            "init", "Setup initial source field [ constant, spherical_harmonic, vortex_rollup (default), solid_body_rotation_wind_magnitude ]"));
        add_option(new SimpleOption<double>("solid_body_rotation.angle", "Angle of solid body rotation (default = 0.)"));
        add_option(new SimpleOption<double>("vortex_rollup.t", "Value that controls vortex rollup (default = 0.5)"));
        add_option(new SimpleOption<double>("constant.value", "Value that is assigned in case init==constant)"));
        add_option(new SimpleOption<long>("spherical_harmonic.n", "total wave number 'n' of a spherical harmonic"));
        add_option(new SimpleOption<long>("spherical_harmonic.m", "zonal wave number 'm' of a spherical harmonic"));
        add_option(new SimpleOption<std::string>("s.grib", "file path to input"));

        add_option(new SimpleOption<double>("missing-value", "Value used with --test-matrix to add as missing value"));
        add_option(new SimpleOption<bool>("gmsh.missing-value", "Treat missing-value in Gmsh output with --test-matrix"));

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

double elapsed_ms(const StopWatch& timer, bool single_rank = false) {
    double timer_elapsed = timer.elapsed() * 1000.;
    if (single_rank) {
        return timer_elapsed;
    }
    mpi::comm().allReduceInPlace(&timer_elapsed, 1, eckit::mpi::sum());
    return timer_elapsed / mpi::size();
}

std::string get_interpolation_method(const AtlasTool::Args& args) {
    std::string method = args.getString("i.type", "finite-element");
    auto method_contains = [method](std::string_view substr) { return method.find(substr) != std::string_view::npos; };
    if (method_contains("cons")) {
        return method_contains("2") ? "cons2" : "cons";
    }
    return method;
}

std::string get_extension(const std::string& str) {
   std::size_t i = str.rfind('.', str.length());
   if (i != std::string::npos) {
      return(str.substr(i+1, str.length() - i));
   }
   return("");
}

std::string get_basename(const std::string& str) {
    std::string basename = str;
    std::string ext = get_extension(str);
    if (ext == "eckit" || ext == "nc") {
        std::size_t lastindex = str.find_last_of(".");
        basename = str.substr(0, lastindex);
    }
    return basename;
}

bool string_contains(std::string_view str, std::string_view substr) {
    return str.find(substr) != std::string_view::npos;
};

bool string_starts_with(std::string_view str, std::string_view substr) {
    return str.rfind(substr, 0) == 0;
};

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


Grid get_sgrid(const AtlasTool::Args& args) {
    ATLAS_ASSERT(args.has("s.grib") or args.has("s.grid"), "Could not detect grid because s.grib or s.grid is not defined");
    if (args.has("s.grid")) {
        return Grid(args.getString("s.grid"));
    }
    else {
        ATLAS_ASSERT(ATLAS_HAVE_GRIB);
#if ATLAS_HAVE_GRIB
        // Read from GRIB
        GribFileReader grib_reader(args.getString("s.grib"));
        ATLAS_ASSERT(grib_reader.count());
        return grib_reader.grib().get_grid();
#else
        return Grid();
#endif
    }
}

Grid get_tgrid(const AtlasTool::Args& args) {
    return Grid(args.getString("t.grid", "O32"));
}

std::string get_matrix_name(const AtlasTool::Args& args) {
    if (args.has("matrix.name")) {
        auto matrix = args.getString("matrix.name");
        return get_basename(matrix);
    }
    auto sgrid = get_sgrid(args).name();
    auto tgrid = args.getString("t.grid");
    auto interpolation_name = get_interpolation_method(args);
    return "remap_" + sgrid + "_" + tgrid + "_" + interpolation_name;
}

bool gridpoints_are_cells(const Grid& g) {
    if (g.type() == "healpix") {
        return true;
    }
    return false;
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
            ATLAS_THROW_EXCEPTION("Could not detect the mesh generator for grid \"" << grid.name() << "\" that treats grid points as nodes");
        }
    }
    return Mesh(grid, dist, option::halo(halo));
}

std::pair<FunctionSpace,FunctionSpace> get_fs(const Grid& sgrid, const Grid& tgrid, const AtlasTool::Args& args) {
    enum class fs_type {
        unset,
        structured,
        cells,
        nodes,
        points
    };
    auto to_fs_type = [](const std::string& str) -> fs_type {
        if (str == "structured") { return fs_type::structured; }
        if (str == "cells")      { return fs_type::cells; }
        if (str == "nodes")      { return fs_type::nodes; }
        if (str == "points")     { return fs_type::points; }
        return fs_type::unset;
    };
    FunctionSpace sfs, tfs;
    std::string interpolation_method = get_interpolation_method(args);
    fs_type sfs_type = fs_type::unset;
    fs_type tfs_type = fs_type::unset;
    bool target_follows_source = true;
    int shalo = 0;
    int thalo = 0;
    if (string_starts_with(interpolation_method, "cons")) {
        int order = string_starts_with(interpolation_method, "cons2") ? 2 : 1;
        target_follows_source = false;
        sfs_type = gridpoints_are_cells(sgrid) ? fs_type::cells : fs_type::nodes;
        tfs_type = gridpoints_are_cells(tgrid) ? fs_type::cells : fs_type::nodes;
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
            ATLAS_THROW_EXCEPTION("Interpolation method " << interpolation_method << " requires that the source grid is a StructuredGrid");
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
            ATLAS_THROW_EXCEPTION("Interpolation method " << interpolation_method << " requires that the source grid is a StructuredGrid");
        }
        if (not StructuredGrid(tgrid)) {
            ATLAS_THROW_EXCEPTION("Interpolation method " << interpolation_method << " requires that the target grid is a StructuredGrid");
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
            Log::warning() << "WARNING: The source will be meshed for reason of creating the distribution of the target grid" << std::endl;
            sfs_type = fs_type::nodes;
        }
    }
    args.get("s.halo", shalo);
    args.get("t.halo", thalo);

    if (target_follows_source && sgrid.type() == "ORCA" && mpi::size() > 1) {
        ATLAS_THROW_EXCEPTION("ORCA grid is not yet supported as source in parallel for interpolation method \"" << interpolation_method << "\"");
    }
    if (not target_follows_source && tgrid.type() == "ORCA" && mpi::size() > 1) {
        ATLAS_THROW_EXCEPTION("ORCA grid is not yet supported as target in parallel for interpolation method \"" << interpolation_method << "\"");
    }
    if (string_starts_with(interpolation_method,"grid-box") && shalo > 0) {
        ATLAS_THROW_EXCEPTION( "Interpolation method \"" << interpolation_method << "\" only supported with --s.halo=0 (default)");
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
            fs = CellColumns(mesh, option::halo(halo));
        }
        else if (type == fs_type::nodes) {
            auto mesh = create_mesh_with_gridpoints_as_nodes(grid, dist, halo);
            fs = NodeColumns(mesh, option::halo(halo));
            auto mask = mesh.nodes().add(fs.createField<int>(option::name("mask")));
            array::make_view<int,1>(mask).assign(1);
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
    return std::make_pair(sfs,tfs);
}


template <typename View>
void gmsh_output(const std::string& name, const Grid& grid, const View& field, const AtlasTool::Args& args) {
    mpi::Scope mpi_scope("self");
    Log::info() << "Storing field '" << name << ".msh'" << std::endl;
    ATLAS_ASSERT(field.size() <= grid.size());
    std::string coords = args.getString("gmsh.coordinates", "lonlat");
    output::Gmsh gmsh(name + ".msh", Config("coordinates", coords) | Config("ghost", "true"));
    auto mesh = Mesh(grid);
    gmsh.write(mesh);
    auto fs = gridpoints_are_cells(grid) ? FunctionSpace(CellColumns(mesh)) : FunctionSpace(NodeColumns(mesh));
    auto f = fs.createField<double>(option::name(name));
    if (args.getBool("gmsh.missing-value",false)) {
        f.metadata().set("missing_value",args.getDouble("missing-value",0.));
        f.metadata().set("missing_value_type","equals");
    }
    auto v = array::make_view<double,1>(f);
    auto g = array::make_view<gidx_t,1>(fs.global_index());
    size_t min_size = std::min<size_t>(f.size(),field.size());
    for(size_t i=0; i<min_size; ++i) {
        v(i) = field[g(i)-1];
    }
    fs.haloExchange(f);
    gmsh.write(f);
}

void test_matrix(const Grid& sgrid, const Grid& tgrid, const Matrix& matrix, const AtlasTool::Args& args) {
    ATLAS_TRACE();

    // Allocate and initialise own memory here to show possibilities
    std::vector<double> sdata(sgrid.size());
    std::vector<double> tdata(tgrid.size());

    ATLAS_TRACE_SCOPE("initialize source") {
        if (args.has("init") || not args.has("s.grib")) {
            ATLAS_DEBUG("using init");
            auto init = get_init(args);
            idx_t n{0};
            for (auto p : sgrid.lonlat()) {
                sdata[n++] = init(p);
            }
        }
        else {
            ATLAS_DEBUG("using grib");

            ATLAS_ASSERT(args.has("s.grib"));
#if ATLAS_HAVE_GRIB
            GribFileReader grib_reader(args.getString("s.grib"));
            grib_reader.grib().get_values(sdata.data(),sdata.size());
#endif
        }
    }

    StopWatch timer_serial_sparse_matrix_multiply;
    timer_serial_sparse_matrix_multiply.start();
    atlas::linalg::sparse_matrix_multiply(
        atlas::linalg::make_host_view<double>(matrix),
        atlas::array::make_view<double,1>(sdata.data(), matrix.cols()),
        atlas::array::make_view<double,1>(tdata.data(), matrix.rows()));
    timer_serial_sparse_matrix_multiply.stop();
    Log::info() << "Serial sparse-matrix-multiply timer  \t: " << elapsed_ms(timer_serial_sparse_matrix_multiply.elapsed(),true) << " [ms]" << std::endl;
    Log::info() << "Serial sparse-matrix non-zero entries\t: " << matrix.nnz() << std::endl;

    {
        double missing_value = args.getDouble("missing-value",0.);
        auto m = atlas::linalg::make_host_view<double>(matrix);
        for (idx_t r=0; r<m.rows(); ++r) {
            if (m.outer()[r+1] - m.outer()[r] == 0) {
                tdata[r] = missing_value;
            }
        }
    }

    if (args.getBool("output-checksum",false)) {
        // The matrix version is possibly not bit-identical with the interpolator version
        // This is because the matrix version uses grid points without halo, wrapping around the globe
        size_t min_size = std::min(tgt_ref.size(), tdata.size());
        for (size_t i=0; i<min_size; ++i) {
            if (std::abs(tgt_ref[i] - tdata[i]) < 1.e-14) {
                tdata[i] = tgt_ref[i];
            }
        }

        std::ofstream checksum_file(args.getString("checksum.file","checksum"),std::ios::app);
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

    // atlas::array::make_view<double,1>(tdata.data(), matrix.rows()).dump(Log::info());
    // Log::info() << std::endl;


    if (args.getBool("output-gmsh",false)) {
        ATLAS_TRACE_SCOPE("Gmsh serial output") {

            if (args.has("s.mask")) {
                std::vector<int> mask(sgrid.size());
                {
                    std::string mask_name = args.getString("s.mask");
                    if (get_mask_format(mask_name) == "scrip") {
                        ScripIO::read_mask(mask_name, mdspan<int,dims<1>>(mask.data(), mask.size()));
                    }
                    else if (get_mask_format(mask_name) == "atlas") {
                        AtlasIO::read_mask(mask_name, mdspan<int,dims<1>>(mask.data(), mask.size()));
                    }
                    else {
                        ATLAS_NOTIMPLEMENTED;
                    }
                }
                double missing_value = args.getDouble("missing-value",0.);
                for(size_t i=0; i<sdata.size(); ++i) {
                    if(mask[i] == 0) {
                        sdata[i] = missing_value;
                    }
                }
            }

            std::string sname = "test_matrix_source_" + get_matrix_name(args);
            gmsh_output(sname, sgrid, atlas::array::make_view<double,1>(sdata.data(), matrix.cols()), args);

            if (args.has("t.mask")) {
                std::vector<int> mask(tgrid.size());
                {
                    std::string mask_name = args.getString("t.mask");
                    if (get_mask_format(mask_name) == "scrip") {
                        ScripIO::read_mask(mask_name, mdspan<int,dims<1>>(mask.data(), mask.size()));
                    }
                    else if (get_mask_format(mask_name) == "atlas") {
                        AtlasIO::read_mask(mask_name, mdspan<int,dims<1>>(mask.data(), mask.size()));
                    }
                    else {
                        ATLAS_NOTIMPLEMENTED;
                    }
                }
                double missing_value = args.getDouble("missing-value",0.);
                for(size_t i=0; i<tdata.size(); ++i) {
                    if(mask[i] == 0) {
                        tdata[i] = missing_value;
                    }
                }
            }
            std::string tname = "test_matrix_target_" + get_matrix_name(args);
            gmsh_output(tname, tgrid, atlas::array::make_view<double,1>(tdata.data(), matrix.rows()), args);
        }
    }
}

int AtlasInterpolations::execute(const AtlasTool::Args& args) {
    ATLAS_TRACE("main");
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

    bool output_checksum = args.getBool("output-checksum",false);
    std::string checksum_file_path = args.getString("checksum.file","checksum");
    if (output_checksum) {
        eckit::PathName checksum_file(checksum_file_path);
        if (checksum_file.exists()) {
            checksum_file.unlink();
        }
    }

    auto sgrid = get_sgrid(args);
    auto tgrid = get_tgrid(args);
    std::string sdata_name = gridpoints_are_cells(sgrid) ? "cell data" : "node data";
    std::string tdata_name = gridpoints_are_cells(tgrid) ? "cell data" : "node data";

    std::string interpolation_method = get_interpolation_method(args);
    Log::info() << "\tsource grid\t:\t" << sgrid.name() << " (" << sdata_name <<")\n\ttarget grid\t:\t" << tgrid.name()
        << " (" << tdata_name <<")\n\tinterpolation\t:\t" << interpolation_method << std::endl;

    std::string matrix_format = get_matrix_format(args);
    std::string matrix_name = get_matrix_name(args);
    bool matrix_tested = false;
    Matrix matrix;
    if (args.getBool("read-matrix",false)) {
        if (mpi::comm().rank() == 0) {
            ATLAS_TRACE_SCOPE("Read matrix on rank 0") {
                timers.global_matrix_read.start();
                matrix = read_matrix(matrix_name, matrix_format);
                timers.global_matrix_read.stop();
                Log::info() << "Global matrix read timer\t: " << elapsed_ms(timers.global_matrix_read,true) << " [ms]"  << std::endl;
                print_matrix(matrix, Log::info());
                Log::info() << std::endl;
                if (output_checksum) {
                    std::ofstream checksum_file(checksum_file_path ,std::ios::app);
                    auto matrix_view = linalg::make_host_view<double,int>(matrix);
                    auto outer_checksum = util::checksum(matrix_view.inner(), matrix_view.inner_size());
                    auto inner_checksum = util::checksum(matrix_view.outer(), matrix_view.outer_size());
                    auto value_checksum = util::checksum(matrix_view.value(), matrix_view.value_size());
                    checksum_file << std::setw(8) << outer_checksum << "    [read-matrix]   checksum of matrix.outer()" << std::endl;
                    checksum_file << std::setw(8) << inner_checksum << "    [read-matrix]   checksum of matrix.inner()" << std::endl;
                    checksum_file << std::setw(8) << value_checksum << "    [read-matrix]   checksum of matrix.value()" << std::endl;
                }
            }
        }
        if (args.getBool("test-matrix",false)) {
            if (mpi::comm().rank() == 0) {
                test_matrix(sgrid, tgrid, matrix, args);
            }
            matrix_tested = true;
        }
    }

    FunctionSpace src_fs;
    FunctionSpace tgt_fs;
    Interpolation interpolator;

    ATLAS_TRACE_SCOPE("Setup source and target") {
        timers.functionspace_setup.start();
        // create_fspaces(config, sgrid, tgrid, src_fs, tgt_fs);
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
        Log::info() << "Global matrix setup timer\t: " << elapsed_ms(timers.global_matrix_setup) << " [ms]"  << std::endl;
    }

    matrix.clear();

    if (args.getBool("interpolate",false) || args.getBool("output-matrix",false) || args.getBool("test-interpolator",false) || (args.getBool("test-matrix",false) && !matrix_tested)) {
;
        if (args.has("s.mask")) {
            Field smask_glb = src_fs.createField<int>(option::name("mask")|option::global());
            {
                std::string mask = args.getString("s.mask");
                if (get_mask_format(mask) == "scrip") {
                    ScripIO::read_mask(mask, array::make_view<int,1>(smask_glb).as_mdspan());
                }
                else if (get_mask_format(mask) == "atlas") {
                    AtlasIO::read_mask(mask, array::make_view<int,1>(smask_glb).as_mdspan());
                }
                else {
                    ATLAS_NOTIMPLEMENTED;
                }
            }
            Field smask = src_fs.mask();
            src_fs.scatter(smask_glb, smask);
            src_fs.haloExchange(smask);

            // smask = src_fs.createField<int>(option::name("smask")|option::global());
            //gmsh_output("smask", sgrid, array::make_view<int,1>(smask), args);
        }

        if (args.has("t.mask")) {
            Field tmask_glb = tgt_fs.createField<int>(option::name("mask")|option::global());
            {
                std::string mask = args.getString("t.mask");
                if (get_mask_format(mask) == "scrip") {
                    ScripIO::read_mask(mask, array::make_view<int,1>(tmask_glb).as_mdspan());
                }
                else if (get_mask_format(mask) == "atlas") {
                    AtlasIO::read_mask(mask, array::make_view<int,1>(tmask_glb).as_mdspan());
                }
                else {
                    ATLAS_NOTIMPLEMENTED;
                }
            }
            Field tmask = tgt_fs.mask();
            tgt_fs.scatter(tmask_glb, tmask);
            tgt_fs.haloExchange(tmask);
        }

        ATLAS_TRACE_SCOPE("Setup interpolator") {
            timers.interpolation_setup.start();
            auto config = get_interpolation_config(sgrid, tgrid, args);
            if (string_starts_with(interpolation_method,"grid-box")) {
                interpolator = Interpolation(config, sgrid, tgrid);
            }
            else {
                interpolator = Interpolation(config, src_fs, tgt_fs, cache);
            }
            timers.interpolation_setup.stop();
            ATLAS_DEBUG_VAR(interpolator.failedInterpolations().size());
        }

        Log::info() << "Grid + FunctionSpace timer\t: " << elapsed_ms(timers.functionspace_setup) << " [ms]"  << std::endl;
        Log::info() << "Interpolation setup timer\t: " << elapsed_ms(timers.interpolation_setup) << " [ms]"  << std::endl;

        if (args.getBool("test-interpolator", false)) {
            auto src_field = interpolator.source().createField<double>();
            auto tgt_field = interpolator.target().createField<double>();
            auto src_lonlat = array::make_view<double, 2>(interpolator.source().lonlat());
            ATLAS_TRACE_SCOPE("initialize source") {
                auto init = get_init(args);
                auto src_field_v = array::make_view<double, 1>(src_field);
                for (idx_t i = 0; i < src_fs.size(); ++i) {
                    src_field_v[i] = init(PointLonLat{src_lonlat(i, 0), src_lonlat(i, 1)});
                }
            }
            src_field.set_dirty(false); // all values are up to date

            timers.interpolation_exe.start();
            interpolator.execute(src_field, tgt_field);
            timers.interpolation_exe.stop();
            Log::info() << "Interpolation execute timer     : " << elapsed_ms(timers.interpolation_exe) << " [ms]"  << std::endl;

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
                tgt_field_global = tgt_fs.createField(tgt_field, option::global());
                auto tgt_field_global_v = array::make_view<double,1>(tgt_field_global);
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
                    checksum_file << std::setw(8) << target_checksum << "    [test-interp]   checksum of target field" << std::endl;

                    // for further comparison with test-matrix
                    tgt_ref.resize(tgt_field_global_v.size());
                    for( idx_t i=0; i<tgt_field_global_v.size(); ++i) {
                        tgt_ref[i] = tgt_field_global_v[i];
                    }
                }
            }

            if (args.getBool("output-gmsh",false)) {
                Mesh tmesh;
                if (auto fs = NodeColumns(tgt_fs)) { tmesh = fs.mesh(); }
                if (auto fs = CellColumns(tgt_fs)) { tmesh = fs.mesh(); }
                if (not tmesh) {
                    Log::warning() << "Target function space is NOT mesh based, hence no Gmsh output will be created for this" << std::endl;
                }
                else {
                    tgt_field.haloExchange();
                    std::string tgt_name = "test_interpolator_target_" + matrix_name;
                    Log::info() << "Storing field '" << tgt_name << ".msh'." << std::endl;
                    std::string coords = args.getString("gmsh.coordinates", "lonlat");
                    output::Gmsh gmsh(tgt_name + ".msh", Config("coordinates", coords) | Config("ghost", "false"));
                    gmsh.write(tmesh);
#if 0
                    auto tview = array::make_view<double,1>(tgt_field);
                    for(size_t i=0; i<tview.shape(0); ++i) {
                        if (tview(i) == 9999.) {
                            tview(i) = 0;
                        }
                    }
#endif
                    gmsh.write(tgt_field);
                }
            }
#if 0
            // This writes the mask to file, can be used to interpolate masks
            std::vector<int> output_mask(tgrid.size());
            auto tview = array::make_view<double,1>(tgt_field);
            for(size_t i=0; i<output_mask.size(); ++i) {
                output_mask[i] = (tview(i) != 0 && tview(i) != 9999.) ? 1 : 0;
            }
            AtlasIO::write_mask("mask_"+tgrid.name()+".atlas", tgrid.name(), mdspan<int,dims<1>>(output_mask.data(), output_mask.size()));
#endif
        }
    }

    if (args.getBool("output-matrix",false) || (args.getBool("test-matrix",false) && !matrix_tested)) {
        matrix = interpolation::assemble_global_matrix(tgrid.size(), sgrid.size(), interpolator);
        if (mpi::comm().rank() == 0) {
            if (args.getBool("output-matrix",false)) {
                write_matrix(matrix, matrix_name, matrix_format);
                print_matrix(matrix, Log::info());
                if (output_checksum) {
                    std::ofstream checksum_file(checksum_file_path ,std::ios::app);
                    auto matrix_view = linalg::make_host_view<double,int>(matrix);
                    auto outer_checksum = util::checksum(matrix_view.inner(), matrix_view.inner_size());
                    auto inner_checksum = util::checksum(matrix_view.outer(), matrix_view.outer_size());
                    auto value_checksum = util::checksum(matrix_view.value(), matrix_view.value_size());
                    checksum_file << std::setw(8) << outer_checksum << "    [output-matrix] checksum of matrix.outer()" << std::endl;
                    checksum_file << std::setw(8) << inner_checksum << "    [output-matrix] checksum of matrix.inner()" << std::endl;
                    checksum_file << std::setw(8) << value_checksum << "    [output-matrix] checksum of matrix.value()" << std::endl;
                }

            }
            if (args.getBool("test-matrix",false) && !matrix_tested) {
                test_matrix(sgrid, tgrid, matrix, args);
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
        eckit_matrix.save(matrix_name+".eckit");
    }
    else if (format == "scrip") {
        Log::info() << "Writing matrix in scrip format to file '" << matrix_name << ".nc'" << std::endl;
        ScripIO::write_matrix(matrix, matrix_name+".nc");
    }
    else {
        ATLAS_NOTIMPLEMENTED;
    }
}

void AtlasInterpolations::print_matrix(const Matrix& matrix, std::ostream& out) {
    out << "    matrix.rows : " << matrix.rows() << '\n';
    out << "    matrix.cols : " << matrix.cols() << '\n';
    out << "    matrix.nnz  : " << matrix.nnz()  << '\n';
}


} // namespace


int main(int argc, char** argv) {
    atlas::AtlasInterpolations tool(argc, argv);
    return tool.start();
}
