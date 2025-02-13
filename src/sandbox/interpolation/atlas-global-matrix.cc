/*
 * (C) Copyright 2025 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "ScripIO.h"

#include "atlas/array.h"
#include "atlas/field/Field.h"
#include "atlas/field/FieldSet.h"
#include "atlas/functionspace/CellColumns.h"
#include "atlas/functionspace/NodeColumns.h"
#include "atlas/functionspace/PointCloud.h"
#include "atlas/functionspace/StructuredColumns.h"
#include "atlas/grid/Grid.h"
#include "atlas/grid/Iterator.h"
#include "atlas/interpolation.h"
#include "atlas/interpolation/AssembleGlobalMatrix.h"
#include "atlas/linalg/sparse.h"
#include "atlas/linalg/sparse/MakeEckitSparseMatrix.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/meshgenerator.h"
#include "atlas/output/Gmsh.h"
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


// USAGE:
//
// 1) To compute and store the global interpolation matrix to the disc:
//
// atlas-global-matrix --sgrid <sgrid> --tgrid <tgrid> --interpolation <intrp> --format <fmt>
//
// 2) To read in a stored global interpolation matrix from the disc:
//
// atlas-global-matrix --sgrid <sgrid> --tgrid <tgrid> --interpolation <intrp> --format <fmt> --read --matrix <matrix>
//
// For available interpolations see AtlasGlobalMatrix::interpolation_config.
// Formats <fmt> can be 'eckit' (binary) or 'scrip' (netcdf).
//
// Using the grid API 'create_fspaces' we can hide interpolation method specific requirements
// such as which functionspace needs to be set-up.


class AtlasGlobalMatrix : public AtlasTool {
    int execute(const AtlasTool::Args& args) override;
    std::string briefDescription() override { return "Assemble global matrix from an Interpolation in distributed parallel run"; }
    std::string usage() override {
        return name() + " [OPTION] ... [--help]";
    }

    int numberOfPositionalArguments() override { return -1; }
    int minimumPositionalArguments() override { return 0; }
    Config interpolation_config(std::string scheme_str);
    Config create_fspaces(const std::string& scheme_str, const Grid& input_grid, const Grid& output_grid,
        FunctionSpace& fs_in, FunctionSpace& fs_out);
    std::string get_matrix_name(std::string sgrid, std::string tgrid, std::string interp = "");
    Matrix read_matrix(std::string matrix_name, std::string matrix_format);
    void write_matrix(const Matrix& mat, std::string matrix_name, std::string matrix_format);

public:
    AtlasGlobalMatrix(int argc, char* argv[]): AtlasTool(argc, argv) {
        add_option(new SimpleOption<std::string>("sgrid", "source grid"));
        add_option(new SimpleOption<std::string>("tgrid", "target grid"));
        add_option(new SimpleOption<std::string>("interpolation", "interpolation methods: linear, cubic, qcubic, cons, cons2, nneighbour, knneighbour, grid-box"));
        add_option(new SimpleOption<bool>("read", "read interpolation matrix"));
        add_option(new SimpleOption<std::string>("matrix", "name of the interpolation matrix"));
        add_option(new SimpleOption<std::string>("format", "format of the matrix output: eckit, SCRIP"));
    }

    struct Timers {
        using StopWatch = atlas::runtime::trace::StopWatch;
        StopWatch fspace_setup;
        StopWatch interpolation_setup;
        StopWatch global_matrix;
    } timers;
};

int AtlasGlobalMatrix::execute(const AtlasTool::Args& args) {
    std::string option;
    std::string sgrid_name = "O8";
    args.get("sgrid", sgrid_name);
    std::string tgrid_name = "O32";
    args.get("tgrid", tgrid_name);

    auto sgrid = Grid{sgrid_name};
    auto tgrid = Grid{tgrid_name};
    Log::info() << "source grid: " << sgrid_name << ", target grid: " << tgrid_name;

    std::string matrix_format = "eckit";
    args.get("format", matrix_format);

    if (args.has("read")) {
        Matrix gmatrix;
        if (mpi::comm().rank() == 0) {
            std::string matrix_name;
            args.get("matrix", matrix_name);
            if (matrix_name == "") {
                matrix_name = get_matrix_name(sgrid_name, tgrid_name);
            }
            Log::info() << "\nreading matrix '" << matrix_name << "' from the disc.\n";

            gmatrix = read_matrix(matrix_name, matrix_format);
            auto eckit_gmatrix = atlas::linalg::make_non_owning_eckit_sparse_matrix(gmatrix);
            eckit_gmatrix.print(Log::info());
            Log::info() << std::endl;
        }

        mpi::comm().barrier();

        std::string interpolation_method = "finite-element";
        args.get("interpolation", interpolation_method);
        FunctionSpace src_fs;
        FunctionSpace tgt_fs;
        auto scheme = create_fspaces(interpolation_method, sgrid, tgrid, src_fs, tgt_fs);
        auto matrix = interpolation::distribute_global_matrix(src_fs, tgt_fs, gmatrix);

        interpolation::MatrixCache cache(std::move(matrix));
        auto interpolator = Interpolation(scheme, src_fs, tgt_fs, cache);

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

        interpolator.execute(src_field, tgt_field);

        ATLAS_TRACE_SCOPE("output from the read-in matrix") {
            std::string tgt_name = "tfield_" + get_matrix_name(sgrid_name, tgrid_name);
            output::Gmsh gmsh(tgt_name + ".msh", Config("coordinates", "lonlat") | Config("ghost", "true"));
            if( functionspace::NodeColumns(tgt_field.functionspace())) {
                Log::info() << "storing distributed remapped field '" << tgt_name << "'." << std::endl;
                gmsh.write(functionspace::NodeColumns(tgt_field.functionspace()).mesh());
                tgt_field.haloExchange();
                gmsh.write(tgt_field);
            }
        }
    }
    else {
        std::string interpolation_method = "finite-element";
        args.get("interpolation", interpolation_method);
        Log::info() << ", interpolation: " << interpolation_method << std::endl;

        FunctionSpace src_fs;
        FunctionSpace tgt_fs;
        Interpolation interpolator;
        timers.fspace_setup.start();
        Config scheme = create_fspaces(interpolation_method, sgrid, tgrid, src_fs, tgt_fs);
        timers.fspace_setup.stop();

        timers.interpolation_setup.start();
        if (interpolation_method == "nearest-neighbour" || interpolation_method == "k-nearest-neighbours" || interpolation_method == "grid-box-average") {
            interpolator = Interpolation{ scheme, sgrid, tgrid };
        }
        else {
            interpolator = Interpolation{ scheme, src_fs, tgt_fs };
        }
        timers.interpolation_setup.stop();

        double timer_fspace_setup = timers.fspace_setup.elapsed() / mpi::comm().size();
        double timer_interp_setup = timers.interpolation_setup.elapsed() / mpi::comm().size();
        mpi::comm().allReduceInPlace(&timer_fspace_setup, 1, eckit::mpi::sum());
        mpi::comm().allReduceInPlace(&timer_interp_setup, 1, eckit::mpi::sum());
        Log::info() << "Grid + FunctionSpace timer      : " << timer_fspace_setup * 1000. << " [ms]"  << std::endl;
        Log::info() << "Interpolation setup timer       : " << timer_interp_setup * 1000. << " [ms]"  << std::endl;

        ATLAS_TRACE_SCOPE("par_output") {
            std::string tgt_name = "par-" + scheme.getString("name");
            auto tgt_field = interpolator.target().createField<double>();
            auto field_in = interpolator.source().createField<double>();
            auto lonlat_in = array::make_view<double,2>(interpolator.source().lonlat());
            auto view_in = array::make_view<double,1>(field_in);
            for(idx_t j = 0; j < field_in.size(); ++j) {
                view_in(j) = util::function::vortex_rollup(lonlat_in(j,0), lonlat_in(j,1), 1.);
            }
            interpolator.execute(field_in, tgt_field);
            output::Gmsh gmsh(tgt_name + ".msh", Config("coordinates", "lonlat") | Config("ghost", "true"));
            if( functionspace::NodeColumns(tgt_field.functionspace())) {
                Log::info() << "storing distributed remapped field '" << tgt_name << "'." << std::endl;
                gmsh.write(functionspace::NodeColumns(tgt_field.functionspace()).mesh());
                tgt_field.haloExchange();
                gmsh.write(tgt_field);
            }
        }

        auto matrix = interpolation::assemble_global_matrix(interpolator);

        if (mpi::comm().rank() == 0) {

            ATLAS_TRACE_SCOPE("store matrix") {
                std::string matrix_name;
                args.get("matrix", matrix_name);
                if (matrix_name == "") {
                    matrix_name = get_matrix_name(sgrid_name, tgrid_name, scheme.getString("name"));
                }
                write_matrix(matrix, matrix_name, matrix_format);
            }

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
            timers.global_matrix.start();
            eckit::linalg::LinearAlgebraSparse::backend().spmv(eckit_matrix, src, tgt);
            timers.global_matrix.stop();
            Log::info() << "Global matrix-multiply timer    : " << 1000.*timers.global_matrix.elapsed()  << " [ms]" << std::endl;
            Log::info() << "Global matrix non-zero entries  : " << matrix.nnz() << std::endl;
            Log::info() << "Global matrix memory            : " << eckit_matrix.footprint() << std::endl;

            ATLAS_TRACE_SCOPE("output from proc 0") {
                mpi::Scope mpi_scope("self"); 
                std::string tgt_name = "tfield_cache_" + get_matrix_name(sgrid_name, tgrid_name);
                auto tgt_field = Field(tgt_name, tgt_data.data(), array::make_shape(tgt_data.size())); // wrap
                output::Gmsh gmsh(tgt_name + ".msh", Config("coordinates", "lonlat") | Config("ghost", "true"));
                gmsh.write(MeshGenerator("structured", util::Config("three_dimensional",true)).generate(tgrid));
                Log::info() << "storing remapped field '" << tgt_field.name() << "'." << std::endl;
                gmsh.write(tgt_field, StructuredColumns(tgrid));
            }
        }
    }

    return success();
}


//-----------------------------------------------------------------------------


Config AtlasGlobalMatrix::interpolation_config(std::string scheme_str) {
    Config scheme;
    scheme.set("type", scheme_str);
    scheme.set("halo", 1);
    if (scheme_str.find("cubic") != std::string::npos) {
        scheme.set("halo", 2);
    }
    if (scheme_str == "k-nearest-neighbours") {
        scheme.set("k-nearest-neighbours", 4);
        scheme.set("halo", 2);
    }
    if (scheme_str == "conservative-spherical-polygon") {
        scheme.set("type", "conservative-spherical-polygon");
        scheme.set("order", 1);
        scheme.set("src_cell_data", false);
        scheme.set("tgt_cell_data", false);
    }
    if (scheme_str == "conservative-spherical-polygon-2") {
        scheme.set("type", "conservative-spherical-polygon");
        scheme.set("order", 2);
        scheme.set("halo", 2);
        scheme.set("src_cell_data", false);
        scheme.set("tgt_cell_data", false);
    }
    scheme.set("name", scheme_str);
    return scheme;
}


std::string AtlasGlobalMatrix::get_matrix_name(std::string sgrid, std::string tgrid, std::string interp) {
    std::stringstream ss;
    ss << mpi::comm().size() << "-" << atlas_omp_get_num_threads();
    std::string grids = "remap_" + ss.str() + "_" + sgrid + "_" + tgrid;
    return (interp != "") ? grids + "_" + interp : grids;
}


Config AtlasGlobalMatrix::create_fspaces(const std::string& scheme_str, const Grid& input_grid, const Grid& output_grid,
        FunctionSpace& fs_in, FunctionSpace& fs_out) {
    Config scheme = interpolation_config(scheme_str);
    auto scheme_type = scheme.getString("type");
    if (scheme_type == "finite-element" || scheme_type == "unstructured-bilinear-lonlat") {
        auto inmesh = Mesh(input_grid);
        fs_in = NodeColumns(inmesh, scheme);
        fs_out = PointCloud(output_grid, grid::MatchingPartitioner(inmesh));
    }
    else if (scheme_type == "conservative-spherical-polygon" || scheme_type == "grid-box-average") {
        bool src_cell_data = scheme.getBool("src_cell_data");
        bool tgt_cell_data = scheme.getBool("tgt_cell_data");
        auto tgt_mesh_config = output_grid.meshgenerator() | option::halo(0);
        auto tgt_mesh = MeshGenerator(tgt_mesh_config).generate(output_grid);
        if (tgt_cell_data) {
            fs_out = functionspace::CellColumns(tgt_mesh, option::halo(0));
        }
        else {
            fs_out = functionspace::NodeColumns(tgt_mesh, option::halo(1));
        }
        auto src_mesh_config = input_grid.meshgenerator() | option::halo(2);
        Mesh src_mesh;
        if (mpi::size() > 1) {
            src_mesh = MeshGenerator(src_mesh_config).generate(input_grid, grid::MatchingPartitioner(tgt_mesh));
        }
        else {
            src_mesh = MeshGenerator(src_mesh_config).generate(input_grid);
        }
        if (src_cell_data) {
            fs_in = functionspace::CellColumns(src_mesh, option::halo(2));
        }
        else {
            fs_in = functionspace::NodeColumns(src_mesh, option::halo(2));
        }
    }
    else if (scheme_type == "nearest-neighbour" || scheme_type == "k-nearest-neighbours") {
        fs_in = PointCloud(input_grid);
        fs_out = functionspace::PointCloud(output_grid, grid::MatchingPartitioner(fs_in));
    }
    else {
        fs_in = StructuredColumns(input_grid, scheme);
        fs_out = functionspace::PointCloud(output_grid, grid::MatchingPartitioner(fs_in), scheme);
    }
    return scheme;
}


Matrix AtlasGlobalMatrix::read_matrix(std::string matrix_name, std::string format) {
    if (format == "eckit") {
        eckit::linalg::SparseMatrix eckit_matrix;
        Log::info() << "reading matrix '" << matrix_name << std::endl;
        eckit_matrix.load(matrix_name);
        return atlas::linalg::make_sparse_matrix_storage(std::move(eckit_matrix));
    }
    else if (format == "scrip") {
        Log::info() << "reading matrix '" << matrix_name << std::endl;
        return ScripIO::read(matrix_name);
    }
    else {
        ATLAS_NOTIMPLEMENTED;
    }
    return Matrix{};
}


void AtlasGlobalMatrix::write_matrix(const Matrix& matrix, std::string matrix_name, std::string format) {
    if (format == "eckit") {
        auto eckit_matrix = linalg::make_non_owning_eckit_sparse_matrix(matrix);
        Log::info() << "storing matrix '" << matrix_name << ".eckit'" << std::endl;
        eckit_matrix.save(matrix_name + ".eckit");
    }
    else if (format == "scrip") {
        Log::info() << "storing matrix '" << matrix_name << ".nc'" << std::endl;
        ScripIO::write(matrix, matrix_name + ".nc");
    }
    else {
        ATLAS_NOTIMPLEMENTED;
    }
}

} // namespace


int main(int argc, char** argv) {
    atlas::AtlasGlobalMatrix tool(argc, argv);
    return tool.start();
}
