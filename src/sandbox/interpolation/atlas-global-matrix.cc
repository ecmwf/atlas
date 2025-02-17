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
// atlas-global-matrix --sgrid <sgrid> --tgrid <tgrid> --interpolation <intrp> --format <fmt>
//
// For available interpolations see AtlasGlobalMatrix::interpoaltion_config.
// Formats can be 'eckit' (binary) or 'scrip' (netcdf).
//
// Using the grid API we can hide interpolation method specific requirements
// such as which functionspace needs to be set-up.


class AtlasGlobalMatrix : public AtlasTool {
    int execute(const AtlasTool::Args& args) override;
    std::string briefDescription() override { return "Assemble global matrix from an Interpolation in distributed parallel run"; }
    std::string usage() override {
        return name() + " [OPTION]... [--help]";
    }

    int numberOfPositionalArguments() override { return -1; }
    int minimumPositionalArguments() override { return 0; }
    Config interpolation_config(std::string scheme_str);
    std::string get_matrix_name(std::string sgrid, std::string tgrid, std::string interp = "");
    Matrix read_matrix(std::string matrix_name, std::string matrix_format);
    void write_matrix(const Matrix& mat, std::string matrix_name, std::string matrix_format);

public:
    AtlasGlobalMatrix(int argc, char* argv[]): AtlasTool(argc, argv) {
        add_option(new SimpleOption<std::string>("sgrid", "source grid"));
        add_option(new SimpleOption<std::string>("tgrid", "target grid"));
        add_option(new SimpleOption<std::string>("interpolation", "interpolation method (default finite-element): linear, cubic, qcubic, cons, cons2, nneighbour, knneighbour, grid-box"));
        add_option(new SimpleOption<bool>("read", "read interpolation matrix (default <sgrid>_<tgrid>.matrix)"));
        add_option(new SimpleOption<std::string>("matrix", "name of interpolation matrix: (default <sgrid>_<tgrid>.matrix)"));
        add_option(new SimpleOption<std::string>("format", "format of the matrix output (eckit, scrip) (default eckit)"));
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
        if (mpi::comm().size() > 1) {
            Log::info() << "\n\n Reading of global matrix has to be done with only one task. Please rerun in serial !\n";
            return failed();
        }
        std::string matrix_name;
        args.get("matrix", matrix_name);
        if (matrix_name == "") {
            matrix_name = get_matrix_name(sgrid_name, tgrid_name);
        }
        Log::info() << "\nreading matrix '" << matrix_name << "' from the disc.\n";

        Matrix matrix = read_matrix(matrix_name, matrix_format);
        auto eckit_matrix = atlas::linalg::make_non_owning_eckit_sparse_matrix(matrix);
        eckit_matrix.print(Log::info());
        Log::info() << std::endl;

        // Allocate and initialise own memory here to show possibilities
        // Note: reading a field from disc is an extra feature
        std::vector<double> src_data(sgrid.size());
        std::vector<double> tgt_data(tgrid.size());
        auto src = eckit::linalg::Vector(src_data.data(), src_data.size());
        auto tgt = eckit::linalg::Vector(tgt_data.data(), tgt_data.size());

        ATLAS_TRACE_SCOPE("initialize source") {
            idx_t n{0};
            for (auto p : sgrid.lonlat()) {
                src_data[n++] = util::function::vortex_rollup(p.lon(), p.lat(), 1.);
            }
        }
        eckit::linalg::LinearAlgebraSparse::backend().spmv(eckit_matrix, src, tgt);

        ATLAS_TRACE_SCOPE("output from read matrix") {
            mpi::Scope mpi_scope("self"); 
            std::string tgt_name = "tfield_" + get_matrix_name(sgrid_name, tgrid_name);
            auto tgt_field = Field(tgt_name, tgt_data.data(), array::make_shape(tgt_data.size())); // wrap
            output::Gmsh gmsh(tgt_name + ".msh", Config("coordinates", "lonlat") | Config("ghost", "true"));
            gmsh.write(MeshGenerator("structured", util::Config("three_dimensional",true)).generate(tgrid));
            Log::info() << "storing remapped field '" << tgt_name << "' to the disc." << std::endl;
            gmsh.write(tgt_field, StructuredColumns(tgrid));
        }
    }
    else {
        std::string interpolation_method = "finite-element";
        args.get("interpolation", interpolation_method);
        Log::info() << ", interpolation: " << interpolation_method << std::endl;

        Interpolation interpolator;
        Config scheme = interpolation_config(interpolation_method);
        auto scheme_type = scheme.getString("type");
        if (scheme_type == "finite-element") {
            timers.fspace_setup.start();
            auto inmesh = Mesh(sgrid);
            auto fs_in = NodeColumns(inmesh, scheme);
            auto fs_out = PointCloud(tgrid, grid::MatchingPartitioner(inmesh));
            timers.fspace_setup.stop();

            timers.interpolation_setup.start();
            interpolator = Interpolation{ scheme, fs_in, fs_out };
            timers.interpolation_setup.stop();
        }
        else if (scheme_type == "conservative-spherical-polygon" || scheme_type == "grid-box-average") {
            timers.fspace_setup.start();
            timers.fspace_setup.stop();

            timers.interpolation_setup.start();
            interpolator = Interpolation{ scheme, sgrid, tgrid };
            timers.interpolation_setup.stop();
        }
        else if (scheme_type == "nearest-neighbour" || scheme_type == "k-nearest-neighbours") {
            timers.fspace_setup.start();
            auto fs_in = PointCloud(sgrid);
            auto fs_out = functionspace::PointCloud(tgrid, grid::MatchingPartitioner(fs_in));
            timers.fspace_setup.stop();

            timers.interpolation_setup.start();
            interpolator = Interpolation{ scheme, sgrid, tgrid };
            timers.interpolation_setup.stop();
        }
        else { 
            timers.fspace_setup.start();
            auto fs_in = StructuredColumns(sgrid, scheme);
            auto fs_out = functionspace::PointCloud(tgrid, grid::MatchingPartitioner(fs_in));
            timers.fspace_setup.stop();

            timers.interpolation_setup.start();
            interpolator = Interpolation{ scheme, fs_in, fs_out };
            timers.interpolation_setup.stop();
        }
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
    if (scheme_str == "linear") {
        scheme.set("type", "structured-linear2D");
        // The stencil does not require any halo, but we set it to 1 for pole treatment!
        scheme.set("halo", 1);
    }
    if (scheme_str == "cubic") {
        scheme.set("type", "structured-cubic2D");
        scheme.set("halo", 2);
    }
    if (scheme_str == "qcubic") {
        scheme.set("type", "structured-quasicubic2D");
        scheme.set("halo", 2);
    }
    if (scheme_str == "cons") {
        scheme.set("type", "conservative-spherical-polygon");
        scheme.set("order", 1);
        scheme.set("src_cell_data", false);
        scheme.set("tgt_cell_data", false);
    }
    if (scheme_str == "cons2") {
        scheme.set("type", "conservative-spherical-polygon");
        scheme.set("order", 2);
        scheme.set("src_cell_data", false);
        scheme.set("tgt_cell_data", false);
    }
    if (scheme_str == "finite-element") {
        scheme.set("type", "finite-element");
        scheme.set("halo", 1);
    }
    if (scheme_str == "nn") {
        scheme.set("type", "nearest-neighbour");
        scheme.set("halo", 1);
    }
    if (scheme_str == "knn") {
        scheme.set("type", "k-nearest-neighbours");
        scheme.set("k-nearest-neighbours", 4);
        scheme.set("halo", 1);
    }
    if (scheme_str == "grid-box") {
        scheme.set("type", "grid-box-average");
        scheme.set("halo", 1);
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
