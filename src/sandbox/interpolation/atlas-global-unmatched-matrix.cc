/*
 * (C) Copyright 2025 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

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
#include "atlas/interpolation/AssembleGlobalMatrix.cc"
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

namespace atlas {


// USAGE:
//
// 1) To compute and store the global interpolation matrix to the disc:
//
// atlas-global-unmatched-matrix --sgrid <sgrid> --tgrid <tgrid> --interpolation <intrp> --format <fmt>
//
// 2) To read in a stored global interpolation matrix from the disc:
//
// atlas-global-unmatched-matrix --sgrid <sgrid> --tgrid <tgrid> --interpolation <intrp> --format <fmt> --read --matrix <matrix>
//
// For available interpolations see AtlasGlobalUnmatchedMatrix::interpolation_config.
// Formats <fmt> can be 'eckit' (binary) or 'scrip' (netcdf).
//
// Using the grid API 'create_fspaces' we can hide interpolation method specific requirements
// such as which functionspace needs to be set-up.

using namespace atlas;
using Matrix = atlas::linalg::SparseMatrixStorage;

class ParInter {
public:
    ParInter(const Matrix& gmat, const FunctionSpace src_fs, const FunctionSpace tgt_fs) : 
        gmat_(gmat), src_fs_(src_fs), tgt_fs_(tgt_fs) {
            extract(rows_, cols_, gcols_, vals_);
    }

    void print() const {
        for (int task = 0; task < mpi::comm().size(); ++task) {
            if (mpi::comm().rank() == task) {
                std::cout << "TASK " << task << std::endl;
                for (int i = 0; i < rows_.size(); ++i) {
                    std::cout << "\t" << rows_[i] << ", " << cols_[i] << ", " << vals_[i];
                    if (i < gcols_.size()) {
                        std::cout << ", " << gcols_[i];
                    }
                }
            }
            mpi::comm().barrier();
        }
    }

private:      
    template <typename Value, typename Index> 
    void extract(std::vector<Index>& local_rows,  std::vector<Index>& local_cols,  std::vector<Index>& local_gcols,
        std::vector<Value>& local_vals, int mpi_root = 0) {
        std::cout << mpi::comm().rank() << " ParInter extract" << std::endl;

        Field field_fs_part_glb = tgt_fs_.createField(tgt_fs_.partition(), option::global(mpi_root));
        ATLAS_TRACE_SCOPE("gather partition") {
            tgt_fs_.gather(tgt_fs_.partition(), field_fs_part_glb);
        }
        std::vector<Index> rows, cols;
        std::vector<Value> vals;
        interpolation::distribute_global_matrix(atlas::linalg::make_host_view<Value, Index>(gmat_), field_fs_part_glb, rows, cols, vals, mpi_root);

        // Log::info() << " rows: " << rows << std::endl;
        // Log::info() << " cols: " << cols << std::endl;
        // Log::info() << " vals: " << vals << std::endl;

        std::unordered_map<gidx_t, idx_t> to_local_rows;
        ATLAS_TRACE_SCOPE("convert to local row indexing") {
            auto fs_gidx_exchanged = tgt_fs_.createField(tgt_fs_.global_index());
            fs_gidx_exchanged.array().copy(tgt_fs_.global_index());
            tgt_fs_.haloExchange(fs_gidx_exchanged);
            const auto fs_global_index = array::make_view<gidx_t, 1>(fs_gidx_exchanged);
            const auto fs_ghost = array::make_view<int,1>(tgt_fs_.ghost());

            for (idx_t r = 0; r < fs_global_index.size(); ++r) {
                auto gr = fs_global_index(r) - 1;
                if (fs_ghost(r) && to_local_rows.find(gr) != to_local_rows.end()) {
                    continue;
                }
                // Log::info() << " to_local_rows: " << gr << " to " << r << std::endl;
                to_local_rows[gr] = r;
            }
        }

        auto find_loc_col_idx = [&local_gcols](const Index& search) {
            auto& v = local_gcols;
            for (int i = 0; i < v.size(); ++i) {
                if (search == v[i]) {
                    return i;
                }
            }
            return -1;
        };

        idx_t loc_col_idx = 0;
        for(size_t idx = 0; idx < rows.size(); ++idx) {
            auto loc_row = to_local_rows[rows[idx]];
            auto glb_col = cols[idx];
            auto val = vals[idx];
            local_rows.emplace_back(loc_row);
            local_vals.emplace_back(val);
            auto found_loc_col_idx = find_loc_col_idx(glb_col);
            if (found_loc_col_idx == -1) { // NOT found
                local_gcols.emplace_back(glb_col);
                local_cols.emplace_back(loc_col_idx++);
                // Log::info() << "put loc_row, cols, gcols, vals : " << loc_row << ", " << loc_col_idx << ", " << col << ", " << val << std::endl;
            }
            else {
                local_cols.emplace_back(found_loc_col_idx);
                // Log::info() << "put loc_row, cols, gcol, vals : " << loc_row << ", " << loc_col_idx << ", " << cols[search_idx] << val << std::endl;
            }
        }
    }

private:
    const Matrix& gmat_;
    const FunctionSpace src_fs_;
    const FunctionSpace tgt_fs_;
    std::vector<eckit::linalg::Index> rows_;
    std::vector<eckit::linalg::Index> cols_;
    std::vector<eckit::linalg::Index> gcols_;
    std::vector<eckit::linalg::Scalar> vals_;

};


class AtlasGlobalUnmatchedMatrix : public AtlasTool {
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
    std::string get_matrix_name(std::string& sgrid, std::string& tgrid, std::string& interp, std::string format = std::string());
    Matrix read_matrix(std::string matrix_name, std::string format);

public:
    AtlasGlobalUnmatchedMatrix(int argc, char* argv[]): AtlasTool(argc, argv) {
        add_option(new SimpleOption<std::string>("sgrid", "source grid"));
        add_option(new SimpleOption<std::string>("tgrid", "target grid"));
        add_option(new SimpleOption<std::string>("interpolation", "interpolation methods"));
        add_option(new SimpleOption<std::string>("matrix", "name of the interpolation matrix"));
        add_option(new SimpleOption<std::string>("format", "format of the matrix output: eckit, SCRIP"));
    }
};

int AtlasGlobalUnmatchedMatrix::execute(const AtlasTool::Args& args) {
    std::string sgrid_name = "O2";
    args.get("sgrid", sgrid_name);
    std::string tgrid_name = "O2";
    args.get("tgrid", tgrid_name);

    std::string interpolation_method = "finite-element";
    args.get("interpolation", interpolation_method);
    Log::info() << "source grid: " << sgrid_name << ", target grid: " << tgrid_name << ", interpolation: " << interpolation_method << std::endl;

    std::string matrix_format = "eckit";
    args.get("format", matrix_format);

    auto sgrid = Grid{sgrid_name};
    auto tgrid = Grid{tgrid_name};

    Matrix gmatrix;
    std::string matrix_name("");
    args.get("matrix", matrix_name);
    if (matrix_name == "") {
        matrix_name = get_matrix_name(sgrid_name, tgrid_name, interpolation_method);
    }
    if (mpi::comm().rank() == 0) {
        gmatrix = read_matrix(matrix_name, matrix_format);
        auto eckit_gmatrix = atlas::linalg::make_non_owning_eckit_sparse_matrix(gmatrix);
        eckit_gmatrix.print(Log::info());
        Log::info() << std::endl;
    }
    mpi::comm().barrier();

    FunctionSpace src_fs;
    FunctionSpace tgt_fs;
    if (interpolation_method.find("nearest-neighbour") != std::string::npos) {
        interpolation_method = "finite-element";
    }
    auto scheme = create_fspaces(interpolation_method, sgrid, tgrid, src_fs, tgt_fs);
    std::cout << mpi::comm().rank() << " src_fs, tgt_fs szie " << src_fs.size() << ", " << tgt_fs.size() << std::endl;

    // auto matrix = interpolation::distribute_global_matrix(src_fs, tgt_fs, gmatrix);
    // std::cout << mpi::comm().rank() << " matrix nnz " << matrix.nnz() << std::endl;

    std::vector<eckit::linalg::Index> rows;
    std::vector<eckit::linalg::Index> cols;
    std::vector<eckit::linalg::Index> gcols;
    std::vector<eckit::linalg::Scalar> vals;

    ParInter parOp(gmatrix, src_fs, tgt_fs);

    parOp.print();

    // FunctionSpace src_fs;
    // FunctionSpace tgt_fs;
    // if (interpolation_method.find("nearest-neighbour") != std::string::npos) {
    //     interpolation_method = "finite-element";
    // }
    // auto scheme = create_fspaces(interpolation_method, sgrid, tgrid, src_fs, tgt_fs);
    // auto matrix = interpolation::distribute_global_matrix(src_fs, tgt_fs, gmatrix);

    // interpolation::MatrixCache cache(std::move(matrix));
    // auto interpolator = Interpolation(scheme, src_fs, tgt_fs, cache);

    // // Allocate and initialise own memory here to show possibilities
    // // Note: reading a field from disc is an extra feature
    // auto src_field = interpolator.source().createField<double>();
    // auto tgt_field = interpolator.target().createField<double>();
    // auto src_lonlat = array::make_view<double, 2>(interpolator.source().lonlat());
    // ATLAS_TRACE_SCOPE("initialize source") {
    //     auto src_field_v = array::make_view<double, 1>(src_field);
    //     for (idx_t i = 0; i < src_fs.size(); ++i) {
    //         src_field_v[i] = util::function::vortex_rollup(src_lonlat(i, 0), src_lonlat(i, 1), 1.);
    //     }
    // }

    // interpolator.execute(src_field, tgt_field);

    // ATLAS_TRACE_SCOPE("output from the read-in matrix") {
    //     std::string tgt_name = "tfield_" + matrix_name;
    //     output::Gmsh gmsh(tgt_name + ".msh", Config("coordinates", "lonlat") | Config("ghost", "true"));
    //     if( functionspace::NodeColumns(tgt_field.functionspace())) {
    //         Log::info() << "storing distributed remapped field '" << tgt_name << "'." << std::endl;
    //         gmsh.write(functionspace::NodeColumns(tgt_field.functionspace()).mesh());
    //         tgt_field.haloExchange();
    //         gmsh.write(tgt_field);
    //     }
    // }

    return success();
}


//-----------------------------------------------------------------------------


Config AtlasGlobalUnmatchedMatrix::interpolation_config(std::string scheme_str) {
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

Config AtlasGlobalUnmatchedMatrix::create_fspaces(const std::string& scheme_str, const Grid& input_grid, const Grid& output_grid,
        FunctionSpace& fs_in, FunctionSpace& fs_out) {
    const Config scheme = interpolation_config(scheme_str);
    auto scheme_type = scheme.getString("type");
    if (scheme_type == "finite-element" || scheme_type == "unstructured-bilinear-lonlat") {
        auto inmesh = Mesh(input_grid);
        if (input_grid.type() == "ORCA") {
            fs_in = functionspace::NodeColumns(inmesh);
        }
        else {
            fs_in = functionspace::NodeColumns(inmesh, option::halo(1));
        }
        auto partitioner = mpi::size() == 1 ? grid::Partitioner("serial") : grid::MatchingPartitioner(inmesh);
        fs_out = functionspace::PointCloud(output_grid, partitioner);
    }
    else if (scheme_type == "conservative-spherical-polygon") {
        bool src_cell_data = scheme.getBool("src_cell_data");
        bool tgt_cell_data = scheme.getBool("tgt_cell_data");
        auto tgt_mesh_config = output_grid.meshgenerator() | option::halo(0);
        auto tgt_mesh = MeshGenerator(tgt_mesh_config).generate(output_grid);
        if (tgt_cell_data) {
            fs_out = functionspace::CellColumns(tgt_mesh, option::halo(0));
        }
        else {
            fs_out = functionspace::NodeColumns(tgt_mesh, option::halo(0));
        }
        auto src_mesh_config = input_grid.meshgenerator() | option::halo(2);
        Mesh src_mesh;
        auto partitioner = mpi::size() == 1 ? grid::Partitioner("serial") : grid::MatchingPartitioner(tgt_mesh);
        src_mesh = MeshGenerator(src_mesh_config).generate(input_grid, partitioner);
        if (src_cell_data) {
            fs_in = functionspace::CellColumns(src_mesh, option::halo(2));
        }
        else {
            fs_in = functionspace::NodeColumns(src_mesh, option::halo(2));
        }
    }
    else if (scheme_type == "nearest-neighbour" || scheme_type == "k-nearest-neighbours" || scheme_type == "grid-box-average") {
        fs_in = PointCloud(input_grid);
        auto inmesh = Mesh(input_grid);
        auto partitioner = mpi::size() == 1 ? grid::Partitioner("serial") : grid::MatchingPartitioner(inmesh);
        fs_out = functionspace::PointCloud(output_grid, partitioner);
    }
    else {
        fs_in = functionspace::StructuredColumns(input_grid, scheme);
        auto partitioner = mpi::size() == 1 ? grid::Partitioner("serial") : grid::MatchingPartitioner(fs_in);
        fs_out = functionspace::PointCloud(output_grid, partitioner, scheme);
    }
    return scheme;
}

Matrix AtlasGlobalUnmatchedMatrix::read_matrix(std::string matrix_name, std::string format) {
    if (format == "eckit") {
        eckit::linalg::SparseMatrix eckit_matrix;
        Log::info() << "reading matrix '" << matrix_name + ".eckit'" << std::endl;
        eckit_matrix.load(matrix_name + ".eckit");
        return atlas::linalg::make_sparse_matrix_storage(std::move(eckit_matrix));
    }
    else {
        ATLAS_NOTIMPLEMENTED;
    }
    return Matrix{};
}


std::string AtlasGlobalUnmatchedMatrix::get_matrix_name(std::string& sgrid, std::string& tgrid,
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

} // namespace


int main(int argc, char** argv) {
    atlas::AtlasGlobalUnmatchedMatrix tool(argc, argv);
    return tool.start();
}
