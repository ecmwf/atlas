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
#include "atlas/interpolation/NonLinear.h"
#include "atlas/linalg/sparse.h"
#include "atlas/linalg/sparse/MakeEckitSparseMatrix.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/meshgenerator.h"
#include "atlas/output/Gmsh.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/parallel/Collect.h"
#include "atlas/runtime/AtlasTool.h"
#include "atlas/util/CoordinateEnums.h"
#include "atlas/util/function/VortexRollup.h"
#include "atlas/util/function/MDPI_functions.h"
#include "atlas/util/Locate.h"
#include "atlas/field/MissingValue.h"

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
    ParInter(const Matrix& gmat, const FunctionSpace src_fs, const FunctionSpace tgt_fs, const grid::Distribution& src_distribution, const grid::Distribution& tgt_distribution) : 
        gmat_(gmat), src_fs_(src_fs), tgt_fs_(tgt_fs), src_distribution_(src_distribution), tgt_distribution_(tgt_distribution) {
            ATLAS_TRACE("Parinter::setup");
            extract(rows_, cols_, gcols_, vals_);
            setup_collect(src_fs_, gcols_);

            std::size_t nr = tgt_fs.size();
            std::size_t nc = collect_size_;
            eckit::linalg::Index index_base = 0;
            bool is_sorted = false;
            lmat_ = linalg::make_sparse_matrix_storage_from_rows_columns_values(nr, nc, rows_, cols_, vals_, index_base, is_sorted);
            find_missing_rows();

            nonlinear_ = interpolation::NonLinear("missing-if-any-missing", util::Config());

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

    void execute(const Field& src, Field& tgt) {
        ATLAS_TRACE("Parinter::execute");
        auto collect_shape = src.shape();
        collect_shape[0] = collect_size_;
        std::unique_ptr<array::Array> collect_src(array::Array::create(src.datatype(),collect_shape));

        collect_.execute<double,1>(const_cast<Field&>(src), *collect_src);

        {
            auto collect_src_view = array::make_view<double,1>(*collect_src);
            auto tgt_view = array::make_view<double,1>(tgt);
            auto matrix = linalg::make_host_view<eckit::linalg::Scalar, eckit::linalg::Index>(lmat_);

            ATLAS_ASSERT(tgt_view.shape(0) == matrix.rows());

            ATLAS_TRACE_SCOPE("sparse_matrix_multiply") {
                // linalg::sparse::current_backend("eckit_linalg");
                if (nonlinear_(src)) {
                    eckit::linalg::SparseMatrix matrix_nl = atlas::linalg::make_eckit_sparse_matrix(matrix);
                    ATLAS_ASSERT(matrix_nl.cols() == collect_src->size());
                    ATLAS_ASSERT(matrix_nl.cols() == collect_src->size());
                    nonlinear_->execute(matrix_nl, src, *collect_src);
                    linalg::sparse_matrix_multiply(matrix_nl, collect_src_view, tgt_view);
                }
                else {
                    linalg::sparse_matrix_multiply(matrix, collect_src_view, tgt_view);
                }
            }

            if (not tgt.metadata().has("missing_value")) {
                field::MissingValue mv_src(src);
                if (mv_src) {
                    ATLAS_DEBUG();
                    mv_src.metadata(tgt);
                    ATLAS_ASSERT(field::MissingValue(tgt));
                }
                else if (not missing_rows_.empty()) {
                    ATLAS_DEBUG();
                    if (not tgt.metadata().has("missing_value")) {
                        tgt.metadata().set("missing_value", 9999.);
                    }
                    tgt.metadata().set("missing_value_type", "equals");
                }
            }
            
            ATLAS_TRACE_SCOPE("mask") {
                ATLAS_DEBUG();
                if (tgt.metadata().has("missing_value")) {
                    double missing_value = tgt.metadata().get<double>("missing_value");
                    for (idx_t r : missing_rows_) {
                        tgt_view(r) = missing_value; 
                    }
                }
            }
        }
    }


private:

    template <typename Vector>
    void setup_collect(FunctionSpace fs, const Vector& global_index) {
            ATLAS_TRACE();
            collect_size_ = global_index.size();

            std::vector<gidx_t> collect_gidx(collect_size_);
            std::vector<int>    collect_partition(collect_size_);
            std::vector<idx_t>  collect_ridx(collect_size_);

            for (size_t i=0; i<global_index.size(); ++i) {
                collect_gidx[i] = global_index[i] + 1;
            }
            idx_t ridx_base = 0;

            util::locate(fs, src_distribution_, collect_gidx, collect_partition, collect_ridx, ridx_base);
            collect_.setup(collect_size_, collect_partition.data(), collect_ridx.data(), ridx_base);
    }

    template <typename Value, typename Index> 
    void extract(std::vector<Index>& local_rows,  std::vector<Index>& local_cols,  std::vector<Index>& local_gcols,
        std::vector<Value>& local_vals, int mpi_root = 0) {
        ATLAS_TRACE();
        // std::cout << mpi::comm().rank() << " ParInter extract" << std::endl;

        // Field field_fs_part_glb = tgt_fs_.createField(tgt_fs_.partition(), option::global(mpi_root));
        // ATLAS_TRACE_SCOPE("gather partition") {
        //     tgt_fs_.gather(tgt_fs_.partition(), field_fs_part_glb);
        // }
        std::vector<Index> rows, cols;
        std::vector<Value> vals;
        interpolation::distribute_global_matrix_as_triplets(
            // array::make_view<int,1>(field_fs_part_glb).data(), 
            tgt_distribution_,
            atlas::linalg::make_host_view<Value, Index>(gmat_), rows, cols, vals, mpi_root);

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

            std::unordered_map<gidx_t,idx_t> map_gcol_to_loc_col_idx;
            auto find_loc_col_idx = [&map_gcol_to_loc_col_idx](const gidx_t gcol) {
                auto it = map_gcol_to_loc_col_idx.find(gcol);
                if (it == map_gcol_to_loc_col_idx.end()) {
                    return -1;
                }
                return it->second;
            };

            idx_t loc_col_idx = 0;
            auto new_loc_col_idx = [&](const gidx_t gcol) {
                idx_t lcol = loc_col_idx;
                map_gcol_to_loc_col_idx[gcol] = lcol;
                ++loc_col_idx;
                return lcol;
            };

            for(size_t idx = 0; idx < rows.size(); ++idx) {
                auto loc_row = to_local_rows[rows[idx]];
                auto glb_col = cols[idx];
                auto val = vals[idx];
                local_rows.emplace_back(loc_row);
                local_vals.emplace_back(val);
                auto found_loc_col_idx = find_loc_col_idx(glb_col);
                if (found_loc_col_idx >= 0) {
                    local_cols.emplace_back(found_loc_col_idx);
                }
                else {
                    local_cols.emplace_back(new_loc_col_idx(glb_col));
                    local_gcols.emplace_back(glb_col);
                }
            }
        }
    }

    void find_missing_rows() {
        auto matrix = linalg::make_host_view<eckit::linalg::Scalar, eckit::linalg::Index>(lmat_);
        for(std::size_t r = 0; r < matrix.rows(); ++r) {
            int cols = matrix.outer()[r + 1] - matrix.outer()[r];
            if (cols == 0) {
                missing_rows_.emplace_back(r);
            }
        }
    }

private:
    const Matrix& gmat_;
    const FunctionSpace src_fs_;
    const FunctionSpace tgt_fs_;
    const grid::Distribution src_distribution_;
    const grid::Distribution tgt_distribution_;
    parallel::Collect collect_;
    std::vector<eckit::linalg::Index> rows_;
    std::vector<eckit::linalg::Index> cols_;
    std::vector<eckit::linalg::Index> gcols_;
    std::vector<eckit::linalg::Scalar> vals_;
    Matrix lmat_;
    std::size_t collect_size_;

    interpolation::NonLinear nonlinear_;
    
    std::vector<idx_t> missing_rows_;

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
        FunctionSpace& fs_in, FunctionSpace& fs_out, grid::Distribution& dist_in, grid::Distribution& dist_out);
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
    ATLAS_TRACE("main");
    std::string sgrid_name = "O2";
    args.get("sgrid", sgrid_name);
    std::string tgrid_name = "O2";
    args.get("tgrid", tgrid_name);

    std::string interpolation_method = "finite-element";
    args.get("interpolation", interpolation_method);

    Log::info() << "Setup:" << std::endl;
    Log::info() << "    source grid:   " << sgrid_name << std::endl;
    Log::info() << "    target grid:   " << tgrid_name << std::endl;
    Log::info() << "    interpolation: " << interpolation_method << std::endl;

    std::string matrix_format = "eckit";
    args.get("format", matrix_format);


    Grid src_grid, tgt_grid;    
    FunctionSpace src_fs, tgt_fs;
    grid::Distribution src_dist, tgt_dist;

    std::unique_ptr<ParInter> interpolator;

    std::string matrix_name("");
    args.get("matrix", matrix_name);
    if (matrix_name == "") {
        matrix_name = get_matrix_name(sgrid_name, tgrid_name, interpolation_method);
    }
    if (not eckit::PathName(matrix_name+".eckit").exists()) {
        Log::error() << "\n ERROR: Matrix with file path " << matrix_name << ".eckit does not exist\n" << std::endl;
        return failed();
    }
    ATLAS_TRACE_SCOPE("Setup") {

        Matrix gmatrix;
        ATLAS_TRACE_SCOPE("Read matrix on rank 0") {
            if (mpi::comm().rank() == 0) {
                gmatrix = read_matrix(matrix_name, matrix_format);
            }
        }
        mpi::comm().barrier();

        ATLAS_TRACE_SCOPE("Setup source and target") {
            if (interpolation_method.find("nearest-neighbour") != std::string::npos) {
                interpolation_method = "finite-element";
            }
            src_grid = Grid{sgrid_name};
            tgt_grid = Grid{tgrid_name};
            create_fspaces(interpolation_method, src_grid, tgt_grid, src_fs, tgt_fs, src_dist, tgt_dist);
            std::cout << "[" << mpi::comm().rank() << "] src_fs.size() = " << src_fs.size() << "\t tgt_fs.size() = " << tgt_fs.size() << std::endl;
        }

        mpi::comm().barrier();

        interpolator = std::make_unique<ParInter>(gmatrix, src_fs, tgt_fs, src_dist, tgt_dist);
    }

    auto src_field = src_fs.createField<double>();
    auto tgt_field = tgt_fs.createField<double>();

    auto src_lonlat = array::make_view<double, 2>(src_fs.lonlat());
    ATLAS_TRACE_SCOPE("initialize source") {
        auto src_field_v = array::make_view<double, 1>(src_field);
        ATLAS_DEBUG_VAR(src_grid.type());
        if (src_grid.type() == "ORCA") {
            std::vector<int> is_water(src_grid.size(), 1);
            if (eckit::PathName("orca_mask").exists()) {
                std::ifstream orca_mask("orca_mask");
                for (int i=0; i<is_water.size(); ++i) {
                    orca_mask >> is_water[i];
                }
            }
            auto glb_idx = array::make_view<gidx_t,1>(src_fs.global_index());
            auto mask = [&](idx_t j) -> bool {
                return not is_water[glb_idx(j)-1];
            };
            double missing_value = 9999.;
            src_field.metadata().set("missing_value", missing_value);
            src_field.metadata().set("missing_value_type", "equals");
            for (idx_t i = 0; i < src_fs.size(); ++i) {
                src_field_v[i] = util::function::MDPI_gulfstream(src_lonlat(i, 0), src_lonlat(i, 1))
                               + util::function::MDPI_vortex(src_lonlat(i, 0), src_lonlat(i, 1));
                if (mask(i)) {
                   src_field_v[i] = missing_value;
                }
            }
        }
        else {
            for (idx_t i = 0; i < src_fs.size(); ++i) {
                src_field_v[i] = util::function::vortex_rollup(src_lonlat(i, 0), src_lonlat(i, 1), 1.);
                // src_field_v[i] = util::function::MDPI_gulfstream(src_lonlat(i, 0), src_lonlat(i, 1));
                // src_field_v[i] = util::function::MDPI_vortex(src_lonlat(i, 0), src_lonlat(i, 1));
            }
        }
    }

    interpolator->execute(src_field, tgt_field);

    ATLAS_TRACE_SCOPE("Output interpolated field") {
        tgt_field.set_dirty(true);
        tgt_field.haloExchange();
        std::string tgt_name = "tfield_" + matrix_name;
        output::Gmsh gmsh(tgt_name + ".msh", Config("coordinates", "lonlat") | Config("ghost", "true"));
        if( functionspace::NodeColumns(tgt_field.functionspace())) {
            Log::info() << "storing distributed remapped field '" << tgt_name << "'." << std::endl;
            gmsh.write(functionspace::NodeColumns(tgt_field.functionspace()).mesh());
            gmsh.write(tgt_field);
        }
        else {
            auto tgt_mesh = Mesh(tgt_grid,tgt_dist);
            gmsh.write(tgt_mesh);
            gmsh.write(tgt_field, functionspace::NodeColumns(tgt_mesh));
        }
    }

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
        FunctionSpace& fs_in, FunctionSpace& fs_out, grid::Distribution& dist_in, grid::Distribution& dist_out) {
    const Config scheme = interpolation_config(scheme_str);
    auto scheme_type = scheme.getString("type");
    //if (scheme_type == "finite-element" || scheme_type == "unstructured-bilinear-lonlat") {
    ATLAS_TRACE_SCOPE("source functionspace") {
        // dist_in = grid::Distribution(input_grid, grid::Partitioner{input_grid.partitioner()});
        dist_in = grid::Distribution(input_grid, grid::Partitioner{"regular_bands"});
        if( input_grid.type() == "ORCA") {
            // fs_in = functionspace::NodeColumns(Mesh(input_grid, dist_in));
            fs_in = functionspace::PointCloud(input_grid, dist_in);
        }
        else {
            fs_in = functionspace::PointCloud(input_grid, dist_in);
        }
    }
    ATLAS_TRACE_SCOPE("target functionspace") {
        //auto partitioner = mpi::size() == 1 ? grid::Partitioner("serial") : grid::MatchingPartitioner(inmesh);
        // auto partitioner = mpi::size() == 1 ? grid::Partitioner("serial") : grid::Partitioner("regular_bands");
        auto partitioner = mpi::size() == 1 ? grid::Partitioner("serial") : grid::Partitioner("equal_regions");
        dist_out = grid::Distribution(output_grid, partitioner);
        fs_out = functionspace::PointCloud(output_grid, dist_out);
        // fs_out = functionspace::NodeColumns(Mesh(output_grid, dist_out));
    }
    //}
#if 0
    else if (scheme_type == "conservative-spherical-polygon") {
        bool src_cell_data = scheme.getBool("src_cell_data");
        bool tgt_cell_data = scheme.getBool("tgt_cell_data");
        auto tgt_mesh_config = output_grid.meshgenerator() | option::halo(0);
        auto dist_out = 
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
        dist_in = grid::Distribution(input_grid, partitioner);
        src_mesh = MeshGenerator(src_mesh_config).generate(input_grid, dist_in);
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
#endif
    return scheme;
}

Matrix AtlasGlobalUnmatchedMatrix::read_matrix(std::string matrix_name, std::string format) {
    ATLAS_TRACE();
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
