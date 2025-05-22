/*
 * (C) Copyright 2025 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */


#include "atlas/coupler/ParInter.h"
#include "atlas/io/atlas-io.h"

#include "tests/AtlasTestEnvironment.h"

using atlas::functionspace::PointCloud;
using atlas::functionspace::NodeColumns;
using atlas::functionspace::StructuredColumns;
using atlas::util::Config;
using Matrix = atlas::linalg::SparseMatrixStorage;

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
        add_option(new SimpleOption<std::string>("source", "atlas-io source file"));
    }
};

int AtlasGlobalUnmatchedMatrix::execute(const AtlasTool::Args& args) {
    ATLAS_TRACE("main");
    std::string sgrid_name = "O2";
    args.get("sgrid", sgrid_name);

    std::vector<double> src_values;
    bool read_source = false;
    if (args.has("source")) {
        read_source = true;
        std::string atlas_io_file;
        args.get("source", atlas_io_file);
        atlas::io::RecordReader atlas_io_reader(atlas_io_file);
        std::string grid;
        atlas_io_reader.read("grid.name",sgrid_name);
        if( atlas::mpi::rank() == 0 ) {
            atlas_io_reader.read("fields[0].array",src_values);
        }
        atlas_io_reader.wait();
    }


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

        double missing_value = 9999.;

        if (read_source) {
            auto field_glb = src_fs.createField(src_field, option::global(0));
            if (mpi::rank() == 0) {
                auto field_glb_view = array::make_view<double,1>(field_glb);
                for (idx_t i=0; i<field_glb_view.size(); ++i) {
                    if( std::abs(src_values[i]-missing_value) < 1.e-4 ) {
                        field_glb_view[i] = missing_value;
                    }
                    else {
                        field_glb_view[i] = std::abs(src_values[i]);
                    }
                }
            }
            src_fs.scatter(field_glb,src_field);
        }
        else {
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
        src_field.metadata().set("missing_value", missing_value);
        src_field.metadata().set("missing_value_type", "equals");
    }

    interpolator->execute(src_field, tgt_field);

    ATLAS_TRACE_SCOPE("Output interpolated field") {
        tgt_field.set_dirty(true);
        tgt_field.haloExchange();
        std::string tgt_name = "tfield_" + matrix_name;
        output::Gmsh gmsh(tgt_name + ".msh", Config("coordinates", "lonlat") | Config("ghost", "true"));
        Log::info() << "storing distributed remapped field '" << tgt_name << ".msh'." << std::endl;
        if( functionspace::NodeColumns(tgt_field.functionspace())) {
            gmsh.write(functionspace::NodeColumns(tgt_field.functionspace()).mesh());
            gmsh.write(tgt_field);
        }
        else if( functionspace::CellColumns(tgt_field.functionspace())) {
            gmsh.write(functionspace::CellColumns(tgt_field.functionspace()).mesh());
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
        if( output_grid.type() == "healpix") {
            auto partitioner = grid::Partitioner{output_grid.partitioner()};
            dist_out = grid::Distribution(output_grid, partitioner);
            fs_out = functionspace::CellColumns(Mesh(output_grid, partitioner));
        }
        else if( output_grid.type() == "ORCA") {
            auto partitioner = grid::Partitioner{output_grid.partitioner()};
            // auto partitioner = mpi::size() == 1 ? grid::Partitioner("serial") : grid::Partitioner("regular_bands");
            dist_out = grid::Distribution(output_grid, partitioner);
        // fs_out = functionspace::NodeColumns(Mesh(output_grid, dist_out));
            fs_out = functionspace::PointCloud(output_grid, dist_out);
        }
        else {
            auto partitioner = grid::Partitioner{output_grid.partitioner()};
            dist_out = grid::Distribution(output_grid, partitioner);
            fs_out = functionspace::PointCloud(output_grid, dist_out);
        }
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
