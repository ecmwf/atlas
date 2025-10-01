/*
 * (C) Copyright 2025 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "HDF5Reader.h"

#include "atlas/array.h"
#include "atlas/field/Field.h"
#include "atlas/runtime/AtlasTool.h"
#include "atlas/runtime/Exception.h"

#include "atlas/grid.h"
#include "atlas/mesh.h"
#include "atlas/output/Gmsh.h"
#include "atlas/functionspace/BlockStructuredColumns.h"
#include "atlas/functionspace/NodeColumns.h"
#include "atlas/interpolation.h"

namespace atlas {


namespace{

template <class ValueType>
void block_copy(const Field sloc, Field loc, const functionspace::BlockStructuredColumns& fs) {
    if (sloc.variables() and sloc.levels()) {
        auto loc_v  = array::make_view<ValueType, 4>(loc);
        auto sloc_v = array::make_view<ValueType, 3>(sloc);
        for (idx_t jblk = 0; jblk < fs.nblks(); ++jblk) {
            auto b = fs.block(jblk);
            for (idx_t jvar = 0; jvar < sloc.shape(2); ++jvar) {
                for (idx_t jlev = 0; jlev < sloc.shape(1); ++jlev) {
                    for (idx_t jrof = 0; jrof < b.size(); ++jrof) {
                        loc_v(jblk, jvar, jlev, jrof) = sloc_v(b.index(jrof), jlev, jvar);
                    }
                }
            }
        }
    }
    else if (not sloc.variables() and sloc.levels()) {
        auto loc_v  = array::make_view<ValueType, 3>(loc);
        auto sloc_v = array::make_view<ValueType, 2>(sloc);
        for (idx_t jblk = 0; jblk < fs.nblks(); ++jblk) {
            auto b = fs.block(jblk);
            for (idx_t jlev = 0; jlev < sloc.shape(1); ++jlev) {
                for (idx_t jrof = 0; jrof < b.size(); ++jrof) {
                    loc_v(jblk, jlev, jrof) = sloc_v(b.index(jrof), jlev);
                }
            }
        }
    }
    else if (sloc.variables() and not sloc.levels()) {
        auto loc_v  = array::make_view<ValueType, 3>(loc);
        auto sloc_v = array::make_view<ValueType, 2>(sloc);
        for (idx_t jblk = 0; jblk < fs.nblks(); ++jblk) {
            auto b = fs.block(jblk);
            for (idx_t jvar = 0; jvar < sloc.shape(1); ++jvar) {
                for (idx_t jrof = 0; jrof < b.size(); ++jrof) {
                    loc_v(jblk, jvar, jrof) = sloc_v(b.index(jrof), jvar);
                }
            }
        }
    }
    else {
        auto loc_v  = array::make_view<ValueType, 2>(loc);
        auto sloc_v = array::make_view<ValueType, 1>(sloc);
        for (idx_t jblk = 0; jblk < fs.nblks(); ++jblk) {
            auto b = fs.block(jblk);
            for (idx_t jrof = 0; jrof < b.size(); ++jrof) {
                loc_v(jblk, jrof) = sloc_v(b.index(jrof));
            }
        }
    }
}

template <class ValueType>
void rev_block_copy(const Field loc, Field sloc, const functionspace::BlockStructuredColumns& fs) {
    if (loc.variables() and loc.levels()) {
        auto loc_v  = array::make_view<ValueType, 4>(loc);
        auto sloc_v = array::make_view<ValueType, 3>(sloc);
        for (idx_t jblk = 0; jblk < fs.nblks(); ++jblk) {
            auto b = fs.block(jblk);
            for (idx_t jvar = 0; jvar < sloc.shape(2); ++jvar) {
                for (idx_t jlev = 0; jlev < sloc.shape(1); ++jlev) {
                    for (idx_t jrof = 0; jrof < b.size(); ++jrof) {
                        sloc_v(b.index(jrof), jlev, jvar) = loc_v(jblk, jvar, jlev, jrof);
                    }
                }
            }
        }
    }
    else if (not loc.variables() and loc.levels()) {
        auto loc_v  = array::make_view<ValueType, 3>(loc);
        auto sloc_v = array::make_view<ValueType, 2>(sloc);
        for (idx_t jblk = 0; jblk < fs.nblks(); ++jblk) {
            auto b = fs.block(jblk);
            for (idx_t jlev = 0; jlev < sloc.shape(1); ++jlev) {
                for (idx_t jrof = 0; jrof < b.size(); ++jrof) {
                    sloc_v(b.index(jrof), jlev) = loc_v(jblk, jlev, jrof);
                }
            }
        }
    }
    else if (loc.variables() and not loc.levels()) {
        auto loc_v  = array::make_view<ValueType, 3>(loc);
        auto sloc_v = array::make_view<ValueType, 2>(sloc);
        for (idx_t jblk = 0; jblk < fs.nblks(); ++jblk) {
            auto b = fs.block(jblk);
            for (idx_t jvar = 0; jvar < sloc.shape(1); ++jvar) {
                for (idx_t jrof = 0; jrof < b.size(); ++jrof) {
                    sloc_v(b.index(jrof), jvar) = loc_v(jblk, jvar, jrof);
                }
            }
        }
    }
    else {
        auto loc_v  = array::make_view<ValueType, 2>(loc);
        auto sloc_v = array::make_view<ValueType, 1>(sloc);
        for (idx_t jblk = 0; jblk < fs.nblks(); ++jblk) {
            auto b = fs.block(jblk);
            for (idx_t jrof = 0; jrof < b.size(); ++jrof) {
                sloc_v(b.index(jrof)) = loc_v(jblk, jrof);
            }
        }
    }
}

void copy_nonblocked_to_blocked(const Field& nonblocked, Field& blocked) {
    auto fs = blocked.functionspace();
    auto kind = nonblocked.datatype().kind();
    if (kind == array::DataType::kind<int>()) {
        block_copy<int>(nonblocked, blocked, fs);
    }
    else if (kind == array::DataType::kind<long>()) {
        block_copy<long>(nonblocked, blocked, fs);
    }
    else if (kind == array::DataType::kind<float>()) {
        block_copy<float>(nonblocked, blocked, fs);
    }
    else if (kind == array::DataType::kind<double>()) {
        block_copy<double>(nonblocked, blocked, fs);
    }
    else {
        throw_Exception("datatype not supported", Here());
    }
}

void copy_blocked_to_nonblocked(const Field& blocked, Field& nonblocked) {
    auto fs = blocked.functionspace();
    auto kind = blocked.datatype().kind();
    if (kind == array::DataType::kind<int>()) {
        rev_block_copy<int>(blocked, nonblocked, fs);
    }
    else if (kind == array::DataType::kind<long>()) {
        rev_block_copy<long>(blocked, nonblocked, fs);
    }
    else if (kind == array::DataType::kind<float>()) {
        rev_block_copy<float>(blocked, nonblocked, fs);
    }
    else if (kind == array::DataType::kind<double>()) {
        rev_block_copy<double>(blocked, nonblocked, fs);
    }
    else {
        throw_Exception("datatype not supported", Here());
    }
}
}



class Program : public AtlasTool {
    int execute(const AtlasTool::Args& args) override;
    std::string briefDescription() override { return "brief description"; }
    std::string usage() override {
        return name() + " [OPTION] ... [--help]";
    }

    int numberOfPositionalArguments() override { return -1; }
    int minimumPositionalArguments() override { return 0; }

public:
    Program(int argc, char* argv[]): AtlasTool(argc, argv) {
        add_option(new SimpleOption<std::string>("file", "input file"));
        add_option(new SimpleOption<std::string>("dataset", "dataset"));
        add_option(new SimpleOption<long>("index", "dataset"));
        add_option(new SimpleOption<std::string>("rad.grid", "target grid to interpolate to"));
        add_option(new SimpleOption<std::string>("rad.nproma", "target grid to nproma interpolate to"));
        add_option(new SimpleOption<std::string>("gmsh.coordinates", "coordinates for gmsh output [lonlat, xyz]"));
    }
};

util::Metadata hdf5_to_metadata(const std::string& file, const std::string& dataset) {
    util::Metadata metadata;
    int error = 0;
    int mpi_root = 0;
    if (mpi::rank() == mpi_root) {
        // Open the HDF5 file
        HDF5_File h5_file(file);
        if (!h5_file.is_open()) {
            error = 1;
            goto endif;
        }

        // Get the dataset from the file
        HDF5_Dataset h5_dataset = h5_file.dataset(dataset);
        if (!h5_dataset.is_open()) {
            error = 1;
            goto endif;
        }
        Log::info() << "Dataset name     : " << h5_dataset.name() << '\n';
        Log::info() << "Dataset rank     : " << h5_dataset.rank() << '\n';
        Log::info() << "Dataset extents  : " << h5_dataset.extents() << '\n';
        Log::info() << "Dataset datatype : " << h5_dataset.datatype().str() << '\n';

        int nproma = h5_dataset.extents()[h5_dataset.rank()-1];
        int nblks  = h5_dataset.extents()[0];
        int nflds  = h5_dataset.rank() == 3 ? h5_dataset.extents()[1] : 1;

        Log::info() << "nproma : " << nproma << '\n';
        Log::info() << "nblks  : " << nblks << '\n';
        Log::info() << "nflds   : " << nflds << '\n';
        Log::info() << "nblk * nproma : " << nblks * nproma << '\n';

        std::map<size_t,std::string> ngptot_to_gridname;
        ngptot_to_gridname.emplace(2048,  "F16");
        ngptot_to_gridname.emplace(6114,  "N32");
        ngptot_to_gridname.emplace(35718, "N80");
        std::string grid;
        for (const auto & [k,v] : ngptot_to_gridname) {
            if (nproma * nblks >= k && nproma * (nblks-1) < k) {
                grid = v;
                break;
            }
        }
        if (grid.empty()) {
            Log::error() << "ERROR: Could not detect grid" << std::endl;
            error = 1;
            goto endif;
        }
        metadata.set("nblks",  nblks);
        metadata.set("nproma", nproma);
        metadata.set("nflds",  nflds);
        metadata.set("grid",   grid);
        metadata.set("kind",   h5_dataset.datatype().kind());
    } endif:
    mpi::comm().broadcast(error, mpi_root);
    if (error) {
        ATLAS_THROW_EXCEPTION("Could not extract metadata from HDF5(file="<<file<<", dataset="<<dataset<<")");
    }
    metadata.broadcast(mpi_root);
    return metadata;
}

void hdf5_read(const std::string& file, const std::string& dataset, Field& field) {
    ATLAS_TRACE();
    int mpi_root = 0;
    functionspace::BlockStructuredColumns fs{field.functionspace()};
    Field field_glb = fs.createField(option::name(field.name()) | option::levels(field.levels()) | option::datatype(field.datatype()) | option::global(mpi_root));

    int error = 0;
    if (mpi::rank() == mpi_root) {
        // Open the HDF5 file
        HDF5_File h5_file(file);
        if (!h5_file.is_open()) {
            error = 1;
            goto endif;
        }

        // Get the dataset from the file
        HDF5_Dataset h5_dataset = h5_file.dataset(dataset);
        if (!h5_dataset.is_open()) {
            error = 1;
            goto endif;
        }

        Field dataset_field(h5_dataset.name(), h5_dataset.datatype(), h5_dataset.extents());
        h5_dataset.read(dataset_field);

        switch( dataset_field.datatype().kind() ) {
            case DataType::KIND_REAL64 : {
                auto glb_view = array::make_view<double,2>(field_glb);
                if (dataset_field.rank() == 2) {
                    auto dataset_view =array::make_view<double,2>(dataset_field);
                    idx_t n = 0;
                    for (idx_t jblk=0; jblk < dataset_view.shape(0); ++jblk) {
                        for (idx_t jrof=0; jrof < dataset_view.shape(1); ++jrof, ++n) {
                            if (n < glb_view.shape(0)) {
                                glb_view(n, 0) = dataset_view(jblk, jrof);
                            }
                        }
                    }
                }
                else if (dataset_field.rank() == 3) {
                    auto dataset_view =array::make_view<double,3>(dataset_field);
                    idx_t n = 0;
                    for (idx_t jblk=0; jblk < dataset_view.shape(0); ++jblk) {
                        for (idx_t jrof=0; jrof < dataset_view.shape(2); ++jrof, ++n) {
                            if (n < glb_view.shape(0)) {
                                for (idx_t jfld=0; jfld < glb_view.shape(1); ++jfld) {
                                    glb_view(n, jfld) = dataset_view(jblk, jfld, jrof);
                                }
                            }
                        }
                    }
                }
                else {
                    ATLAS_NOTIMPLEMENTED;
                }
                break;
            }
            default:
                ATLAS_THROW_EXCEPTION("Datatype " << dataset_field.datatype().str() << " not implemented");
        }
    } endif:
    mpi::comm().broadcast(error, mpi_root);
    if (error) {
        ATLAS_THROW_EXCEPTION("Could not extract field data from HDF5(file="<<file<<", dataset="<<dataset<<")");
    }
    fs.scatter(field_glb, field);
}


int Program::execute(const AtlasTool::Args& args) {
    ATLAS_TRACE("main");
    std::string input_file   = args.getString("file");
    std::string dataset = args.getString("dataset");

    util::Metadata metadata = hdf5_to_metadata(input_file, dataset);
    Log::info() << metadata << std::endl;

    int nproma = metadata.getLong("nproma");
    int nflds  = metadata.getLong("nflds");
    DataType datatype(metadata.getLong("kind"));
    Grid IFS_grid(metadata.getString("grid"));

    functionspace::BlockStructuredColumns IFS_fs(IFS_grid, util::Config("nproma",nproma));
    Field IFS_f = IFS_fs.createField(option::name(dataset)|option::datatype(datatype)|option::levels(nflds));
    hdf5_read(input_file, dataset, IFS_f);

    ATLAS_DEBUG_VAR(IFS_f);

    // Now this needs to go to different grid with different nproma
    auto rad_grid   = Grid(args.getString("rad.grid",IFS_grid.name()));
    auto rad_nproma = args.getLong("rad.nproma",nproma);
    functionspace::BlockStructuredColumns rad_fs(rad_grid, grid::MatchingPartitioner(IFS_fs), util::Config("nproma",rad_nproma));
    Field rad_f = rad_fs.createField(IFS_f);

    ATLAS_TRACE_SCOPE("Remap IFS with nproma to rad with nproma")
    {
#if 0
        auto src_mesh = Mesh(IFS_grid);
        auto int_src_fs = functionspace::NodeColumns{src_mesh, option::halo(1)};
        auto int_src_fs = functionspace::StructuredColumns{IFS_grid, option::halo(2)};
        auto tgt_mesh = Mesh(rad_grid, grid::MatchingPartitioner(IFS_fs));
        auto int_tgt_fs = functionspace::NodeColumns{tgt_mesh}; // This could be another fs
        auto interpolation = Interpolation(option::type("finite-element"),  int_src_fs, int_tgt_fs);
#endif

#if 1
        auto rad_mesh = Mesh(rad_grid, grid::MatchingPartitioner(IFS_fs));
        auto int_IFS_fs = functionspace::StructuredColumns{IFS_grid, option::halo(2)};
        auto int_rad_fs = functionspace::NodeColumns{rad_mesh}; // This could be another fs
        auto interpolation = Interpolation(option::type("structured-bicubic"),  int_IFS_fs, int_rad_fs);
#endif

#if 0
        auto IFS_mesh = Mesh(IFS_grid);
        auto rad_mesh = Mesh(rad_grid, grid::MatchingPartitioner(IFS_fs));
        auto int_IFS_fs = functionspace::NodeColumns{IFS_mesh, option::halo(1)};
        auto int_rad_fs = functionspace::NodeColumns{rad_mesh};
        auto interpolation = Interpolation(option::type("conservative-spherical-polygon"),  int_IFS_fs, int_rad_fs);
#endif

        auto int_IFS_f = int_IFS_fs.createField(IFS_f);
        copy_blocked_to_nonblocked(IFS_f, int_IFS_f);
        int_IFS_f.haloExchange();
        auto int_rad_f = int_rad_fs.createField(int_IFS_f);

        ATLAS_TRACE_SCOPE("interpolation.execute")
        interpolation.execute(int_IFS_f, int_rad_f);

        // Output Gmsh, only works if we have a mesh for radiation
        if (true) {
            int_rad_f.haloExchange();
            auto rad_gmsh = output::Gmsh("rad.msh", util::Config("coordinates",args.getString("gmsh.coordinates","lonlat")));
            rad_gmsh.write(rad_mesh);
            rad_gmsh.write(int_rad_f);
        }

        copy_nonblocked_to_blocked(int_rad_f, rad_f);
    }
    ATLAS_DEBUG_VAR(rad_f);
    return success();
}

//-----------------------------------------------------------------------------

} // namespace


int main(int argc, char** argv) {
    atlas::Program tool(argc, argv);
    return tool.start();
}
