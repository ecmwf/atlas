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
#include <numeric>
#include <cstring>

#include "eckit/filesystem/PathName.h"

#include "atlas/array.h"
#include "atlas/field/Field.h"
#include "atlas/field/FieldSet.h"
#include "atlas/grid.h"
#include "atlas/interpolation.h"
#include "atlas/functionspace/NodeColumns.h"


#include "atlas/mesh/Mesh.h"
#include "atlas/meshgenerator.h"
#include "atlas/output/Gmsh.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/runtime/AtlasTool.h"
#include "atlas/util/CoordinateEnums.h"
#include "atlas/util/Checksum.h"


#include "atlas/io/atlas-io.h"

#include "atlas/runtime/Exception.h"
#include "atlas/runtime/Trace.h"
#include "atlas/runtime/Log.h"

using atlas::util::Config;

namespace atlas {

std::vector<double> tgt_ref;


void read_mask(const std::string& mask_name, std::vector<int>& mask) {
    ATLAS_TRACE();
    // std::size_t size;
    atlas::io::RecordReader reader(mask_name);
    // reader.read("size",size).wait();
    // ATLAS_ASSERT(mask.size() == size);
    reader.read("mask", mask);
    reader.wait();
}

std::string read_grid_name(const std::string& mask_name) {
    ATLAS_TRACE();
    // std::size_t size;
    atlas::io::RecordReader reader(mask_name);
    std::string grid_name;
    // reader.read("size",size).wait();
    // ATLAS_ASSERT(mask.size() == size);
    reader.read("grid", grid_name);
    reader.wait();
    return grid_name;
}


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

class AtlasMaskDiff : public AtlasTool {
    int execute(const AtlasTool::Args& args) override;
    std::string briefDescription() override {
        return "briefDescription\n"; }
    std::string usage() override {
        return name() + " [OPTION] ... [--help]";
    }
    std::string longDescription() override {
        return "longDescription\n"
               "\n";
    }

    int numberOfPositionalArguments() override { return -1; }
    int minimumPositionalArguments() override { return 0; }

public:
    AtlasMaskDiff(int argc, char* argv[]): AtlasTool(argc, argv) {
        // Configuration
        add_option(new SimpleOption<std::string>("mask1", ""));
        add_option(new SimpleOption<std::string>("mask2", ""));
        add_option(new SimpleOption<bool>("show", ""));
    }
};


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

int AtlasMaskDiff::execute(const AtlasTool::Args& args) {
    ATLAS_TRACE("main");

    std::string mask_name1 = args.getString("mask1");
    std::string mask_name2 = args.getString("mask2");
    std::vector<int> mask1, mask2, mask_diff;
    std::string grid_name1 = read_grid_name(mask_name1);
    std::string grid_name2 = read_grid_name(mask_name2);
    
    Log::info() << "Reading mask1 from file '" << mask_name1 << "' containing grid " << grid_name1 << std::endl;
    read_mask(args.getString("mask1"), mask1);
    Log::info() << "Reading mask2 from file '" << mask_name2 << "' containing grid " << grid_name2 << std::endl;
    read_mask(args.getString("mask2"), mask2);

    if( grid_name1 != grid_name2 ) {
        Log::warning() << "Masks are defined on different grids: " << grid_name1 << " != " << grid_name2 << std::endl;
        // return failed();
    }
    if( mask1.size() != mask2.size() ) {
        Log::error() << "Masks have different sizes: " << mask1.size() << " != " << mask2.size() << std::endl;
        return failed();
    }

    int n_diff = 0;
    for( size_t i=0; i<mask1.size(); ++i ) {
        if( mask1[i] != mask2[i] ) {
            ++n_diff;
        }
    }
    Log::info() << "Number of different mask points: " << n_diff << " out of " << mask1.size() << std::endl;
    
    if (args.getBool("show",false)) {
        Grid grid(grid_name2);
        
        Log::info() << std::setprecision(4) << std::setw(5);
        auto it = grid.lonlat().begin();
        for( size_t i=0; i<mask1.size(); ++i, ++it ) {
            if( mask1[i] != mask2[i] ) {
                Log::info() << " index " << i << ", point " << *it << " : " << std::setprecision(4) << std::setw(5) << mask1[i] << " != " << std::setprecision(4) << std::setw(5) << mask2[i] << std::endl;
            }
        }
    }

    mask_diff.resize(mask1.size());
    for (size_t i=0; i<mask1.size(); ++i) {
        if (mask1[i] && mask2[i]) {
            mask_diff[i] = 1; // overlap
        }
        if (mask1[i] && !mask2[i]) {
            mask_diff[i] = 3; // only in mask1
        }
        if (!mask1[i] && mask2[i]) {
            mask_diff[i] = 2; // only in mask2
        }
    }

    
    if (n_diff) {
        Log::info() << "Writing mask_diff to file 'mask_diff.msh'" << std::endl;
        auto mesh = Mesh( Grid(grid_name2) );
        auto fs = functionspace::NodeColumns(mesh);

        // auto mask_diff_glb = fs.createField<debug_t>(option::name("mask_diff")|option::global());
        auto mask_diff_field = fs.createField<int>(option::name("mask_diff"));
        auto mask_diff_view = array::make_view<int,1>(mask_diff_field);
        for( idx_t i=0; i<mask_diff_view.size(); ++i ) {
            if (i < mask_diff.size()) {
                mask_diff_view(i) = mask_diff[i];
            }
            else {
                mask_diff_view(i) = -1; // unassigned!?
            }
        }
        fs.haloExchange(mask_diff_field);
        auto gmsh = output::Gmsh("mask_diff.msh", util::Config("coordinates","xyz"));
        gmsh.write(mesh);
        gmsh.write(mask_diff_field);
    }


    // auto mesh = Mesh( Grid(grid_name2) );
    // auto fs = functionspace::NodeColumns(mesh);

    // using debug_t = int;

    // // auto mask_diff_glb = fs.createField<debug_t>(option::name("mask_diff")|option::global());
    // auto mask_diff_field = fs.createField<int>(option::name("mask_diff"));
    // auto mask_diff_view = array::make_view<int,1>(mask_diff_field);
    // for( idx_t i=0; i<mask_diff.size(); ++i ) {
    //     mask_diff_view(i) = mask_diff[i];
    // }
    // fs.haloExchange(mask_diff_field);

    // auto field2 = fs.createField<int>(option::name("mask2"));
    // auto field_diff = fs.createField<int>(option::name("diff"));
    // auto field_idx = fs.createField<int>(option::name("idx"));

    // std::cout << "field1_glb : " << field1_glb << std::endl;
    // std::cout << "field1 : " << field1 << std::endl;

    // auto view1_glb = array::make_view<debug_t,1>(field1_glb);
    // auto view1 = array::make_view<int,1>(field1);
    // auto view2 = array::make_view<int,1>(field2);
    // auto diff = array::make_view<int,1>(field_diff);
    // auto gidx = array::make_view<gidx_t,1>(fs.global_index());
    // auto ridx = array::make_view<idx_t,1>(fs.remote_index());
    // auto idx = array::make_view<int,1>(field_idx);
    // auto ghost = array::make_view<int,1>(fs.ghost());
    // auto lonlat = array::make_view<double,2>(fs.lonlat());
    // size_t lim = 1733000;
    // for( size_t i=0; i<mask1.size(); ++i ) {
    //     if (i < lim) {
    //         idx(i) = 0;
    //     } else {
    //         idx(i) = ridx(i) - lim;
    //     }
    //     // view1_glb(i) = lonlat(i,1);
    //     // if (ghost(i)) {
    //     //     view1_glb(i) = -180.;
    //     // }
    //     if (ghost(i) == 0) {
    //         view1_glb(i) = std::max<debug_t>(0,i+1 - 20000);
    //         view1_glb(i) = std::max<debug_t>(0,gidx(i) - 20000);
    //     }
    //     else {
    //         view1_glb(i) = 0;
    //     }
    //     view1(i) = mask1[i];
    //     view2(i) = mask2[i];
    //     diff(i) = view2(i) - view1(i);
    // }
    // auto field1_scat = fs.createField<debug_t>(option::name("mask1_scat"));
    // array::make_view<debug_t,1>(field1_scat).assign(-1);
    // fs.scatter(field1_glb, field1_scat);
    // field1_scat.set_dirty(false);
    // auto field1_scat_exch = fs.createField<debug_t>(option::name("mask1_scat_exch"));
    // field1_scat_exch.array().copy(field1_scat.array());
    // fs.haloExchange(field1_scat_exch);
    // // field1_scat.set_dirty(true);
    // // field1_scat.haloExchange();
    
    // field1.set_dirty(false);
    // field2.set_dirty(false);
    // field_diff.set_dirty(false);
    // field_idx.set_dirty(false);
    // fs.global_index().set_dirty(false);
    // fs.remote_index().set_dirty(false);
    // // fs.haloExchange(field1);
    // // fs.haloExchange(field2);

    // #if 0
    // std::string coordinates = "xy";
    // if (mesh.nodes().has_field("ij")) {
    //     coordinates = "ij";
    // }
    // output::Gmsh gmsh("masks.msh", Config("coordinates",coordinates) | Config("ghost","true"));
    // gmsh.write(mesh);
    // // gmsh.write(field1);
    // gmsh.write(field1_glb);
    // gmsh.write(field1_scat);
    // gmsh.write(field1_scat_exch);
    // // gmsh.write(field2);
    // // gmsh.write(field_diff);
    // // gmsh.write(fs.ghost());
    // // gmsh.write(fs.remote_index());
    // // gmsh.write(fs.global_index());
    // // gmsh.write(field_idx);
    // #endif
    return success();
}

//-----------------------------------------------------------------------------

} // namespace


int main(int argc, char** argv) {
    atlas::AtlasMaskDiff tool(argc, argv);
    return tool.start();
}
