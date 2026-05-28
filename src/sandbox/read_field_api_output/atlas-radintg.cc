/*
 * (C) Copyright 2025- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */
#define NDEBUG

#include <chrono>
#include <iomanip>

#include "HDF5Reader.h"
#include "transpositions.h"
#include "interpolate.h"
#include "output_gmsh.h"

#include "atlas/array.h"
#include "atlas/field/Field.h"
#include "atlas/runtime/AtlasTool.h"
#include "atlas/runtime/Exception.h"

#include "atlas/grid.h"
#include "atlas/mesh.h"
#include "atlas/output/Gmsh.h"
#include "atlas/functionspace/BlockStructuredColumns.h"
#include "atlas/functionspace/NodeColumns.h"
#include "atlas/functionspace/PointCloud.h"

#include "atlas/field/MultiField.h"
namespace atlas {

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
        add_option(new SimpleOption<std::string>("config-file", "configuration file"));
        add_option(new SimpleOption<std::string>("data-file", "input file (hdf5)"));
        add_option(new SimpleOption<std::string>("rad.grid", "target grid to interpolate to"));
        add_option(new SimpleOption<std::string>("rad.nproma", "target grid to nproma interpolate to"));
        add_option(new SimpleOption<std::string>("gmsh.coordinates", "coordinates for gmsh output [lonlat, xyz]"));
        add_option(new SimpleOption<std::string>("interpolation", "coordinates for gmsh output [lonlat, xyz]"));
        add_option(new SimpleOption<long>("niter", "Number of iterations"));
        add_option(new SimpleOption<bool>("output.gmsh", "coordinates for gmsh output [lonlat, xyz]"));
        add_option(new SimpleOption<bool>("multifield", "Use a multifield for the input fieldset"));
        add_option(new SimpleOption<bool>("on-device", "Do the interpolation on device"));
    }
};

auto to_upper = [](const std::string& str) {
    std::string upper{str};
    for (auto & c: upper) c = std::toupper(c);
    return upper;
};


util::Config get_subconfiguration(const eckit::Configuration& config, const std::string& key) {
    util::Config c;
    config.get(key, c);
    return c;
}

void read_field_with_timing(const std::string& data_file, Field& field) {
    auto start = std::chrono::steady_clock::now();
    HDF5Reader{data_file, field.name()}.read_field(field);
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start);
    auto& out = Log::info();
    out << "    Field " << std::left << std::setw(10) << field.name() << std::right << " took " << elapsed.count() << " ms" << std::endl;
}

std::tuple<util::Config, std::vector<std::string>, std::map<std::string,util::Metadata>>
extract_variables(const util::Config& zrgp_in, const std::string& data_file) {
    ATLAS_TRACE("Reading metadata");
    Log::info() << "Extracting variable information from data file " << data_file << " :" << std::endl;
    auto zrgp_in_vars = zrgp_in.getSubConfiguration("variables").getSubConfigurations();
    std::vector<util::Config> filtered_vars;
    std::vector<std::string> var_names;
    std::map<std::string, util::Metadata> var_metadata;
    for (auto var : zrgp_in_vars) {
        std::string name = to_upper(var.getString("name"));
        auto metadata = HDF5Reader{data_file, name}.read_metadata();
        auto size = metadata.getLong("size");

        std::vector<std::string> dims;
        var.get("dim", dims);
        std::string dim = dims[0];
        for (int i = 1; i < dims.size(); ++i) {
            dim += "*" + dims[i];
        }
        dims = {dim};
        dims.emplace(dims.begin(), "nblk");
        dims.emplace_back("nproma");
        Log::info() << "    " << name << dims;
        if (size == 0) {
            Log::info() << "   (empty --> ignored)" << std::endl;
            continue;
        }
        Log::info() << std::endl;

        var_names.emplace_back(name);
        var_metadata.emplace(name, metadata);
        filtered_vars.emplace_back(var);
    }
    util::Config zrgp_in_filtered = zrgp_in;
    zrgp_in_filtered.set("variables", filtered_vars);
    return std::make_tuple(zrgp_in_filtered, var_names, var_metadata);
}

class IFS {
    public:
    IFS(const util::Config& zrgp_in, const std::string& data_file, bool multifield) {
        util::Config multifield_config; // A copy of zrgp_in but only with the non-empty variables, and with some additional parameters set (e.g. nproma, nlev, etc.) based on the metadata of the variables
        std::tie(multifield_config, var_names, var_metadata) = extract_variables(zrgp_in, data_file);
        if (var_names.empty()) {
            ATLAS_THROW_EXCEPTION("No non-empty variables found in configuration");
        }
        DataType datatype(var_metadata[var_names[0]].getString("datatype"));
        IFS_grid = Grid(var_metadata[var_names[0]].getString("grid"));
        int ngptotg = IFS_grid.size();
        int nproma = var_metadata[var_names[0]].getLong("nproma");
        int nblk   = var_metadata[var_names[0]].getLong("nblk");
        IFS_blocked_fs = functionspace::BlockStructuredColumns(IFS_grid, util::Config("nproma",nproma));

        auto get_nfld = [&](const std::string& name) -> long { return (var_metadata.find(name) == var_metadata.end()) ? -1 : var_metadata[name].getLong("nfld");};
        int nlev = get_nfld("ISWA");
        int nrftotal_radgrid = get_nfld("IPERT");
        int nlwemiss = get_nfld("IEMISS");
        int nsw = get_nfld("IALD");
        int nprogaer = (get_nfld("IPROGAERO") >= 0 ) ? get_nfld("IPROGAERO") / nlev : -1 ;
        Log::info() << "Extracted variable information:" << std::endl;
        Log::info() << "    ngptotg  : " << ngptotg << "    (grid=" << IFS_grid.name() << ")" << std::endl;
        Log::info() << "    nproma   : " << nproma << std::endl;
        Log::info() << "    nblk     : " << nblk << std::endl;
        Log::info() << "    nlev     : " << nlev << std::endl;
        Log::info() << "    nsw      : " << nsw << std::endl;
        Log::info() << "    nlwemiss : " << nlwemiss << std::endl;
        Log::info() << "    nrftotal_radgrid   : " << nrftotal_radgrid << std::endl;
        Log::info() << "    nprogaer           : " << nprogaer << std::endl;
        
        multifield_config.set("type","MultiFieldCreatorRad");
        multifield_config.set("nproma",nproma);
        multifield_config.set("nlev",nlev);
        multifield_config.set("nrftotal_radgrid",nrftotal_radgrid);
        multifield_config.set("nlwemiss",nlwemiss);
        multifield_config.set("nsw",nsw);
        multifield_config.set("nprogaer", nprogaer);
        multifield_config.set("ngptot",IFS_blocked_fs.size());

        IFS_blocked_fs = functionspace::BlockStructuredColumns(IFS_grid, util::Config("nproma",nproma));

        Log::info() << "Creating IFS fields (multifield=" << std::boolalpha << multifield << "), read data on one rank, and scatter:" << std::endl;
        ATLAS_TRACE("Reading data");
        if (multifield) {
            zrgp_fields = field::MultiField(multifield_config);
            for( auto& f : zrgp_fields) {
                f.set_functionspace(IFS_blocked_fs);
                read_field_with_timing(data_file, f);
            }
        }
        else {
            auto nblks = IFS_blocked_fs.nblks();
            for( const auto& name : var_names ) {
                auto dim = get_nfld(name);
                Field f;
                if (dim == 1) {
                    f = zrgp_fields.add(Field(name, datatype, array::make_shape(nblks, nproma)));
                }
                else {
                    f = zrgp_fields.add(Field(name, datatype, array::make_shape(nblks, dim, nproma)));
                    f.set_levels(dim);
                }
                f.set_functionspace(IFS_blocked_fs);
                read_field_with_timing(data_file, f);
            }
        }
        nproma_ = nproma;
        datatype_ = datatype;
    }

    Grid grid() const {
        return IFS_grid;
    }
    FieldSet fields() const {
        return zrgp_fields;
    }
    functionspace::BlockStructuredColumns functionspace() const {
        return IFS_blocked_fs;
    }
    int nproma() const {
        return nproma_;
    }
    DataType datatype() const { return datatype_; }
    std::vector<std::string> var_names;
    std::map<std::string,util::Metadata> var_metadata;
    FieldSet zrgp_fields;
    Grid IFS_grid;
    functionspace::BlockStructuredColumns IFS_blocked_fs;
    int nproma_;
    DataType datatype_{DataType::KIND_REAL64};
};

class Radiation {
public:
    Radiation(Grid grid, long nproma, functionspace::BlockStructuredColumns _IFS_blocked_fs) : grid_(grid), nproma_(nproma), IFS_blocked_fs(_IFS_blocked_fs) {
        if (grid_ == IFS_blocked_fs.grid() && nproma_ == IFS_blocked_fs.nproma()) {
            rad_blocked_fs = IFS_blocked_fs;
        }
        else {
            rad_blocked_fs = functionspace::BlockStructuredColumns(grid_, grid::MatchingPartitioner(IFS_blocked_fs), util::Config("nproma",nproma_));
        }
    }
    Grid grid() const { return grid_; }
    int nproma() const { return nproma_; }
    functionspace::BlockStructuredColumns& functionspace() { return rad_blocked_fs; }
    void interpolate_from_IFS(const std::string& interpolation_method, bool on_device, const FieldSet& ifs_fields, FieldSet& rad_fields) {
        if (grid() == IFS_blocked_fs.grid()) {
            if (nproma() == IFS_blocked_fs.nproma()) {
                return; // perhaps assert that rad_fields equals ifs_fields still!, otherwise clone
            }
            ATLAS_TRACE("copy_blocked_to_blocked");
            copy_blocked_to_blocked(ifs_fields, rad_fields, on_device);
        }
        else {
            if (interpolation_from_IFS.count(interpolation_method) == 0) {
                ATLAS_TRACE("create_interpolation");
                interpolation_from_IFS.emplace(interpolation_method, create_interpolation(interpolation_method, on_device, IFS_blocked_fs, rad_blocked_fs));
            }
            ATLAS_TRACE("interpolate from IFS");
            interpolate(interpolation_from_IFS[interpolation_method], on_device, ifs_fields, rad_fields);
        }
    }
    void interpolate_to_IFS(const std::string& interpolation_method, bool on_device, const FieldSet& rad_fields, FieldSet& ifs_fields) {
        if (grid() == IFS_blocked_fs.grid()) {
            if (nproma() == IFS_blocked_fs.nproma()) {
                return; // perhaps assert that rad_fields equals ifs_fields still!, otherwise clone
            }
            ATLAS_TRACE("copy_blocked_to_blocked");
            copy_blocked_to_blocked(rad_fields, ifs_fields, on_device);
        }
        else {
            if (interpolation_to_IFS.count(interpolation_method) == 0) {
                ATLAS_TRACE("create_interpolation");
                interpolation_to_IFS.emplace(interpolation_method, create_interpolation(interpolation_method, on_device, rad_blocked_fs, IFS_blocked_fs));
            }
            ATLAS_TRACE("interpolate to IFS");
            interpolate(interpolation_to_IFS[interpolation_method], on_device, rad_fields, ifs_fields);
        }
    }
private:
    Grid grid_;
    long nproma_;
    functionspace::BlockStructuredColumns IFS_blocked_fs;
    functionspace::BlockStructuredColumns rad_blocked_fs;
    std::map<std::string, Interpolation> interpolation_from_IFS;
    std::map<std::string, Interpolation> interpolation_to_IFS;
};

int Program::execute(const AtlasTool::Args& args) {
    ATLAS_TRACE("main");

    IFS ifs(
        util::Config(args.getString("config-file")).getSubConfigurations()[0],
        args.getString("data-file"),
        args.getBool("multifield",false));

    FieldSet ifs_input_fields = ifs.fields();

    bool with_gmsh_output = args.getBool("output.gmsh",false);
    if (with_gmsh_output) {
        output_gmsh(ifs_input_fields, "ifs_input.msh", get_subconfiguration(args, "gmsh"));
    }

    Radiation rad{Grid{args.getString("rad.grid", ifs.grid().name())}, args.getLong("rad.nproma", ifs.nproma()), ifs.functionspace()};
    std::string interpolation_method = args.getString("interpolation","bicubic");
    bool on_device = args.getBool("on-device",false);
    int niter = args.getInt("niter", 1);

    if (on_device) {
        ifs.fields().updateDevice();
    }

    FieldSet ifs_output_fields = ifs.fields(); // to be created here!

    for ( int n=0; n<niter; ++n) {
        ATLAS_TRACE("iteration");
        Log::info() << "iteration " << n+1 << "/" << niter << std::endl;
        FieldSet rad_input_fields;
        FieldSet rad_output_fields;
        {
            if (rad.grid() == ifs.grid() && rad.nproma() == ifs.nproma()) {
                rad_input_fields  = ifs_input_fields;
            }
            else {
                for(auto ifs_input_field : ifs_input_fields) {
                    rad_input_fields.add(rad.functionspace().createField(ifs_input_field));
                }
                for(auto ifs_output_field : ifs_output_fields) {
                    rad_output_fields.add(rad.functionspace().createField(ifs_output_field));
                }
                if (on_device) {
                    rad_input_fields.allocateDevice();
                    rad_output_fields.allocateDevice();
                }
            }
        }
        {
            ATLAS_TRACE("remap IFS -> rad");
            rad.interpolate_from_IFS(interpolation_method, on_device, ifs_input_fields, rad_input_fields);
        }
        {
            ATLAS_TRACE("CALL radiation(rad_input_fields, rad_output_fields)");

            if (with_gmsh_output && n == 0) {
                output_gmsh(rad_input_fields, "rad_input.msh", get_subconfiguration(args,"gmsh"));
            }
            rad_output_fields = rad_input_fields;
            // Add content to rad_output_fields
        }
        {
            ATLAS_TRACE("remap rad -> IFS");
            rad.interpolate_to_IFS(interpolation_method, on_device, rad_output_fields, ifs_output_fields);
        }
    }
    if (with_gmsh_output) {
        output_gmsh(ifs_output_fields, "ifs_output.msh", get_subconfiguration(args,"gmsh"));
    }

    return success();
}

//-----------------------------------------------------------------------------

} // namespace


int main(int argc, char** argv) {
    atlas::Program tool(argc, argv);
    return tool.start();
}
