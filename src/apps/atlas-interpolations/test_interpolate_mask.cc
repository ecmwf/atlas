/*
 * (C) Copyright 2026 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "atlas/grid.h"
#include "atlas/io/atlas-io.h"
#include "atlas/library.h"

#include "MaskIO.h"

int main(int argc, char** argv) {
    atlas::initialize(argc, argv);
    if (argc != 4) {
        std::cerr << "Usage: " << argv[0] << " <write|check|write-field|check-extracted> <path> <grid>" << std::endl;
        atlas::finalize();
        return 1;
    }

    const std::string action = argv[1];
    const std::string path   = argv[2];
    const atlas::Grid grid(argv[3]);

    if (action == "write") {
        std::vector<int> mask(grid.size());
        for (size_t i = 0; i < mask.size(); ++i) {
            mask[i] = static_cast<int>(i % 2);
        }
        atlas::MaskIO::write(path, atlas::MaskData{grid.name(), mask});
    }
    else if (action == "write-field") {
        std::vector<double> values(grid.size(), 1.);
        values[1] = 9999.;
        atlas::io::RecordWriter writer;
        writer.set("grid.name", grid.name());
        writer.set("fields[0].array", atlas::io::ref(values));
        writer.write(path, atlas::io::Mode::write);
    }
    else if (action == "check-extracted") {
        const atlas::MaskData data = atlas::MaskIO::read(path);
        if (data.grid != grid.name() || data.values.size() != grid.size()) {
            std::cerr << "Unexpected extracted mask grid or size" << std::endl;
            atlas::finalize();
            return 1;
        }
        for (size_t i = 0; i < data.values.size(); ++i) {
            const int expected = i == 1 ? 0 : 1;
            if (data.values[i] != expected) {
                std::cerr << "Unexpected extracted mask value at index " << i << ": " << data.values[i]
                          << std::endl;
                atlas::finalize();
                return 1;
            }
        }
    }
    else if (action == "check") {
        std::string format;
        int version{0};
        std::string grid_spec;
        atlas::io::RecordReader reader(path);
        reader.read("format", format);
        reader.read("version", version);
        reader.read("grid", grid_spec);
        reader.wait();

        std::istringstream stream(grid_spec);
        const atlas::Grid encoded_grid{atlas::util::Config(stream)};
        if (format != "atlas-mask" || version != 1 || encoded_grid != grid) {
            std::cerr << "Unexpected mask metadata: format=" << format << ", version=" << version
                      << ", grid=" << grid_spec << std::endl;
            atlas::finalize();
            return 1;
        }

        std::vector<int> mask(grid.size());
        atlas::MaskIO::read(path, atlas::mdspan<int, atlas::dims<1>>(mask.data(), mask.size()));
        for (size_t i = 0; i < mask.size(); ++i) {
            if (mask[i] != static_cast<int>(i % 2)) {
                std::cerr << "Unexpected mask value at index " << i << ": " << mask[i] << std::endl;
                atlas::finalize();
                return 1;
            }
        }
        const std::string gmsh_path = path.substr(0, path.size() - std::string(".atlas").size()) + ".msh";
        if (not std::ifstream(gmsh_path).good()) {
            std::cerr << "Gmsh mask output was not created: " << gmsh_path << std::endl;
            atlas::finalize();
            return 1;
        }
    }
    else {
        std::cerr << "Unknown action: " << action << std::endl;
        atlas::finalize();
        return 1;
    }

    atlas::finalize();
    return 0;
}