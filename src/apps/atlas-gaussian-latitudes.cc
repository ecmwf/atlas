/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <cassert>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <vector>

#include "eckit/config/Resource.h"
#include "eckit/runtime/Main.h"
#include "eckit/runtime/Tool.h"

#include "atlas/grid/detail/spacing/GaussianSpacing.h"
#include "atlas/library.h"
#include "atlas/runtime/Exception.h"

//------------------------------------------------------------------------------------------------------

using namespace eckit;
using namespace atlas;
using atlas::grid::spacing::GaussianSpacing;

//------------------------------------------------------------------------------------------------------

class AtlasGaussianLatitudes : public eckit::Tool {
    void run() override;

public:
    AtlasGaussianLatitudes(int argc, char** argv): eckit::Tool(argc, argv) {
        atlas::initialise(argc,argv);

        do_run = false;

        bool help = Resource<bool>("--help", false);

        std::string help_str =
            "NAME\n"
            "       atlas-gaussian-latitudes - Compute gaussian latitudes, given "
            "the N number\n"
            "\n"
            "SYNOPSIS\n"
            "       atlas-gaussian-latitudes [--help] [-N] [OPTION]...\n"
            "\n"
            "DESCRIPTION\n"
            "       Compute gaussian latitudes, given the N number.\n"
            "       Latitudes start at the pole (+90), and decrease in value.\n"
            "\n"
            "       -N          Number of points between pole and equator\n"
            "\n"
            "       --full      If set, all latitudes will be given, otherwise "
            "only\n"
            "                   between North pole and equator.\n"
            "\n"
            "       --format    \"table\" (default)\n"
            "                   \"C\"\n"
            "\n"
            "       --compact   Write 5 latitudes per line if the format supports "
            "it\n"
            "\n"
            "AUTHOR\n"
            "       Written by Willem Deconinck.\n"
            "\n"
            "ECMWF                        December 2014";

        N = Resource<int>("-N", 0);

        full = Resource<bool>("--full", false);

        compact = Resource<bool>("--compact", false);

        format = Resource<std::string>("--format", std::string("table"));

        if (N > 0) {
            do_run = true;
        }

        if (do_run == false) {
            if (help) {
                Log::info() << help_str << std::endl;
            }
            else {
                Log::info() << "usage: atlas-gaussian-latitudes [--help] [-N] [OPTION]..." << std::endl;
            }
        }
    }

private:
    int N;
    bool full;
    bool compact;
    std::string format;
    bool do_run;
};

//------------------------------------------------------------------------------

void AtlasGaussianLatitudes::run() {
    if (!do_run) {
        return;
    }

    GaussianSpacing lats(2 * N);

    int end = full ? 2 * N : N;

    if (format == "table") {
        for (int jlat = 0; jlat < end; ++jlat) {
            std::cout << std::setw(4) << jlat + 1 << std::setw(17) << std::setprecision(12) << std::fixed << lats[jlat]
                      << std::endl;
        }
    }
    if (format == "C") {
        std::cout << "double lat[] = {" << std::endl;
        for (int jlat = 0; jlat < end; ++jlat) {
            std::cout << std::setw(16) << std::setprecision(12) << std::fixed << lats[jlat];
            if (jlat != end - 1) {
                std::cout << ",";
            }
            if ((compact && (jlat + 1) % 5 == 0) || !compact || jlat == end - 1) {
                std::cout << std::endl;
            }
        }
        std::cout << "};" << std::endl;
    }
    atlas::finalise();
}

//------------------------------------------------------------------------------

int main(int argc, char** argv) {
    AtlasGaussianLatitudes tool(argc, argv);
    return tool.start();
}
