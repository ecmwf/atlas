/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "atlas/library.h"
#include "atlas/runtime/Log.h"
#include "eckit/config/Resource.h"
#include "eckit/runtime/Tool.h"

using namespace eckit;
namespace atlas {

//----------------------------------------------------------------------------------------------------------------------

class Version : public Tool {
public:
    Version(int argc, char** argv): Tool(argc, argv) {}

    ~Version() override = default;

    void run() override;
};

void Version::run() {
    if (Resource<bool>("--version", false)) {
        Log::info() << atlas::Library::instance().version() << std::endl;
        return;
    }
    else if (Resource<bool>("--git", false)) {
        Log::info() << atlas::Library::instance().gitsha1(12) << std::endl;
        return;
    }
    else if (Resource<bool>("--info", false)) {
        Log::info() << atlas::Library::instance().information() << std::endl;
        return;
    }
    else if (Resource<bool>("--init", false)) {
        Log::info() << "+ atlas::initialize()" << std::endl;
        atlas::initialize();
        Log::info() << "+ atlas::finalize()" << std::endl;
        atlas::finalize();
        return;
    }
    else if (Resource<bool>("--help", false)) {
        Log::info() << "NAME\n"
                       "       atlas - Framework for parallel flexible data structures on "
                       "the sphere\n"
                       "\n"
                       "SYNOPSIS\n"
                       "       atlas [--help] [--version] [--git] [--info]\n"
                       "\n"
                       "DESCRIPTION\n"
                       "       Framework for parallel flexible data structures on the "
                       "sphere.\n"
                       "\n"
                       "OPTIONS\n"
                       "       --help\n"
                       "           Print this help\n"
                       "\n"
                       "       --version\n"
                       "           Print short version string 'MAJOR.MINOR.PATCH' \n"
                       "\n"
                       "       --info\n"
                       "           Print build configuration and features\n"
                       "\n"
                       "       --init\n"
                       "           Initialise and finalise atlas library, useful for printing debug information (environment ATLAS_DEBUG=1)\n"
                       "\n"
                       "AUTHOR\n"
                       "       Written by Willem Deconinck.\n"
                       "\n"
                       "ECMWF                        December 2014"
                    << std::endl;
        return;
    }
    else {
        Log::info() << "usage: atlas [--help] [--version] [--git] [--info] [--init]" << std::endl;
    }
}

}  // namespace atlas

//----------------------------------------------------------------------------------------------------------------------

using namespace atlas;

int main(int argc, char** argv) {
    Version tool(argc, argv);
    return tool.start();
}
