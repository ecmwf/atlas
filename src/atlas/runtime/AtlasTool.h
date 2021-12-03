/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#pragma once

#include <iostream>
#include <vector>

#include "eckit/option/CmdArgs.h"
#include "eckit/option/Separator.h"
#include "eckit/option/SimpleOption.h"
#include "eckit/option/VectorOption.h"
#include "eckit/runtime/Tool.h"

#include "atlas/library/Library.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/runtime/Log.h"
#include "atlas/util/Config.h"

//--------------------------------------------------------------------------------

using atlas::util::Config;
using eckit::option::CmdArgs;
using eckit::option::Option;
using eckit::option::Separator;
using eckit::option::SimpleOption;
using eckit::option::VectorOption;

namespace atlas {

class AtlasTool : public eckit::Tool {
public:
    typedef std::vector<eckit::option::Option*> Options;
    typedef eckit::option::CmdArgs Args;

protected:
    virtual bool serial() { return false; }
    virtual std::string indent() { return "      "; }
    virtual std::string briefDescription() { return ""; }
    virtual std::string longDescription() { return ""; }
    virtual std::string usage() { return name() + " [OPTION]... [--help,-h] [--debug]"; }

    void add_option(eckit::option::Option* option);

    virtual void help(std::ostream& out = Log::info());

    virtual int numberOfPositionalArguments() { return -1; }
    virtual int minimumPositionalArguments() { return 0; }

    bool handle_help();

    std::string get_positional_arg(const Args&, size_t pos) const;

public:
    AtlasTool(int argc, char** argv);

    int start();

    virtual void run();  // unused

    virtual int execute(const Args&) = 0;

    static constexpr int success() { return 0; }
    static constexpr int failed() { return 1; }

private:
    void setupLogging();

private:
    Options options_;
};

}  // namespace atlas
