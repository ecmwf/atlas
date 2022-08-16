/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <iostream>
#include <string>
#include <vector>

#include "eckit/option/CmdArgs.h"
#include "eckit/option/Separator.h"
#include "eckit/option/SimpleOption.h"
#include "eckit/option/VectorOption.h"
#include "eckit/runtime/Tool.h"

#include "eckit/filesystem/PathName.h"

#include "atlas_io/Exceptions.h"
#include "atlas_io/RecordPrinter.h"
#include "atlas_io/print/Bytes.h"

//--------------------------------------------------------------------------------

using eckit::Log;

class AtlasIOTool : public eckit::Tool {
public:
    using Options = std::vector<eckit::option::Option*>;
    using Args    = eckit::option::CmdArgs;

protected:
    virtual std::string indent() { return "      "; }
    virtual std::string briefDescription() { return ""; }
    virtual std::string longDescription() { return ""; }
    virtual std::string usage() { return name() + " [OPTION]... [--help,-h]"; }

    void add_option(eckit::option::Option* option) { options_.push_back(option); }

    virtual void help(std::ostream& out = Log::info()) {
        auto indented = [&](const std::string& s) -> std::string {
            std::string str = indent() + s;
            size_t pos      = 0;
            while ((pos = str.find('\n', pos)) != std::string::npos) {
                str.replace(pos, 1, '\n' + indent());
                ++pos;
            }
            return str;
        };

        out << "NAME\n" << indented(name());
        std::string brief = briefDescription();
        if (brief.size()) {
            out << " - " << brief << '\n';
        }

        std::string usg = usage();
        if (usg.size()) {
            out << '\n';
            out << "SYNOPSIS\n" << indented(usg) << '\n';
        }
        std::string desc = longDescription();
        if (desc.size()) {
            out << '\n';
            out << "DESCRIPTION\n" << indented(desc) << '\n';
        }
        out << '\n';
        out << "OPTIONS\n";
        for (Options::const_iterator it = options_.begin(); it != options_.end(); ++it) {
            std::stringstream s;
            s << **it;
            out << indented(s.str()) << "\n\n";
        }
        out << std::flush;
    }

    virtual int numberOfPositionalArguments() { return -1; }
    virtual int minimumPositionalArguments() { return 0; }

    bool handle_help() {
        for (int i = 1; i < argc(); ++i) {
            if (argv(i) == "--help" || argv(i) == "-h") {
                help(std::cout);
                return true;
            }
        }
        return false;
    }

public:
    AtlasIOTool(int argc, char** argv): eckit::Tool(argc, argv) {
        add_option(new eckit::option::SimpleOption<bool>("help", "Print this help"));
    }

    int start() {
        try {
            if (handle_help()) {
                return success();
            }

            if (argc() - 1 < minimumPositionalArguments()) {
                Log::error() << "Usage: " << usage() << std::endl;
                return failed();
            }

            Options opts                                  = options_;
            std::function<void(const std::string&)> dummy = [](const std::string&) {};
            Args args(dummy, opts, numberOfPositionalArguments(), minimumPositionalArguments() > 0);

            int err_code = execute(args);
            return err_code;
        }
        catch (eckit::Exception& e) {
            Log::error() << "** " << e.what() << " Caught in " << Here() << std::endl;
            Log::error() << "** Exception terminates " << name() << std::endl;
        }
        catch (std::exception& e) {
            Log::error() << "** " << e.what() << " Caught in " << Here() << std::endl;
            Log::error() << "** Exception terminates " << name() << std::endl;
        }
        return failed();
    }

    void run() final {}  // unused

    virtual int execute(const Args&) = 0;

    static constexpr int success() { return 0; }
    static constexpr int failed() { return 1; }

private:
    Options options_;
};


//----------------------------------------------------------------------------------------------------------------------

struct AtlasIOList : public AtlasIOTool {
    std::string briefDescription() override { return "Inspection of atlas-io files"; }
    std::string usage() override { return name() + " <file> [OPTION]... [--help,-h]"; }
    std::string longDescription() override {
        return "Inspection of atlas-io files\n"
               "\n"
               "       <file>: path to atlas-io file";
    }

    AtlasIOList(int argc, char** argv): AtlasIOTool(argc, argv) {
        add_option(new eckit::option::SimpleOption<std::string>("format", "Output format"));
        add_option(new eckit::option::SimpleOption<bool>("version", "Print version of records"));
        add_option(new eckit::option::SimpleOption<bool>("details", "Print detailed information"));
    }
    int execute(const Args& args) override {
        auto return_code = success();

        using namespace atlas;

        // User sanity checks
        if (args.count() == 0) {
            Log::error() << "No file specified." << std::endl;
            help(Log::error());
            return failed();
        }

        // Configuration
        eckit::LocalConfiguration config;
        config.set("format", args.getString("format", "table"));
        config.set("details", args.getBool("details", false));

        // Loop over files
        for (size_t f = 0; f < args.count(); ++f) {
            eckit::PathName file(args(f));
            if (!file.exists()) {
                Log::error() << "File does not exist: " << file << std::endl;
                return failed();
            }
            auto filesize = size_t(file.size());

            io::Session session;

            std::uint64_t pos = 0;
            try {
                while (pos < filesize) {
                    auto uri    = io::Record::URI{file, pos};
                    auto record = io::RecordPrinter{uri, config};

                    std::stringstream out;
                    out << "\n# " << uri.path << " [" << uri.offset << "]    "
                        << "{ size: " << atlas::io::Bytes{record.size()}.str(0) << ",    version: " << record.version()
                        << ",    created: " << record.time() << " }";
                    out << '\n' << (config.getString("format") == "table" ? "" : "---") << '\n';
                    out << record << std::endl;

                    std::cout << out.str();

                    pos += record.size();
                }
            }
            catch (const io::Exception& e) {
                Log::error() << "    ATLAS-IO-ERROR: " << e.what() << std::endl;
                return_code = failed();
            }
        }
        return return_code;
    }
};

//------------------------------------------------------------------------------------------------------

int main(int argc, char** argv) {
    return AtlasIOList{argc, argv}.start();
}
