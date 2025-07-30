/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <algorithm>
#include <chrono>
#include <thread>

#include "eckit/config/LibEcKit.h"
#include "eckit/filesystem/LocalPathName.h"
#include "eckit/log/FileTarget.h"
#include "eckit/log/PrefixTarget.h"

#include "atlas/library/config.h"
#include "atlas/runtime/AtlasTool.h"

namespace atlas {
static void usage(const std::string& name) {
    Log::info() << "dummy usage" << std::endl;
}
}  // namespace atlas

namespace {

int digits(int number) {
    int d = 0;
    while (number) {
        number /= 10;
        d++;
    }
    return d;
}

static std::string debug_prefix(const std::string& libname) {
    std::string s = libname;
    std::transform(s.begin(), s.end(), s.begin(), ::toupper);
    s += "_DEBUG";
    return s;
}

void debug_addTarget(eckit::LogTarget* target) {
    for (std::string libname : eckit::system::Library::list()) {
        const eckit::system::Library& lib = eckit::system::Library::lookup(libname);
        if (lib.debug()) {
            lib.debugChannel().addTarget(new eckit::PrefixTarget(debug_prefix(libname), target));
        }
    }
    if (eckit::Log::debug()) {
        eckit::Log::debug().addTarget(target);
    }
}

void debug_setTarget(eckit::LogTarget* target) {
    for (std::string libname : eckit::system::Library::list()) {
        const eckit::system::Library& lib = eckit::system::Library::lookup(libname);
        if (lib.debug()) {
            lib.debugChannel().setTarget(new eckit::PrefixTarget(debug_prefix(libname), target));
        }
    }
    if (eckit::Log::debug()) {
        eckit::Log::debug().setTarget(target);
    }
}

void debug_reset() {
    for (std::string libname : eckit::system::Library::list()) {
        const eckit::system::Library& lib = eckit::system::Library::lookup(libname);
        if (lib.debug()) {
            lib.debugChannel().reset();
        }
    }
    if (eckit::Log::debug()) {
        eckit::Log::debug().reset();
    }
}

bool getEnv(const std::string& env, bool default_value) {
    const char* cenv = ::getenv(env.c_str());
    if (cenv != nullptr) {
        return eckit::Translator<std::string, bool>()(cenv);
    }
    return default_value;
}

int getEnv(const std::string& env, int default_value) {
    const char* cenv = ::getenv(env.c_str());
    if (cenv != nullptr) {
        return eckit::Translator<std::string, int>()(cenv);
    }
    return default_value;
}

std::string getEnv(const std::string& env, const std::string& default_value) {
    const char* cenv = ::getenv(env.c_str());
    if (cenv != nullptr) {
        return cenv;
    }
    return default_value;
}

void setEnv(const std::string& env, bool value) {
    constexpr int DO_NOT_REPLACE_IF_EXISTS = 0;
    ::setenv(env.c_str(), eckit::Translator<bool, std::string>()(value).c_str(), DO_NOT_REPLACE_IF_EXISTS);
}


}  // namespace

namespace atlas {

static bool use_logfile;
static std::string logfile_name;
static std::string workdir;

[[noreturn]] void atlas_terminate() {
    // This routine is called for uncaught exceptions.
    // It can be set with std::set_terminate( &atlas_terminate )

    Log::flush();

    if (not use_logfile and mpi::size() > 1) {
        eckit::LogTarget* logfile = new eckit::FileTarget(logfile_name);
        Log::error().addTarget(logfile);
    }
    if (std::exception_ptr eptr = std::current_exception()) {
        std::ostream& out = Log::error();
        try {
            std::rethrow_exception(eptr);  // throw to recognise the type
        }
        catch (const eckit::Abort& exception) {
            out << "\n"
                << "=========================================\n"
                << "[" << mpi::rank() << "] Aborting " << eckit::Main::instance().displayName() << "\n"
                << "-----------------------------------------\n"
                << exception.what() << "\n";
            if (exception.location()) {
                out << "-----------------------------------------\n"
                    << "LOCATION: " << exception.location() << "\n";
            }
            out << "-----------------------------------------\n"
                << "BACKTRACE\n"
                << "-----------------------------------------\n"
                << exception.callStack() << "\n"
                << "=========================================\n"
                << std::endl;
        }
        catch (const eckit::Exception& exception) {
            out << "\n"
                << "=========================================\n"
                << "[" << mpi::rank() << "] TERMINATING " << eckit::Main::instance().displayName() << "\n"
                << "-----------------------------------------\n"
                << exception.what() << "\n"
                << "-----------------------------------------\n";

            if (exception.location()) {
                out << "LOCATION: " << exception.location() << "\n"
                    << "-----------------------------------------\n";
            }

            out << "BACKTRACE\n"
                << "-----------------------------------------\n"
                << exception.callStack() << "\n"
                << "=========================================\n"
                << std::endl;
        }
        catch (const std::exception& exception) {
            out << "\n"
                << "=========================================\n"
                << "[" << mpi::rank() << "] TERMINATING " << eckit::Main::instance().displayName() << "\n"
                << "-----------------------------------------\n"
                << exception.what() << "\n"
                << "-----------------------------------------\n"
                << "BACKTRACE\n"
                << "-----------------------------------------\n"
                << backtrace() << "\n"
                << "=========================================\n"
                << std::endl;
        }
        catch (...) {
            out << "\n"
                << "=========================================\n"
                << "[" << mpi::rank() << "] TERMINATING " << eckit::Main::instance().displayName() << "\n"
                << "-----------------------------------------\n"
                << "BACKTRACE\n"
                << "-----------------------------------------\n"
                << backtrace() << "\n"
                << "=========================================" << std::endl;
        }
    }

    eckit::LibEcKit::instance().abort();

    // Just in case we end up here, as last resort, exit immediately without
    // cleanup.
    std::_Exit(EXIT_FAILURE);
}

}  // namespace atlas


void atlas::AtlasTool::add_option(eckit::option::Option* option) {
    options_.push_back(option);
}

void atlas::AtlasTool::help(std::ostream& out) {
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

bool atlas::AtlasTool::handle_help() {
    for (int i = 1; i < argc(); ++i) {
        if (argv(i) == "--help" || argv(i) == "-h") {
            if (taskID() == 0) {
                help(std::cout);
            }
            return true;
        }
    }
    return false;
}

std::string atlas::AtlasTool::get_positional_arg(const Args& args, size_t p) const {
    int c = 0;
    for (int i = 0; i < args.count(); ++i) {
        if (args(i)[0] == '-') {
            i++;
        }
        else {
            if (c == p) {
                return args(i);
            }
            c++;
        }
    }
    return std::string{};
}

atlas::AtlasTool::AtlasTool(int argc, char** argv): eckit::Tool(argc, argv) {
    eckit::LibEcKit::instance().setAbortHandler([] {
        std::cerr << "[" << atlas::mpi::rank() << "] "
                  << "calling MPI_Abort";
        if (not use_logfile and mpi::size() > 1) {
            std::cerr << ", logfile: " << logfile_name;
        }
        std::cerr << std::endl;
        std::this_thread::sleep_for(std::chrono::milliseconds(3000));
        atlas::mpi::comm().abort(1);
    });
    std::set_terminate(&atlas_terminate);
    setEnv("ECKIT_EXCEPTION_IS_SILENT", true);
    setEnv("ECKIT_ASSERT_FAILED_IS_SILENT", true);
    setEnv("ATLAS_FPE", true);
    setEnv("ATLAS_SIGNAL_HANDLER", true);
    add_option(new SimpleOption<bool>("help", "Print this help"));
    add_option(new SimpleOption<long>("debug", "Debug level"));
    taskID(eckit::mpi::comm("world").rank());
    workdir = getEnv("ATLAS_WORKDIR", eckit::PathName(eckit::LocalPathName::cwd()).fullName());

    add_option(new SimpleOption<std::string>("display-name", "Name to give to program"));
}

int atlas::AtlasTool::start() {
    try {
        if (handle_help()) {
            return success();
        }

        if (argc() - 1 < minimumPositionalArguments()) {
            if (taskID() == 0) {
                std::cout << "Usage: " << usage() << std::endl;
            }
            return failed();
        }
        atlas::initialize();
        setupLogging();

        Options opts = options_;
        Args args(&atlas::usage, opts, numberOfPositionalArguments(), minimumPositionalArguments() > 0);

        int err_code = execute(args);
        atlas::finalize();
        atlas::mpi::finalize();
        return err_code;
    }
    catch (...) {
        atlas_terminate();
    }
}

void atlas::AtlasTool::run() {}

void atlas::AtlasTool::setupLogging() {
    int log_rank = getEnv("ATLAS_LOG_RANK", 0);
    use_logfile  = getEnv("ATLAS_LOG_FILE", false);

    int d               = digits(mpi::size());
    std::string rankstr = std::to_string(taskID());
    for (int i = rankstr.size(); i < d; ++i) {
        rankstr = "0" + rankstr;
    }

    logfile_name = workdir + "/" + displayName() + ".log.p" + rankstr;

    if (use_logfile) {
        eckit::LogTarget* logfile = new eckit::FileTarget(logfile_name);

        if (int(mpi::rank()) == log_rank) {
            if (Log::info()) {
                Log::info().addTarget(logfile);
            }
            if (Log::warning()) {
                Log::warning().addTarget(logfile);
            }
            if (Log::error()) {
                Log::error().addTarget(logfile);
            }
            debug_addTarget(logfile);
        }
        else {
            if (Log::info()) {
                Log::info().setTarget(logfile);
            }
            if (Log::warning()) {
                Log::warning().setTarget(logfile);
            }
            if (Log::error()) {
                Log::error().setTarget(logfile);
            }
            debug_setTarget(logfile);
        }
    }
    else {
        if (int(mpi::rank()) != log_rank) {
            if (Log::info()) {
                Log::info().reset();
            }
            if (Log::warning()) {
                Log::warning().reset();
            }
            if (Log::error()) {
                Log::error().reset();
            }
            debug_reset();
        }
    }
    // codechecker_false_positive [NewDeleteLeaks] Potential leak of memory pointed to by 'logfile'
}
