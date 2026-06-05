/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "eckit/os/BackTrace.h"

#include "atlas/library/Library.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/runtime/Log.h"



// for MacOS backtrace
#ifdef __APPLE__
#include <execinfo.h>
#include <dlfcn.h>
#include <string>
#include <string_view>
#include <sstream>
#include <iomanip>
namespace atlas {
namespace {
inline std::string hex_str(void* addr) {
    std::stringstream ss;
    ss << "0x" << std::hex << reinterpret_cast<uintptr_t>(addr);
    return ss.str();
};
inline std::vector<std::string> atos_command_args(const Dl_info& info, void* runtime_addr) {
    // Arguments to "atos" (Apple's address-to-symbol utility) to resolve function names and
    // source locations, including line numbers, from the runtime addresses obtained from backtrace.
    return std::vector<std::string> {
        "-fullPath",                   // Print the full path of the source files
        "-i",                          // Display inlined symbols
        "-o", info.dli_fname,          // Path to the binary image containing the address
        "-l", hex_str(info.dli_fbase), // Load address of the binary image
        hex_str(runtime_addr)};        // The address to be resolved
}
inline std::string full_command(const std::string& cmd, const std::vector<std::string>& args) {
    std::stringstream ss;
    ss << cmd;
    for (const auto& arg : args) {
        ss << " " << arg;
    }
    return ss.str();
}
inline std::string run_command(const std::string& cmd, const std::vector<std::string>& args) {
    // This may not be all that safe in a multi-threaded context, or a signal handling context
    std::string full_cmd = full_command(cmd, args);
    std::string output;
    FILE* fp = ::popen(full_cmd.c_str(), "r");
    if (fp) {
        char* line_buffer = nullptr;
        size_t line_len = 0;

        // getline dynamically allocates memory and handles ANY string size safely
        if (getline(&line_buffer, &line_len, fp) != -1) {
            output = std::string(line_buffer, line_len);
        }
        // POSIX requirement: Free the memory allocated by getline
        std::free(line_buffer);
        ::pclose(fp);
    }
    return output;
}

// This function is a custom backtrace implementation for MacOS that uses the atos command
// (Apple's address-to-symbol utility) to resolve function names and source locations,
// including line numbers, from the raw addresses obtained from backtrace.
// It also filters out frames from Apple internal system libraries for better readability.
std::string macos_backtrace() {
    void* callstack[128];
    int frames = ::backtrace(callstack, 128);

    std::stringstream out;
    int i = 0;
    for (int f = 0; f < frames; ++f) {
        Dl_info info;

        // Match the raw pointer address to its underlying loaded image binary
        if (dladdr(callstack[f], &info) && info.dli_fname) {
            // SKIP frame processing if it belongs to Apple internal system libraries
            if (std::string_view{info.dli_fname}.find("/usr/lib/") != std::string_view::npos ||
                std::string_view{info.dli_fname}.find("/System/Library/") != std::string_view::npos) {
                continue;
            }
            auto cmd_output = run_command("atos", atos_command_args(info, callstack[f]));
            if (cmd_output.empty()) {
                out << "[" << i++ << "]\t" << "[FRAME PROCESSING FAILED]: " << full_command("atos", atos_command_args(info, callstack[f])) << '\n';
            }
            else {
                bool in_this_file = cmd_output.find("atlas/runtime/Log.cc") != std::string_view::npos;
                if (not in_this_file) {
                    out << "[" << i++ << "]\t" << cmd_output;
                }
            }
        }
        else {
            out << "[" << i++ << "]\tUnknown Address: " << hex_str(callstack[f]) << '\n';
        }
    }
    return out.str();
}
}
}
#endif

#if ATLAS_HAVE_FORTRAN
#include "fckit/Log.h"
#endif

namespace atlas {

std::string backtrace() {
#if defined(__APPLE__)
    return atlas::macos_backtrace();
#else
    return eckit::BackTrace::dump();
#endif
}

namespace detail {

void debug_parallel_here(const eckit::CodeLocation& here) {
    const auto& comm = mpi::comm();
    comm.barrier();
    Log::info() << "DEBUG_PARALLEL() @ " << here << std::endl;
}

void debug_parallel_what(const eckit::CodeLocation& here, const std::string& what) {
    const auto& comm = mpi::comm();
    comm.barrier();
    Log::info() << "DEBUG_PARALLEL(" << what << ") @ " << here << std::endl;
}

}  // namespace detail

Log::Channel& Log::info() {
    return atlas::Library::instance().infoChannel();
}

Log::Channel& Log::warning() {
    return atlas::Library::instance().warningChannel();
}

Log::Channel& Log::trace() {
    return atlas::Library::instance().traceChannel();
}

Log::Channel& Log::debug() {
    return atlas::Library::instance().debugChannel();
}

void Log::addFortranUnit(int unit, Style style, const char* prefix) {
#if ATLAS_HAVE_FORTRAN
    fckit::Log::addFortranUnit(unit, fckit::Log::Style(style), prefix);
#else
/*NOTIMP*/
#endif
}

void Log::setFortranUnit(int unit, Style style, const char* prefix) {
#if ATLAS_HAVE_FORTRAN
    fckit::Log::setFortranUnit(unit, fckit::Log::Style(style), prefix);
#else
/*NOTIMP*/
#endif
}

}  // namespace atlas
