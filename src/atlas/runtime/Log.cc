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
#include "eckit/config/Resource.h"

#include <chrono>
#include <thread>
#include <unistd.h>
#include "eckit/io/FDataSync.h"

#include "eckit/config/Resource.h"

#include "atlas/library/Library.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/runtime/Log.h"
#include "atlas/runtime/trace/StopWatch.h"
#include "atlas/runtime/Exception.h"



// for MacOS backtrace
#ifdef __APPLE__
#include <algorithm>
#include <execinfo.h>
#include <dlfcn.h>
#include <cstdlib>
#include <vector>
#include <string>
#include <string_view>
#include <sstream>
#include <iomanip>
#include <unordered_map>
#include <regex>
namespace atlas {
namespace {
inline std::string hex_str(void* addr) {
    std::stringstream ss;
    ss << "0x" << std::hex << reinterpret_cast<uintptr_t>(addr);
    return ss.str();
};
inline std::vector<std::string> atos_command_args(const Dl_info& info, const std::vector<void*>& runtime_addrs) {
    std::vector<std::string> args{
        "-fullPath",                   // Print the full path of the source files
        "-i",                          // Display inlined symbols
        "-o", info.dli_fname,          // Path to the binary image containing the addresses
        "-l", hex_str(info.dli_fbase)  // Load address of the binary image
    };
    args.reserve(args.size() + runtime_addrs.size());
    for (void* runtime_addr : runtime_addrs) {
        args.push_back(hex_str(runtime_addr));
    }
    return args;
}
inline std::string full_command(const std::string& cmd, const std::vector<std::string>& args) {
    std::stringstream ss;
    ss << cmd;
    for (const auto& arg : args) {
        ss << " " << arg;
    }
    return ss.str();
}
inline std::string trim_ascii_whitespace(std::string_view text) {
    const auto begin = text.find_first_not_of(" \t\r\n");
    if (begin == std::string_view::npos) {
        return {};
    }
    const auto end = text.find_last_not_of(" \t\r\n");
    return std::string(text.substr(begin, end - begin + 1));
}
inline std::vector<std::string> split_lines(std::string_view text) {
    std::vector<std::string> lines;
    size_t begin = 0;
    while (begin < text.size()) {
        const auto end = text.find('\n', begin);
        const auto line = trim_ascii_whitespace(text.substr(begin, end == std::string_view::npos ? text.size() - begin : end - begin));
        if (not line.empty()) {
            lines.push_back(line);
        }
        if (end == std::string_view::npos) {
            break;
        }
        begin = end + 1;
    }
    return lines;
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
        ssize_t chars_read = 0;
        while ((chars_read = getline(&line_buffer, &line_len, fp)) != -1) {
            const auto line = trim_ascii_whitespace(
                std::string_view(line_buffer, static_cast<size_t>(chars_read)));
            if (line.empty()) {
                continue;
            }
            if (not output.empty()) {
                output.push_back('\n');
            }
            output += line;
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

    auto simplify = [](std::string frame) {
        auto replace = [&frame] (const std::string& from, const std::string& to) {
            size_t start_pos = 0;
            while((start_pos = frame.find(from, start_pos)) != std::string::npos) {
                frame.replace(start_pos, from.length(), to);
                start_pos += to.length(); // Handles case where 'to' is a substring of 'from'
            }
            return frame;
        };
        auto simplify_std_string = [&replace]() {
            replace("std::basic_string<char, std::char_traits<char>, std::allocator<char>>", "std::string");
        };
        auto simplify_std_vector = [&frame]() {
            // Match std::vector with its type and remove the redundant allocator argument
            // Matches: std::vector<TYPE, std::allocator<TYPE>>
            // Captures TYPE in group 1
            std::regex allocatorRegex(R"(std::vector<([^,>]+),\s*std::allocator<\1\s*>>)");
            
            // Replace the entire match with just std::vector<TYPE>
            frame = std::regex_replace(frame, allocatorRegex, "std::vector<$1>");
        };
        replace("std::__1::", "std::");
        simplify_std_string();
        simplify_std_vector();
        return frame;
    };

    void* callstack[128];
    int frames = ::backtrace(callstack, 128);

    struct Frame {
        int frame_index;
        void* runtime_addr;
        Dl_info info;
    };
    struct FrameBatch {
        Dl_info info;
        std::vector<Frame> frames;
    };
    struct FrameBatchKey {
        std::string image_path;
        void* image_base;

        bool operator==(const FrameBatchKey& other) const {
            return image_base == other.image_base && image_path == other.image_path;
        }

        FrameBatchKey(const Frame& frame) {
            image_path = frame.info.dli_fname ? frame.info.dli_fname : "";
            image_base = frame.info.dli_fbase;
        }
    };
    struct FrameBatchKeyHash {
        size_t operator()(const FrameBatchKey& key) const {
            return std::hash<std::string>{}(key.image_path) ^ (std::hash<void*>{}(key.image_base) << 1);
        }
    };

    std::vector<Frame> pending_frames;
    pending_frames.reserve(frames);
    std::vector<std::string> formatted_frames(frames);
    std::vector<bool> emit_frame(frames, false);

    for (int f = 0; f < frames; ++f) {
        Dl_info info;

        // Match the raw pointer address to its underlying loaded image binary
        if (dladdr(callstack[f], &info) && info.dli_fname) {
            // SKIP frame processing if it belongs to Apple internal system libraries
            if (std::string_view{info.dli_fname}.find("/usr/lib/") != std::string_view::npos ||
                std::string_view{info.dli_fname}.find("/System/Library/") != std::string_view::npos) {
                continue;
            }
            pending_frames.push_back({f, callstack[f], info});
        }
        else {
            emit_frame[f] = true;
            formatted_frames[f] = "Unknown Address: " + hex_str(callstack[f]);
        }
    }

    std::vector<FrameBatch> batches;
    batches.reserve(pending_frames.size());
    std::unordered_map<FrameBatchKey, size_t, FrameBatchKeyHash> batch_index_by_image;
    batch_index_by_image.reserve(pending_frames.size());
    for (const auto& pending_frame : pending_frames) {
        auto [it, inserted] = batch_index_by_image.emplace(pending_frame, batches.size());
        if (inserted) {
            batches.push_back({pending_frame.info, {}});
        }
        batches[it->second].frames.push_back(pending_frame);
    }

    std::vector<std::string> resolved_frames(frames);
    std::vector<bool> frame_resolution_failed(frames, false);
    std::vector<std::string> failed_commands(frames);
    for (const auto& batch : batches) {
        std::vector<void*> runtime_addrs;
        runtime_addrs.reserve(batch.frames.size());
        for (const auto& pending_frame : batch.frames) {
            runtime_addrs.push_back(pending_frame.runtime_addr);
        }

        auto args = atos_command_args(batch.info, runtime_addrs);
        auto outputs = split_lines(run_command("atos", args));
        if (outputs.size() != batch.frames.size()) {
            const auto command = full_command("atos", args);
            for (const auto& pending_frame : batch.frames) {
                frame_resolution_failed[pending_frame.frame_index] = true;
                failed_commands[pending_frame.frame_index] = command;
            }
            continue;
        }

        for (size_t idx = 0; idx < batch.frames.size(); ++idx) {
            resolved_frames[batch.frames[idx].frame_index] = outputs[idx];
        }
    }

    for (const auto& pending_frame : pending_frames) {
        if (frame_resolution_failed[pending_frame.frame_index]) {
            emit_frame[pending_frame.frame_index] = true;
            formatted_frames[pending_frame.frame_index] = "[FRAME PROCESSING FAILED]: " + failed_commands[pending_frame.frame_index];
            continue;
        }

        const auto& cmd_output = resolved_frames[pending_frame.frame_index];
        if (cmd_output.empty()) {
            emit_frame[pending_frame.frame_index] = true;
            formatted_frames[pending_frame.frame_index] =
                "[FRAME PROCESSING FAILED]: " + full_command("atos", atos_command_args(pending_frame.info, std::vector<void*>{pending_frame.runtime_addr}));
            continue;
        }

        bool in_this_file = cmd_output.find("atlas/runtime/Log.cc") != std::string_view::npos;
        if (not in_this_file) {
            emit_frame[pending_frame.frame_index] = true;
            formatted_frames[pending_frame.frame_index] = simplify(cmd_output);
        }
    }

    int nb_emitted_frames = std::count(emit_frame.begin(), emit_frame.end(), true);
    std::stringstream out;
    for (int i = 0, f = 0; f < frames; ++f) {
        if (emit_frame[f]) {
            std::string frame_index_str = "#" + std::to_string(i++);
            out << std::left << std::setw(6) << frame_index_str << formatted_frames[f];
            if (i < nb_emitted_frames) {
                out << '\n';
            }
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

bool ATLAS_BACKTRACE_RESOLVE() {
    static bool resolve = eckit::Resource<bool>("$ATLAS_BACKTRACE_RESOLVE", true);
    return resolve;
}

std::string backtrace() {
#if defined(__APPLE__)
    if (ATLAS_BACKTRACE_RESOLVE()) {
        return atlas::macos_backtrace();
    }
#endif
    return eckit::BackTrace::dump();
}

namespace detail {

namespace {

double ATLAS_MPI_BARRIER_TIMEOUT() {
    static double timeout = eckit::Resource<double>("$ATLAS_MPI_BARRIER_TIMEOUT", 3.);
    return timeout;
}

enum class BarrierTimeoutAction
{
    ABORT,
    THROW,
    CONTINUE
};

std::string_view to_string_view(BarrierTimeoutAction action) {
    switch (action) {
        case BarrierTimeoutAction::ABORT: return "ABORT";
        case BarrierTimeoutAction::THROW: return "THROW";
        case BarrierTimeoutAction::CONTINUE: return "CONTINUE";
    }
    return "UNKNOWN";
}

BarrierTimeoutAction to_BarrierTimeoutAction(const std::string& str) {
    if (str == "ABORT") return BarrierTimeoutAction::ABORT;
    if (str == "THROW") return BarrierTimeoutAction::THROW;
    if (str == "CONTINUE") return BarrierTimeoutAction::CONTINUE;
    ATLAS_THROW_EXCEPTION("Invalid value for ATLAS_MPI_BARRIER_TIMEOUT_ACTION: " << str << "."
                          << "Valid options are: ABORT, THROW, CONTINUE.");
}

BarrierTimeoutAction ATLAS_MPI_BARRIER_TIMEOUT_ACTION() {
    static BarrierTimeoutAction action = to_BarrierTimeoutAction(eckit::Resource<std::string>("$ATLAS_MPI_BARRIER_TIMEOUT_ACTION", "THROW"));
    return action;
}

bool mpi_barrier_timeout(const mpi::Comm& comm, double seconds) {
    if (seconds <= 0.) {
        return false;
    }
    auto req = comm.iBarrier();
    const auto test_interval = std::chrono::milliseconds(10);
    const auto deadline = std::chrono::steady_clock::now() + std::chrono::duration<double>(seconds);
    while (not req.test()) {
        std::this_thread::sleep_for(test_interval);
        if (std::chrono::steady_clock::now() > deadline) {
            return not req.test();
        }
    }
    return false;
}

std::vector<int> get_mpi_ranks_that_timed_out(const mpi::Comm& comm) {
    int local_timed_out = 0;
    const int size = comm.size();
    const int rank = comm.rank();
    const int tag  = 0;

    std::vector<int> all_timed_out(size);
    all_timed_out[rank] = local_timed_out;

    std::vector<int> peers;
    peers.reserve(size > 0 ? size - 1 : 0);
    std::vector<eckit::mpi::Request> recv_requests;
    std::vector<eckit::mpi::Request> send_requests;
    recv_requests.reserve(size > 0 ? size - 1 : 0);
    send_requests.reserve(size > 0 ? size - 1 : 0);

    for (int peer = 0; peer < size; ++peer) {
        if (peer == rank) {
            continue;
        }
        peers.push_back(peer);
        recv_requests.push_back(comm.iReceive(all_timed_out[peer], peer, tag));
    }

    for (int peer : peers) {
        send_requests.push_back(comm.iSend(local_timed_out, peer, tag));
    }

    std::vector<bool> recv_completed(recv_requests.size(), false);
    std::vector<bool> send_completed(send_requests.size(), false);
    size_t recv_pending = recv_requests.size();
    size_t send_pending = send_requests.size();

    const auto deadline = std::chrono::steady_clock::now() + std::chrono::milliseconds(100);
    while (recv_pending > 0 || send_pending > 0) {
        for (size_t i = 0; i < recv_requests.size(); ++i) {
            if (not recv_completed[i] && recv_requests[i].test()) {
                recv_completed[i] = true;
                --recv_pending;
            }
        }
        for (size_t i = 0; i < send_requests.size(); ++i) {
            if (not send_completed[i] && send_requests[i].test()) {
                send_completed[i] = true;
                --send_pending;
            }
        }

        if (recv_pending == 0 && send_pending == 0) {
            break;
        }

        std::this_thread::sleep_for(std::chrono::milliseconds(100));
        if (std::chrono::steady_clock::now() > deadline) {
            for (size_t i = 0; i < recv_requests.size(); ++i) {
                if (not recv_completed[i]) {
                    all_timed_out[peers[i]] = 1;
                }
            }
            break;
        }
    }

    std::vector<int> timed_out_ranks;
    for (int i = 0; i < all_timed_out.size(); ++i) {
        if (all_timed_out[i]) {
            timed_out_ranks.push_back(i);
        }
    }
    return timed_out_ranks;
}

void debug_mpi_barrier(const eckit::CodeLocation& here, const mpi::Comm& comm,
                                     std::string_view what = {}) {
    if (mpi_barrier_timeout(comm, ATLAS_MPI_BARRIER_TIMEOUT())) {
        std::ostringstream out;
        out << "DEBUG_SYNC[" << comm.name() << "]";
        if (not what.empty()) {
            out << "(" << what << ")";
        }
        out << " @ " << here << "\n";
        out << "ATLAS_DEBUG_SYNC MPI barrier timed out on ranks " << get_mpi_ranks_that_timed_out(comm)
            << " in communicator '" << comm.name() << "' (${ATLAS_MPI_BARRIER_TIMEOUT}=" << ATLAS_MPI_BARRIER_TIMEOUT()
            << ", ${ATLAS_MPI_BARRIER_TIMEOUT_ACTION}=" << to_string_view(ATLAS_MPI_BARRIER_TIMEOUT_ACTION()) << ").\n"
            << "-----------------------------------------\n"
            << "BACKTRACE [rank=" << comm.rank() << "]\n"
            << "-----------------------------------------\n"
            << backtrace() << "\n"
            << "-----------------------------------------";
        switch (ATLAS_MPI_BARRIER_TIMEOUT_ACTION()) {
            case BarrierTimeoutAction::ABORT:
                Log::error() << out.str() << " Calling MPI_Abort..." << std::endl;
                comm.abort();
                break;
            case BarrierTimeoutAction::THROW:
                ATLAS_THROW_EXCEPTION("ATLAS_MPI_BARRIER_TIMEOUT: \n" << out.str());
                break;
            case BarrierTimeoutAction::CONTINUE:
                Log::warning() << out.str() << " Continuing execution despite barrier timeout..." << std::endl;
                break;
        }
    }
}

void flush_and_sync_output_streams() {
    for( auto& channel : {
        &Log::info(),
        &Log::warning(),
        &Log::trace(),
        &Log::debug(),
        &Log::error()} ) {
        if (channel) {
            channel->flush();
        }
    }
    std::cout.flush();
    std::cerr.flush();

    eckit::fdatasync(STDOUT_FILENO);
    eckit::fdatasync(STDERR_FILENO);
}

}  // namespace

void debug_sync(const eckit::CodeLocation& here, const mpi::Comm& comm) {
    flush_and_sync_output_streams();
    debug_mpi_barrier(here, comm);
    Log::info() << "DEBUG_SYNC["<<comm.name()<<"] @ " << here << std::endl;
}

void debug_sync(const eckit::CodeLocation& here) {
    debug_sync(here, mpi::comm());
}

void debug_sync(const eckit::CodeLocation& here, const mpi::Comm& comm, std::string_view what) {
    flush_and_sync_output_streams();
    debug_mpi_barrier(here, comm, what);
    Log::info() << "DEBUG_SYNC["<<comm.name()<<"](" << what << ") @ " << here << std::endl;
}

void debug_sync(const eckit::CodeLocation& here, std::string_view what) {
    debug_sync(here, mpi::comm(), what);
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
