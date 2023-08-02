/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "atlas/library/Library.h"

#include <sstream>
#include <string>

#include "eckit/config/LibEcKit.h"
#include "eckit/config/Resource.h"
#include "eckit/eckit.h"
#include "eckit/filesystem/LocalPathName.h"
#include "eckit/filesystem/PathExpander.h"
#include "eckit/filesystem/PathName.h"
#include "eckit/log/Log.h"
#include "eckit/log/OStreamTarget.h"
#include "eckit/log/PrefixTarget.h"
#include "eckit/runtime/Main.h"
#include "eckit/system/SystemInfo.h"
#include "eckit/types/Types.h"
#include "eckit/utils/Translator.h"
#include "eckit/system/LibraryManager.h"

#if ATLAS_ECKIT_HAVE_ECKIT_585
#include "eckit/linalg/LinearAlgebraDense.h"
namespace {
static bool feature_MKL() {
    return eckit::linalg::LinearAlgebraDense::hasBackend("mkl");
}
}  // namespace
#else
#include "eckit/linalg/LinearAlgebra.h"
namespace {
static bool feature_MKL() {
    return eckit::linalg::LinearAlgebra::hasBackend("mkl");
}
}  // namespace
#endif

#include "atlas_io/Trace.h"

#include "atlas/library/FloatingPointExceptions.h"
#include "atlas/library/Plugin.h"
#include "atlas/library/config.h"
#include "atlas/library/git_sha1.h"
#include "atlas/library/version.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/parallel/omp/omp.h"
#include "atlas/runtime/Exception.h"
#include "atlas/runtime/Log.h"
#include "atlas/runtime/Trace.h"
#include "atlas/util/Config.h"

#if ATLAS_HAVE_TRANS
#if ATLAS_HAVE_ECTRANS
#include "ectrans/transi.h"
#else
#include "transi/trans.h"
#endif
#endif

using eckit::LocalPathName;
using eckit::Main;
using eckit::PathName;

namespace atlas {

//----------------------------------------------------------------------------------------------------------------------

namespace {
std::string str(bool v) {
    return v ? "ON" : "OFF";
}

std::string str(const eckit::system::Library& lib) {
    std::string gitsha1 = lib.gitsha1();
    std::stringstream ss;
    ss << lib.name() << " version (" << lib.version() << "),";
    if (lib.gitsha1() != "not available") {
        ss << "  git-sha1 " << lib.gitsha1(7);
    }
    return ss.str();
}

std::string getEnv(const std::string& env, const std::string& default_value = "") {
    char const* val = ::getenv(env.c_str());
    return val == nullptr ? default_value : std::string(val);
}

bool getEnv(const std::string& env, bool default_value) {
    if (::getenv(env.c_str())) {
        return eckit::Translator<std::string, bool>()(::getenv(env.c_str()));
    }
    return default_value;
}

int getEnv(const std::string& env, int default_value) {
    if (::getenv(env.c_str())) {
        return eckit::Translator<std::string, int>()(::getenv(env.c_str()));
    }
    return default_value;
}

static void add_tokens(std::vector<std::string>& tokens, const std::string& str, const std::string& sep) {
    eckit::Tokenizer tokenize{sep};
    std::vector<std::string> tokenized;
    tokenize(str, tokenized);
    for (auto& t : tokenized) {
        if (not t.empty()) {
            tokens.push_back(eckit::PathExpander::expand(t));
        }
    }
};

static void init_data_paths(std::vector<std::string>& data_paths) {
    ATLAS_ASSERT(eckit::Main::ready());
    add_tokens(data_paths, eckit::LibResource<std::string, Library>("atlas-data-path;$ATLAS_DATA_PATH", ""), ":");
    add_tokens(data_paths, "~atlas/share", ":");
}

}  // namespace

//----------------------------------------------------------------------------------------------------------------------

void initialise(int argc, char** argv) {
    Library::instance().initialise(argc, argv);
}
void initialize(int argc, char** argv) {
    Library::instance().initialise(argc, argv);
}
void initialise() {
    Library::instance().initialise();
}
void initialize() {
    Library::instance().initialise();
}
void finalise() {
    Library::instance().finalise();
}
void finalize() {
    Library::instance().finalise();
}

static Library libatlas;

Library::Library():
    eckit::system::Library(std::string("atlas")),
    debug_(eckit::system::Library::debug()),
    info_(getEnv("ATLAS_INFO", true)),
    warning_(getEnv("ATLAS_WARNING", true)),
    trace_(getEnv("ATLAS_TRACE", false)),
    trace_memory_(getEnv("ATLAS_TRACE_MEMORY", false)),
    trace_barriers_(getEnv("ATLAS_TRACE_BARRIERS", false)),
    trace_report_(getEnv("ATLAS_TRACE_REPORT", false)),
    atlas_io_trace_hook_(::atlas::io::TraceHookRegistry::invalidId()) {
    std::string ATLAS_PLUGIN_PATH = getEnv("ATLAS_PLUGIN_PATH");
#if ATLAS_ECKIT_VERSION_AT_LEAST(1, 25, 0) || ATLAS_ECKIT_DEVELOP
    eckit::system::LibraryManager::addPluginSearchPath(ATLAS_PLUGIN_PATH);
#else
    if (ATLAS_PLUGIN_PATH.size()) {
        std::cout << "WARNING: atlas::Library discovered environment variable ATLAS_PLUGIN_PATH. "
                  << "Currently used version of eckit (" << eckit_version_str() << " [" << eckit_git_sha1() << "]) "
                  << "does not support adding plugin search paths. "
                  << "When using latest eckit develop branch, please rebuild Atlas with "
                  << "CMake argument -DENABLE_ECKIT_DEVELOP=ON\n"
                  << "Alternatively, use combination of environment variables 'PLUGINS_MANIFEST_PATH' "
                  << "and 'LD_LIBRARY_PATH (for UNIX) / DYLD_LIBRARY_PATH (for macOS)' (colon-separated lists)\n" << std::endl;
    }
#endif
}

void Library::registerPlugin(eckit::system::Plugin& plugin) {
    plugins_.push_back(&plugin);
}

void Library::deregisterPlugin(eckit::system::Plugin& plugin) {
    auto it = std::find(plugins_.begin(), plugins_.end(), &plugin);
    ATLAS_ASSERT(it != plugins_.end());
    plugins_.erase(it);
}

std::string Library::cachePath() const {
    auto resource = []() -> std::string {
        return eckit::LibResource<std::string, Library>("atlas-cache-path;$ATLAS_CACHE_PATH", "/tmp/cache");
    };
    static std::string ATLAS_CACHE_PATH = eckit::PathExpander::expand(resource());
    return ATLAS_CACHE_PATH;
}

void Library::registerDataPath(const std::string& path) {
    if (data_paths_.empty()) {
        init_data_paths(data_paths_);
    }
    add_tokens(data_paths_, path, ":");
}


std::string Library::dataPath() const {
    if (data_paths_.empty()) {
        init_data_paths(data_paths_);
    }
    std::vector<std::string> paths = data_paths_;
    auto join                      = [](const std::vector<std::string>& v, const std::string& sep) -> std::string {
        std::stringstream joined;
        for (size_t i = 0; i < v.size(); ++i) {
            if (i > 0) {
                joined << sep;
            }
            joined << v[i];
        }
        return joined.str();
    };
    return join(paths, ":");
}

std::string atlas::Library::linalgSparseBackend() const {
    auto resource = []() -> std::string {
        return eckit::LibResource<std::string, Library>("atlas-linalg-sparse-backend;$ATLAS_LINALG_SPARSE_BACKEND", "");
    };
    static std::string ATLAS_LINALG_SPARSE_BACKEND = resource();
    return ATLAS_LINALG_SPARSE_BACKEND;
}

std::string atlas::Library::linalgDenseBackend() const {
    auto resource = []() -> std::string {
        return eckit::LibResource<std::string, Library>("atlas-linalg-dense-backend;$ATLAS_LINALG_DENSE_BACKEND", "");
    };
    static std::string ATLAS_LINALG_DENSE_BACKEND = resource();
    return ATLAS_LINALG_DENSE_BACKEND;
}


Library& Library::instance() {
    return libatlas;
}


const void* Library::addr() const {
    return this;
}

std::string Library::version() const {
    return atlas::library::version();
}

std::string Library::gitsha1(unsigned int count) const {
    return atlas::library::git_sha1(count);
}

void Library::initialise(int argc, char** argv) {
    if (not Main::ready()) {
        Main::initialise(argc, argv);
        Main::instance().taskID(eckit::mpi::comm("world").rank());
        if (Main::instance().taskID() != 0) {
            eckit::Log::warning().reset();
            eckit::Log::info().reset();
            eckit::Log::debug().reset();
            atlas::Log::debug().reset();
        }
        Log::debug() << "Atlas initialised eckit::Main.\n";
        if (eckit::mpi::comm("world").size() > 1) {
            Log::debug() << "--> Only MPI rank 0 is logging. Please initialise eckit::Main \n"
                            "    before to avoid this behaviour.\n";
        }
    }
    initialise();
}


void Library::initialise(const eckit::Parametrisation& config) {
    ATLAS_ASSERT(eckit::Main::ready());

    if (initialized_) {
        return;
    }
    initialized_ = true;
    if (config.has("log")) {
        config.get("log.info", info_);
        config.get("log.trace", trace_);
        config.get("log.warning", warning_);
        config.get("log.debug", debug_);
    }
    if (config.has("trace")) {
        config.get("trace.barriers", trace_barriers_);
        config.get("trace.report", trace_report_);
        config.get("trace.memory", trace_memory_);
    }

    if (not debug_) {
        debug_channel_.reset();
    }
    if (not trace_) {
        trace_channel_.reset();
    }
    if (not info_) {
        info_channel_.reset();
    }
    if (not warning_) {
        warning_channel_.reset();
    }

    auto& out = [&]() -> eckit::Channel& {
        if (getEnv("ATLAS_LOG_RANK", 0) == int(mpi::rank())) {
            return Log::debug();
        }
        static eckit::Channel sink;
        return sink;
    }();

    library::enable_floating_point_exceptions();
    library::enable_atlas_signal_handler();

    if (data_paths_.empty()) {
        init_data_paths(data_paths_);
    }

    atlas_io_trace_hook_ =
        atlas::io::TraceHookRegistry::add([](const eckit::CodeLocation& loc, const std::string& title) {
            struct Adaptor : public atlas::io::TraceHook {
                Adaptor(const eckit::CodeLocation& loc, const std::string& title): trace{loc, title} {}
                atlas::Trace trace;
            };
            return std::make_unique<Adaptor>(loc, title);
        });


    // Summary
    if (getEnv("ATLAS_LOG_RANK", 0) == int(mpi::rank())) {
        out << "Executable        [" << Main::instance().name() << "]\n";
        out << " \n";
        out << "  current dir     [" << PathName(LocalPathName::cwd()).fullName() << "]\n";
        out << " \n";
        out << "  MPI\n";
        out << "    communicator  [" << mpi::comm() << "] \n";
        out << "    size          [" << mpi::size() << "] \n";
        out << "    rank          [" << mpi::rank() << "] \n";
        out << "  OMP\n";
        out << "    max_threads   [" << atlas_omp_get_max_threads() << "] \n";
        out << " \n";
        out << "  log.info        [" << str(info_) << "] \n";
        out << "  log.trace       [" << str(trace()) << "] \n";
        out << "  log.debug       [" << str(debug()) << "] \n";
        out << "  trace.barriers  [" << str(traceBarriers()) << "] \n";
        out << "  trace.report    [" << str(trace_report_) << "] \n";
        out << "  trace.memory    [" << str(trace_memory_) << "] \n";
        out << " \n";
        out << atlas::Library::instance().information();
        out << std::flush;
    }
}


void Library::initialise() {
    initialise(util::NoConfig());
}

void Library::finalise() {
    if( atlas_io_trace_hook_ != atlas::io::TraceHookRegistry::invalidId() ) {
        atlas::io::TraceHookRegistry::disable(atlas_io_trace_hook_);
        atlas_io_trace_hook_ = atlas::io::TraceHookRegistry::invalidId();
    }

    if (ATLAS_HAVE_TRACE && trace_report_) {
        Log::info() << atlas::Trace::report() << std::endl;
    }

    if (getEnv("ATLAS_FINALISES_MPI", false)) {
        Log::debug() << "ATLAS_FINALISES_MPI is set: calling atlas::mpi::finalize()" << std::endl;
        mpi::finalise();
    }

    // Make sure that these specialised channels that wrap Log::info() are
    // destroyed before eckit::Log::info gets destroyed.
    // Just in case someone still tries to log, we reset to empty channels.
    trace_channel_.reset(new eckit::Channel());

    Log::debug() << "Atlas finalised" << std::endl;

    Log::flush();

    if (debugChannel()) {
        debug_channel_.reset(new eckit::Channel(new eckit::PrefixTarget("ATLAS_DEBUG")));
    }
    if (infoChannel()) {
        info_ = false;
        info_channel_.reset(new eckit::Channel(new eckit::PrefixTarget("ATLAS_INFO")));
    }
    if (warningChannel()) {
        warning_ = false;
        warning_channel_.reset(new eckit::Channel(new eckit::PrefixTarget("ATLAS_WARNING")));
    }
    initialized_ = false;
}

eckit::Channel& Library::infoChannel() const {
    if (info_) {
        return eckit::Log::info();
    }
    else if (!info_channel_) {
        info_channel_.reset(new eckit::Channel());
    }
    return *info_channel_;
}


eckit::Channel& Library::warningChannel() const {
    if (warning_) {
        return eckit::Log::warning();
    }
    else if (!warning_channel_) {
        warning_channel_.reset(new eckit::Channel());
    }
    return *warning_channel_;
}


eckit::Channel& Library::traceChannel() const {
    if (trace_channel_) {
        return *trace_channel_;
    }
    if (trace_) {
        trace_channel_.reset(
            new eckit::Channel(new eckit::PrefixTarget("ATLAS_TRACE", new eckit::OStreamTarget(eckit::Log::info()))));
    }
    else {
        trace_channel_.reset(new eckit::Channel());
    }
    return *trace_channel_;
}

eckit::Channel& Library::debugChannel() const {
    if (debug_channel_) {
        return *debug_channel_;
    }
    if (debug_) {
        debug_channel_.reset(new eckit::Channel(new eckit::PrefixTarget("ATLAS_DEBUG")));
    }
    else {
        debug_channel_.reset(new eckit::Channel());
    }
    return *debug_channel_;
}

//----------------------------------------------------------------------------------------------------------------------

void Library::Information::print(std::ostream& out) const {
    out << "atlas version (" << atlas::Library::instance().version() << "), "
        << "git-sha1 " << atlas::Library::instance().gitsha1(7) << '\n';
    out << " \n";
    out << "  Build:" << std::endl;
    out << "    build type      : " << ATLAS_BUILD_TYPE << '\n'
        << "    timestamp       : " << ATLAS_BUILD_TIMESTAMP << '\n'
        << "    source dir      : " << ATLAS_DEVELOPER_SRC_DIR << '\n'
        << "    build dir       : " << ATLAS_DEVELOPER_BIN_DIR << '\n'
        << "    op. system      : " << ATLAS_OS_NAME << " (" << ATLAS_OS_STR << ")" << '\n'
        << "    processor       : " << ATLAS_SYS_PROCESSOR << std::endl
        << "    c compiler      : " << ATLAS_C_COMPILER_ID << " " << ATLAS_C_COMPILER_VERSION << '\n'
        << "      flags         : " << ATLAS_C_FLAGS << '\n'
        << "    c++ compiler    : " << ATLAS_CXX_COMPILER_ID << " " << ATLAS_CXX_COMPILER_VERSION << '\n'
        << "      flags         : " << ATLAS_CXX_FLAGS << '\n'
#ifdef ATLAS_Fortran_COMPILER
        << "    fortran compiler: " << ATLAS_Fortran_COMPILER_ID << " " << ATLAS_Fortran_COMPILER_VERSION << '\n'
        << "      flags         : " << ATLAS_Fortran_FLAGS << '\n'
#else
        << "    fortran         : NO " << '\n'
#endif
        << " \n";

    bool feature_fortran(ATLAS_HAVE_FORTRAN);
    bool feature_OpenMP(ATLAS_HAVE_OMP);
    bool feature_ecTrans(ATLAS_HAVE_ECTRANS);
    bool feature_FFTW(ATLAS_HAVE_FFTW);
    bool feature_Eigen(ATLAS_HAVE_EIGEN);
    bool feature_Tesselation(ATLAS_HAVE_TESSELATION);
    bool feature_PROJ(ATLAS_HAVE_PROJ);
    bool feature_BoundsChecking(ATLAS_ARRAYVIEW_BOUNDS_CHECKING);
    bool feature_Init_sNaN(ATLAS_INIT_SNAN);
    bool feature_MPI(false);
#if ATLAS_HAVE_MPI
    feature_MPI = true;
#endif
    std::string array_data_store = "Native";
#if ATLAS_HAVE_GRIDTOOLS_STORAGE
    array_data_store = "Gridtools-host";
#if ATLAS_GRIDTOOLS_STORAGE_BACKEND_CUDA
    array_data_store = "GridTools-CUDA";
#endif
#endif
    out << "  Features:" << '\n'
        << "    Fortran        : " << str(feature_fortran) << '\n'
        << "    MPI            : " << str(feature_MPI) << '\n'
        << "    OpenMP         : " << str(feature_OpenMP) << '\n'
        << "    BoundsChecking : " << str(feature_BoundsChecking) << '\n'
        << "    Init_sNaN      : " << str(feature_Init_sNaN) << '\n'
        << "    ecTrans        : " << str(feature_ecTrans) << '\n'
        << "    FFTW           : " << str(feature_FFTW) << '\n'
        << "    Eigen          : " << str(feature_Eigen) << '\n'
        << "    MKL            : " << str(feature_MKL()) << '\n'
        << "    Tesselation    : " << str(feature_Tesselation) << '\n'
        << "    PROJ           : " << str(feature_PROJ) << '\n'
        << "    ArrayDataStore : " << array_data_store << '\n'
        << "    idx_t          : " << ATLAS_BITS_LOCAL << " bit integer" << '\n'
        << "    gidx_t         : " << ATLAS_BITS_GLOBAL << " bit integer" << '\n';

    auto& plugins = Library::instance().plugins();
    if (!plugins.empty()) {
        out << "    \n  Plugins: \n";
        for (auto& plugin : plugins) {
            out << "    " << str(*plugin) << '\n';
        }
    }

    out << "    \n  Dependencies: \n";
    out << "    ecbuild version (" << ECBUILD_VERSION << ")" << '\n';
    if (Library::exists("eckit")) {
        out << "    " << str(Library::lookup("eckit")) << '\n';
    }
    if (Library::exists("fckit")) {
        out << "    " << str(Library::lookup("fckit")) << '\n';
    }

#if ATLAS_HAVE_TRANS
#if ATLAS_HAVE_ECTRANS
    out << "    ectrans version (" << ectrans_version() << "), "
        << "git-sha1 " << ectrans_git_sha1_abbrev(7) << '\n';
    out << "    fiat version (" << ectrans_fiat_version() << "), "
        << "git-sha1 " << ectrans_fiat_git_sha1_abbrev(7) << '\n';
#else
    out << "    transi version (" << transi_version() << "), "
        << "git-sha1 " << transi_git_sha1_abbrev(7) << '\n';
    out << "    trans version (" << trans_version() << "), "
        << "git-sha1 " << trans_git_sha1_abbrev(7) << '\n';
#endif
#endif
}

}  // namespace atlas
