/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "atlas/library/FloatingPointExceptions.h"

#include <cfenv>
#include <cstring>
#include <iomanip>

#include "eckit/config/LibEcKit.h"
#include "eckit/config/Resource.h"
#include "eckit/runtime/Main.h"
#include "eckit/utils/StringTools.h"
#include "eckit/utils/Translator.h"

#include "atlas/library/config.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/runtime/Exception.h"
#include "atlas/runtime/Log.h"
#include "atlas/runtime/Trace.h"

namespace {
int getEnv(const std::string& env, int default_value) {
    const char* cenv = ::getenv(env.c_str());
    if (cenv != nullptr) {
        return eckit::Translator<std::string, int>()(cenv);
    }
    return default_value;
}
}  // namespace

static int atlas_feenableexcept(int excepts) {
#if ATLAS_HAVE_FEENABLEEXCEPT
    return ::feenableexcept(excepts);
#else
    return 0;
#endif
}

static int atlas_fedisableexcept(int excepts) {
#if ATLAS_HAVE_FEDISABLEEXCEPT
    return ::fedisableexcept(excepts);
#else
    return 0;
#endif
}

namespace atlas {
namespace library {

static std::map<std::string, int> str_to_except = {{"FE_INVALID", FE_INVALID},     {"FE_INEXACT", FE_INEXACT},
                                                   {"FE_DIVBYZERO", FE_DIVBYZERO}, {"FE_OVERFLOW", FE_OVERFLOW},
                                                   {"FE_UNDERFLOW", FE_UNDERFLOW}, {"FE_ALL_EXCEPT", FE_ALL_EXCEPT}};
static std::map<int, std::string> except_to_str = {{FE_INVALID, "FE_INVALID"},     {FE_INEXACT, "FE_INEXACT"},
                                                   {FE_DIVBYZERO, "FE_DIVBYZERO"}, {FE_OVERFLOW, "FE_OVERFLOW"},
                                                   {FE_UNDERFLOW, "FE_UNDERFLOW"}, {FE_ALL_EXCEPT, "FE_ALL_EXCEPT"}};

static std::map<std::string, int> str_to_signal = {{"SIGINT", SIGINT},  {"SIGILL", SIGILL},   {"SIGABRT", SIGABRT},
                                                   {"SIGFPE", SIGFPE},  {"SIGKILL", SIGKILL}, {"SIGSEGV", SIGSEGV},
                                                   {"SIGTERM", SIGTERM}};
static std::map<int, std::string> signal_to_str = {{SIGINT, "SIGINT"},  {SIGILL, "SIGILL"},   {SIGABRT, "SIGABRT"},
                                                   {SIGFPE, "SIGFPE"},  {SIGKILL, "SIGKILL"}, {SIGSEGV, "SIGSEGV"},
                                                   {SIGTERM, "SIGTERM"}};

// ------------------------------------------------------------------------------------

class Signal {
    // Not sure if this should be made public (in header file) just yet

public:
    Signal();

    Signal(int signum);

    Signal(int signum, signal_action_t);

    Signal(int signum, signal_handler_t signal_handler);

    operator int() const { return signum_; }
    int signum() const { return signum_; }

    const std::string& code() const { return signal_to_str[signum_]; }

    std::string str() const { return str_; }
    const signal_handler_t& handler() const { return signal_action_.sa_handler; }
    const struct sigaction* action() const { return &signal_action_; }

private:
    friend std::ostream& operator<<(std::ostream&, const Signal&);

    int signum_;
    std::string str_;
    struct sigaction signal_action_;
};

// ------------------------------------------------------------------------------------

class Signals {
    // Not sure if this should be made public (in header file) just yet

private:
    Signals();

public:
    static Signals& instance();
    void setSignalHandlers();
    void setSignalHandler(const Signal&);
    void restoreSignalHandler(int signum);
    void restoreAllSignalHandlers();
    const Signal& signal(int signum) const;

private:
    using registered_signals_t = std::map<int, Signal>;
    registered_signals_t registered_signals_;
    eckit::Channel& out_;
};

// ------------------------------------------------------------------------------------


// ------------------------------------------------------------------------------------

[[noreturn]] void atlas_signal_handler(int signum, siginfo_t* si, void* /*unused*/) {
    Signal signal = Signals::instance().signal(signum);

    std::string signal_code;
    if (signum == SIGFPE) {
        switch (si->si_code) {
            case FPE_FLTDIV:
                signal_code = " [FE_DIVBYZERO]";
                break;
            case FPE_FLTINV:
                signal_code = " [FE_INVALID]";
                break;
            case FPE_FLTOVF:
                signal_code = " [FE_OVERFLOW]";
                break;
            case FPE_FLTUND:
                signal_code = " [FE_UNDERFLOW]";
                break;
            case FPE_FLTRES:
                signal_code = " [FE_INEXACT]";
                break;
        }
    }

    std::ostream& out = Log::error();
    out << "\n"
        << "=========================================\n"
        << signal << signal_code << " (signal intercepted by atlas)\n";
    out << "-----------------------------------------\n"
        << "BACKTRACE\n"
        << "-----------------------------------------\n"
        << backtrace() << "\n"
        << "=========================================\n"
        << std::endl;

    Signals::instance().restoreSignalHandler(signum);
    eckit::LibEcKit::instance().abort();

    // Just in case we end up here, which normally we shouldn't.
    std::_Exit(EXIT_FAILURE);
}

//----------------------------------------------------------------------------------------------------------------------


Signals::Signals():
    out_([&]() -> eckit::Channel& {
        if (getEnv("ATLAS_LOG_RANK", 0) == int(mpi::rank())) {
            return Log::debug();
        }
        static eckit::Channel sink;
        return sink;
    }()) {}

Signals& Signals::instance() {
    static Signals signals;
    return signals;
}

void Signals::restoreSignalHandler(int signum) {
    if (registered_signals_.find(signum) != registered_signals_.end()) {
        out_ << "\n";
        std::signal(signum, SIG_DFL);
        out_ << "Atlas restored default signal handler for signal " << std::setw(7) << std::left
             << registered_signals_[signum].code() << " [" << registered_signals_[signum] << "]\n";
        out_ << std::endl;
        registered_signals_.erase(signum);
    }
}

void Signals::restoreAllSignalHandlers() {
    out_ << "\n";
    for (registered_signals_t::const_iterator it = registered_signals_.begin(); it != registered_signals_.end(); ++it) {
        std::signal(it->first, SIG_DFL);
        out_ << "Atlas restored default signal handler for signal " << std::setw(7) << std::left << it->second.code()
             << " [" << it->second.str() << "]\n";
    }
    out_ << std::endl;
    registered_signals_.clear();
}

const Signal& Signals::signal(int signum) const {
    return registered_signals_.at(signum);
}

std::ostream& operator<<(std::ostream& out, const Signal& signal) {
    out << signal.str();
    return out;
}

void Signals::setSignalHandlers() {
    setSignalHandler(SIGINT);
    setSignalHandler(SIGILL);
    setSignalHandler(SIGABRT);
    setSignalHandler(SIGFPE);
    setSignalHandler(SIGSEGV);
    setSignalHandler(SIGTERM);
}

void Signals::setSignalHandler(const Signal& signal) {
    registered_signals_[signal] = signal;
    sigaction(signal, signal.action(), nullptr);
    out_ << "Atlas registered signal handler for signal " << std::setw(7) << std::left << signal.code() << " ["
         << signal << "]" << std::endl;
}


Signal::Signal(): signum_(0), str_() {
    signal_action_.sa_handler = SIG_DFL;
}

Signal::Signal(int signum): Signal(signum, atlas_signal_handler) {}

Signal::Signal(int signum, signal_handler_t signal_handler): signum_(signum), str_(strsignal(signum)) {
    memset(&signal_action_, 0, sizeof(signal_action_));
    sigemptyset(&signal_action_.sa_mask);
    signal_action_.sa_handler = signal_handler;
    signal_action_.sa_flags   = 0;
}

Signal::Signal(int signum, signal_action_t signal_action): signum_(signum), str_(strsignal(signum)) {
    memset(&signal_action_, 0, sizeof(signal_action_));
    sigemptyset(&signal_action_.sa_mask);
    signal_action_.sa_sigaction = signal_action;
    signal_action_.sa_flags     = SA_SIGINFO;
}


void enable_floating_point_exceptions() {
    auto& out = [&]() -> eckit::Channel& {
        if (getEnv("ATLAS_LOG_RANK", 0) == int(mpi::rank())) {
            return Log::debug();
        }
        static eckit::Channel sink;
        return sink;
    }();

    // Following line gives runtime errors with Cray 8.6 due to compiler bug (but works with Cray 8.5 and Cray 8.7)
    //   std::vector<std::string> floating_point_exceptions = eckit::Resource<std::vector<std::string>>( "atlasFPE;$ATLAS_FPE", {"false"} );
    // Instead, manually access environment
    std::vector<std::string> floating_point_exceptions{"false"};
    const char* ATLAS_FPE = ::getenv("ATLAS_FPE");
    if (ATLAS_FPE != nullptr) {
        std::vector<std::string> tmp = eckit::Translator<std::string, std::vector<std::string>>()(ATLAS_FPE);
        floating_point_exceptions    = tmp;
        // Above trick with "tmp" is what avoids the Cray 8.6 compiler bug
    }
    else {
        if (eckit::Main::ready()) {
            floating_point_exceptions = eckit::Resource<std::vector<std::string>>("atlasFPE", {"false"});
        }
    }
    {
        bool _enable = false;
        int _excepts = 0;
        auto enable  = [&](int except) {
            _excepts |= except;
            _enable = true;
            out << "Atlas enabled floating point exception " << except_to_str[except] << std::endl;
        };
        bool skip_map = false;
        if (floating_point_exceptions.size() == 1) {
            std::string s = eckit::StringTools::lower(floating_point_exceptions[0]);
            if (s == "no" || s == "off" || s == "false" || s == "0") {
                _enable  = false;
                skip_map = true;
            }
            else if (s == "yes" || s == "on" || s == "true" || s == "1") {
                enable(FE_INVALID);
                enable(FE_DIVBYZERO);
                enable(FE_OVERFLOW);
                skip_map = true;
            }
        }
        if (not skip_map) {
            for (auto& s : floating_point_exceptions) {
                if (str_to_except.find(s) == str_to_except.end()) {
                    throw eckit::UserError(
                        s + " is not a valid floating point exception code. "
                            "Valid codes: [FE_INVALID,FE_INEXACT,FE_DIVBYZERO,FE_OVERFLOW,FE_ALL_EXCEPT]",
                        Here());
                }
                enable(str_to_except[s]);
            }
        }
        if (_enable) {
            atlas_feenableexcept(_excepts);
        }
    }
}

bool enable_floating_point_exception(int except) {
    auto check_flag = [](int flags, int bits) -> bool { return (flags & bits) == bits; };
    int previous = atlas_feenableexcept(except);
    return !check_flag(previous,except);
}

bool disable_floating_point_exception(int except) {
    auto check_flag = [](int flags, int bits) -> bool { return (flags & bits) == bits; };
    int previous = atlas_fedisableexcept(except);
    return !check_flag(previous,except);
}

bool enable_floating_point_exception(const std::string& floating_point_exception) {
    auto it = str_to_except.find(floating_point_exception);
    if (it == str_to_except.end()) {
        throw eckit::UserError(
            floating_point_exception + " is not a valid floating point exception code. "
                "Valid codes: [FE_INVALID,FE_INEXACT,FE_DIVBYZERO,FE_OVERFLOW,FE_ALL_EXCEPT]",
            Here());
    }
    return enable_floating_point_exception(it->second);
}

bool disable_floating_point_exception(const std::string& floating_point_exception) {
    auto it = str_to_except.find(floating_point_exception);
    if (it == str_to_except.end()) {
        throw eckit::UserError(
            floating_point_exception + " is not a valid floating point exception code. "
                "Valid codes: [FE_INVALID,FE_INEXACT,FE_DIVBYZERO,FE_OVERFLOW,FE_ALL_EXCEPT]",
            Here());
    }
    return disable_floating_point_exception(it->second);
}


void enable_atlas_signal_handler() {
    bool enable = false;
    const char* ATLAS_SIGNAL_HANDLER = ::getenv("ATLAS_SIGNAL_HANDLER");
    if (ATLAS_SIGNAL_HANDLER != nullptr) {
        bool tmp = eckit::Translator<std::string, bool>()(ATLAS_SIGNAL_HANDLER);
        enable   = tmp;
        // Above trick with "tmp" is what avoids the Cray 8.6 compiler bug
    }
    else {
        if (eckit::Main::ready()) {
            enable = eckit::Resource<bool>("atlasSignalHandler", false);
        }
    }

    if (enable) {
        Signals::instance().setSignalHandlers();
    }
}

}  // namespace library
}  // namespace atlas
