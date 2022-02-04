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

#include <algorithm>
#include <chrono>
#include <exception>
#include <iomanip>
#include <sstream>
#include <string>
#include <thread>

#include "eckit/config/LibEcKit.h"
#include "eckit/config/Resource.h"
#include "eckit/eckit.h"
#include "eckit/log/FileTarget.h"
#include "eckit/log/PrefixTarget.h"
#include "eckit/mpi/Comm.h"
#include "eckit/runtime/Main.h"
#include "eckit/testing/Test.h"
#include "eckit/types/Types.h"

#include "atlas/library/Library.h"
#include "atlas/runtime/Exception.h"
#include "atlas/runtime/Log.h"
#include "atlas/runtime/Trace.h"
#include "atlas/runtime/trace/StopWatch.h"
#include "atlas/util/Config.h"
#include "atlas/util/Point.h"

namespace atlas {
namespace test {

using eckit::types::is_approximately_equal;

class Test;
static Test* current_test_{nullptr};

static size_t ATLAS_MAX_FAILED_EXPECTS() {
    static size_t v = size_t(eckit::Resource<long>("$ATLAS_MAX_FAILED_EXPECTS", 100));
    return v;
}

class Test {
    struct Failure {
        std::string message;
        eckit::CodeLocation location;
    };

public:
    Test(const std::string& description, const eckit::CodeLocation& location):
        description_(description), location_(location) {
        current_test_ = this;
    }
    ~Test() { current_test_ = nullptr; }
    void expect_failed(const std::string& message, const eckit::CodeLocation& location) {
        failures_.emplace_back(Failure{message, location});
        eckit::Log::error() << message << std::endl;
        if (failures_.size() == ATLAS_MAX_FAILED_EXPECTS()) {
            std::stringstream msg;
            msg << "Maximum number of allowed EXPECTS have failed (${ATLAS_MAX_FAILED_EXPECTS}="
                << ATLAS_MAX_FAILED_EXPECTS() << ").";
            throw eckit::testing::TestException(msg.str(), location_);
        }
    }
    bool failed() const { return failures_.size() > 0; }
    void throw_on_failed_expects() {
        if (failed()) {
            std::stringstream msg;
            msg << failures_.size() << " EXPECTS have failed";
            throw eckit::testing::TestException(msg.str(), location_);
        }
    }
    const std::string& description() const { return description_; }

private:
    std::vector<Failure> failures_;
    std::string description_;
    eckit::CodeLocation location_;
};

Test& current_test() {
    ATLAS_ASSERT(current_test_);
    return *current_test_;
}

//----------------------------------------------------------------------------------------------------------------------

#ifdef MAYBE_UNUSED
#elif defined(__GNUC__)
#define MAYBE_UNUSED __attribute__((unused))
#else
#define MAYBE_UNUSED
#endif

// Redefine macro's defined in "eckit/testing/Test.h" to include trace
// information
#undef CASE
#define CASE(description)                                                                                              \
    void UNIQUE_NAME2(test_, __LINE__)(std::string&, int&, int);                                                       \
    static const eckit::testing::TestRegister UNIQUE_NAME2(test_registration_, __LINE__)(                              \
        description, &UNIQUE_NAME2(test_, __LINE__));                                                                  \
    void UNIQUE_NAME2(traced_test_, __LINE__)(std::string & _test_subsection, int& _num_subsections, int _subsection); \
    void UNIQUE_NAME2(test_, __LINE__)(std::string & _test_subsection, int& _num_subsections, int _subsection) {       \
        Test UNIQUE_NAME2(testobj_, __LINE__)(description, Here());                                                    \
        ATLAS_TRACE(description);                                                                                      \
        UNIQUE_NAME2(traced_test_, __LINE__)                                                                           \
        (_test_subsection, _num_subsections, _subsection);                                                             \
        current_test().throw_on_failed_expects();                                                                      \
        if (atlas::test::barrier_timeout(atlas::test::ATLAS_MPI_BARRIER_TIMEOUT())) {                                  \
            atlas::Log::warning() << "\nWARNING: Test \"" << description                                               \
                                  << "\" failed with MPI deadlock.  (${ATLAS_MPI_BARRIER_TIMEOUT}="                    \
                                  << atlas::test::ATLAS_MPI_BARRIER_TIMEOUT() << ").\nCalling MPI_Abort..."            \
                                  << std::endl;                                                                        \
            eckit::mpi::comm().abort();                                                                                \
        }                                                                                                              \
    }                                                                                                                  \
    void UNIQUE_NAME2(traced_test_, __LINE__)(MAYBE_UNUSED std::string & _test_subsection,                             \
                                              MAYBE_UNUSED int& _num_subsections, MAYBE_UNUSED int _subsection)

#undef SECTION
#define SECTION(name)                                                                            \
    _num_subsections += 1;                                                                       \
    _test_subsection = (name);                                                                   \
    if ((_num_subsections - 1) == _subsection) {                                                 \
        atlas::Log::info() << "Running section \"" << _test_subsection << "\" ..." << std::endl; \
    }                                                                                            \
    if ((_num_subsections - 1) == _subsection)                                                   \
    ATLAS_TRACE_SCOPE(_test_subsection)

#ifdef EXPECT_EQ
#undef EXPECT_EQ
#endif
#ifdef EXPECT_APPROX_EQ
#undef EXPECT_APPROX_EQ
#endif

#ifdef EXPECT
#undef EXPECT
#endif

#define REQUIRE(expr)                                                                       \
    do {                                                                                    \
        if (!(expr)) {                                                                      \
            throw eckit::testing::TestException("EXPECT condition failed: " #expr, Here()); \
        }                                                                                   \
    } while (false)

#define EXPECT(expr)                                                                 \
    do {                                                                             \
        if (!(expr)) {                                                               \
            current_test().expect_failed("EXPECT condition failed: " #expr, Here()); \
        }                                                                            \
    } while (false)

template <typename Value>
struct Printer {
    static void print(std::ostream& out, const Value& v) { out << v; }
};

template <>
struct Printer<double> {
    static void print(std::ostream& out, const double& v) { out << std::fixed << std::setprecision(12) << v; }
};

template <>
struct Printer<PointLonLat> {
    static void print(std::ostream& out, const PointLonLat& v) { out << std::fixed << std::setprecision(12) << v; }
};

template <>
struct Printer<eckit::CodeLocation> {
    static void print(std::ostream& out, const eckit::CodeLocation& location) {
        out << eckit::PathName{location.file()}.baseName() << " +" << location.line();
    }
};

template <typename Value>
struct PrintValue {
    const Value& value;
    PrintValue(const Value& v): value(v) {}
    void print(std::ostream& out) const { Printer<Value>::print(out, value); }
    friend std::ostream& operator<<(std::ostream& out, const PrintValue& v) {
        v.print(out);
        return out;
    }
};

template <typename Value>
PrintValue<Value> print(const Value& v) {
    return PrintValue<Value>(v);
}

bool approx_eq(const float& v1, const float& v2) {
    return is_approximately_equal(v1, v2);
}
bool approx_eq(const float& v1, const float& v2, const float& t) {
    return is_approximately_equal(v1, v2, t);
}
bool approx_eq(const double& v1, const double& v2) {
    return is_approximately_equal(v1, v2);
}
bool approx_eq(const double& v1, const double& v2, const double& t) {
    return is_approximately_equal(v1, v2, t);
}
bool approx_eq(const Point2& v1, const Point2& v2) {
    return approx_eq(v1[0], v2[0]) && approx_eq(v1[1], v2[1]);
}
bool approx_eq(const Point2& v1, const Point2& v2, const double& t) {
    return approx_eq(v1[0], v2[0], t) && approx_eq(v1[1], v2[1], t);
}

template <typename T1, typename T2>
std::string expect_message(const std::string& condition, const T1& lhs, const T2& rhs, const eckit::CodeLocation& loc) {
    std::stringstream msg;
    msg << eckit::Colour::red << condition << " FAILED @ " << print(loc) << eckit::Colour::reset << "\n"
        << eckit::Colour::red << " --> lhs = " << print(lhs) << eckit::Colour::reset << "\n"
        << eckit::Colour::red << " --> rhs = " << print(rhs) << eckit::Colour::reset;
    return msg.str();
}

#define EXPECT_EQ(lhs, rhs)                                                                                            \
    do {                                                                                                               \
        if (!(lhs == rhs)) {                                                                                           \
            current_test().expect_failed(expect_message("EXPECT_EQ( " #lhs ", " #rhs " )", lhs, rhs, Here()), Here()); \
        }                                                                                                              \
    } while (false)

#define __EXPECT_APPROX_EQ(lhs, rhs)                                                                                 \
    do {                                                                                                             \
        if (!(approx_eq(lhs, rhs))) {                                                                                \
            current_test().expect_failed(expect_message("EXPECT_APPROX_EQ( " #lhs ", " #rhs " )", lhs, rhs, Here()), \
                                         Here());                                                                    \
        }                                                                                                            \
    } while (false)

#define __EXPECT_APPROX_EQ_TOL(lhs, rhs, tol)                                                                  \
    do {                                                                                                       \
        if (!(approx_eq(lhs, rhs, tol))) {                                                                     \
            current_test().expect_failed(                                                                      \
                expect_message("EXPECT_APPROX_EQ( " #lhs ", " #rhs ", " #tol " )", lhs, rhs, Here()), Here()); \
        }                                                                                                      \
    } while (false)

#define EXPECT_APPROX_EQ(...) __ATLAS_SPLICE(__EXPECT_APPROX_EQ__, __ATLAS_NARG(__VA_ARGS__))(__VA_ARGS__)
#define __EXPECT_APPROX_EQ__2 __EXPECT_APPROX_EQ
#define __EXPECT_APPROX_EQ__3 __EXPECT_APPROX_EQ_TOL


//----------------------------------------------------------------------------------------------------------------------


static double ATLAS_MPI_BARRIER_TIMEOUT() {
    static double v = eckit::Resource<double>("$ATLAS_MPI_BARRIER_TIMEOUT", 3.);
    return v;
}

static int barrier_timeout(double seconds) {
    auto req = eckit::mpi::comm().iBarrier();
    runtime::trace::StopWatch watch;
    while (not req.test()) {
        watch.start();
        std::this_thread::sleep_for(std::chrono::milliseconds(100));
        watch.stop();
        if (watch.elapsed() > seconds) {
            return 1;
        }
    }
    return 0;
}


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
    if (eckit::Log::debug())
        eckit::Log::debug().addTarget(target);
}

void debug_setTarget(eckit::LogTarget* target) {
    for (std::string libname : eckit::system::Library::list()) {
        const eckit::system::Library& lib = eckit::system::Library::lookup(libname);
        if (lib.debug()) {
            lib.debugChannel().setTarget(new eckit::PrefixTarget(debug_prefix(libname), target));
        }
    }
    if (eckit::Log::debug())
        eckit::Log::debug().setTarget(target);
}

void debug_reset() {
    for (std::string libname : eckit::system::Library::list()) {
        const eckit::system::Library& lib = eckit::system::Library::lookup(libname);
        if (lib.debug()) {
            lib.debugChannel().reset();
        }
    }
    if (eckit::Log::debug())
        eckit::Log::debug().reset();
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

void setEnv(const std::string& env, bool value) {
    constexpr int DO_NOT_REPLACE_IF_EXISTS = 0;
    ::setenv(env.c_str(), eckit::Translator<bool, std::string>()(value).c_str(), DO_NOT_REPLACE_IF_EXISTS);
}

}  // namespace

struct AtlasTestEnvironment {
    using Config = util::Config;

    AtlasTestEnvironment(int argc, char* argv[]) {
        eckit::Main::initialise(argc, argv);
        eckit::Main::instance().taskID(eckit::mpi::comm("world").rank());

        long log_rank    = getEnv("ATLAS_LOG_RANK", 0);
        bool use_logfile = getEnv("ATLAS_LOG_FILE", false);

        auto rank_str = []() {
            int d           = digits(eckit::mpi::comm().size());
            std::string str = std::to_string(eckit::Main::instance().taskID());
            for (int i = str.size(); i < d; ++i)
                str = "0" + str;
            return str;
        };

        if (use_logfile) {
            eckit::LogTarget* logfile =
                new eckit::FileTarget(eckit::Main::instance().displayName() + ".log.p" + rank_str());

            if (eckit::Main::instance().taskID() == log_rank) {
                eckit::Log::info().addTarget(logfile);
                eckit::Log::warning().addTarget(logfile);
                eckit::Log::error().setTarget(logfile);
                debug_addTarget(logfile);
            }
            else {
                eckit::Log::info().setTarget(logfile);
                eckit::Log::warning().setTarget(logfile);
                eckit::Log::error().setTarget(logfile);
                debug_setTarget(logfile);
            }
        }
        else {
            if (eckit::Main::instance().taskID() != log_rank) {
                eckit::Log::info().reset();
                eckit::Log::warning().reset();
                debug_reset();
            }
            eckit::Log::error().reset();
        }
        Log::error().addTarget(new eckit::PrefixTarget("[" + std::to_string(eckit::mpi::comm().rank()) + "]"));

        eckit::LibEcKit::instance().setAbortHandler([] {
            Log::error() << "Calling MPI_Abort" << std::endl;
            eckit::mpi::comm().abort();
        });

        setEnv("ATLAS_FPE", true);
        setEnv("ATLAS_SIGNAL_HANDLER", true);

        atlas::initialize();
        eckit::mpi::comm().barrier();
    }

    ~AtlasTestEnvironment() {
        atlas::finalize();
        atlas::mpi::finalize();
    }
};


//----------------------------------------------------------------------------------------------------------------------

template <typename Environment>
int run(int argc, char* argv[]) {
    Environment env(argc, argv);
    int errors = eckit::testing::run_tests(argc, argv, false);
    if (eckit::mpi::comm().size() > 1) {
        if (barrier_timeout(ATLAS_MPI_BARRIER_TIMEOUT())) {
            eckit::Log::warning() << "\nWARNING: Tests failed with MPI deadlock "
                                     "(${ATLAS_MPI_BARRIER_TIMEOUT}="
                                  << ATLAS_MPI_BARRIER_TIMEOUT() << ").\nCalling MPI_Abort..." << std::endl;
            eckit::mpi::comm().abort();
        }
        eckit::mpi::comm().allReduceInPlace(errors, eckit::mpi::max());
    }
    return errors;
}

int run(int argc, char* argv[]) {
    return run<atlas::test::AtlasTestEnvironment>(argc, argv);
}

//----------------------------------------------------------------------------------------------------------------------

}  // namespace test
}  // namespace atlas
