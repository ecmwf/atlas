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
#include "eckit/log/PrefixTarget.h"
#include "eckit/runtime/Main.h"
#include "eckit/testing/Test.h"
#include "eckit/types/Types.h"

#include "atlas_io/atlas-io.h"
#include "atlas_io/Trace.h"

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
    ATLAS_IO_ASSERT(current_test_);
    return *current_test_;
}

//----------------------------------------------------------------------------------------------------------------------

#ifdef MAYBE_UNUSED
#elif defined(__GNUC__)
#define MAYBE_UNUSED __attribute__((unused))
#else
#define MAYBE_UNUSED
#endif

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
//bool approx_eq(const Point2& v1, const Point2& v2) {
//    return approx_eq(v1[0], v2[0]) && approx_eq(v1[1], v2[1]);
//}
//bool approx_eq(const Point2& v1, const Point2& v2, const double& t) {
//    return approx_eq(v1[0], v2[0], t) && approx_eq(v1[1], v2[1], t);
//}

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

#define EXPECT_APPROX_EQ(...) __ATLAS_IO_SPLICE(__EXPECT_APPROX_EQ__, __ATLAS_IO_NARG(__VA_ARGS__))(__VA_ARGS__)
#define __EXPECT_APPROX_EQ__2 __EXPECT_APPROX_EQ
#define __EXPECT_APPROX_EQ__3 __EXPECT_APPROX_EQ_TOL


//----------------------------------------------------------------------------------------------------------------------


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

struct TestEnvironment {
    TestEnvironment(int argc, char* argv[]) { eckit::Main::initialise(argc, argv); }

    ~TestEnvironment() {}
};


//----------------------------------------------------------------------------------------------------------------------


template <typename Environment>
int run(int argc, char* argv[]) {
    Environment env(argc, argv);
    int errors = eckit::testing::run_tests(argc, argv, false);
    return errors;
}

int run(int argc, char* argv[]) {
    return run<atlas::test::TestEnvironment>(argc, argv);
}

//----------------------------------------------------------------------------------------------------------------------

}  // namespace test
}  // namespace atlas
