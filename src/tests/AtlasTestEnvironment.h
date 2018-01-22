/*
 * (C) Copyright 1996-2017 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#pragma once

#include "atlas/library/Library.h"
#include "eckit/runtime/Main.h"
#include "eckit/mpi/Comm.h"
#include "atlas/runtime/Log.h"
#include "atlas/runtime/Trace.h"
#include "eckit/config/Resource.h"
#include "eckit/testing/Test.h"

namespace atlas {
namespace test {

//----------------------------------------------------------------------------------------------------------------------

// Redefine macro's defined in "eckit/testing/Test.h" to include trace information

#ifdef CASE
#undef CASE
#endif

#define CASE(description) \
void UNIQUE_NAME2(test_, __LINE__) (std::string& _test_subsection); \
static eckit::testing::TestRegister UNIQUE_NAME2(test_registration_, __LINE__)(description, &UNIQUE_NAME2(test_, __LINE__)); \
void UNIQUE_NAME2(traced_test_, __LINE__)(std::string& _test_subsection); \
void UNIQUE_NAME2(test_, __LINE__) (std::string& _test_subsection) { \
    ATLAS_TRACE(description); \
    UNIQUE_NAME2(traced_test_, __LINE__)(_test_subsection); \
} \
void UNIQUE_NAME2(traced_test_, __LINE__)(std::string& _test_subsection)

#ifdef SECTION
#undef SECTION
#endif
#define SECTION(name) \
    _test_num += 1; \
    _test_count = _test_num; \
    _test_subsection = name; \
    if ((_test_num - 1) == _test) ATLAS_TRACE_SCOPE(name)

//----------------------------------------------------------------------------------------------------------------------

struct AtlasTestEnvironment {

    AtlasTestEnvironment(int argc, char * argv[]) {
        eckit::Main::initialise(argc, argv);
        eckit::Main::instance().taskID( eckit::mpi::comm("world").rank() );
        if( eckit::mpi::comm("world").size() != 1 ) {
            long logtask = eckit::Resource<long>("--logtask;$ATLAS_TEST_LOGTASK", 0);
            if( eckit::Main::instance().taskID() != logtask ) Log::reset();
        }
        atlas::Library::instance().initialise();
    }

    ~AtlasTestEnvironment() {
        bool report = eckit::Resource<bool>("--report;$ATLAS_TEST_REPORT", false);
        if( report ) {
#if ATLAS_HAVE_TRACE
          Log::info() << atlas::Trace::report() << std::endl;
#else
          Log::warning() << "Atlas cannot generate report as ATLAS_HAVE_TRACE is not defined." << std::endl;
#endif
        }
        atlas::Library::instance().finalise();
    }
};

//----------------------------------------------------------------------------------------------------------------------

template< typename Environment >
int run(int argc, char* argv[]) {
    Environment env( argc, argv );
    return eckit::testing::run_tests(argc,argv,false);
}

int run(int argc, char* argv[]) {
    return run<atlas::test::AtlasTestEnvironment>( argc, argv );
}

//----------------------------------------------------------------------------------------------------------------------

}  // end namespace test
}  // end namespace atlas
