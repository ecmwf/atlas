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

#include <chrono>
#include <exception>
#include <thread>

#include "eckit/config/LibEcKit.h"
#include "eckit/config/Resource.h"
#include "eckit/eckit_config.h"
#include "eckit/log/PrefixTarget.h"
#include "eckit/mpi/Comm.h"
#include "eckit/runtime/Main.h"
#include "eckit/testing/Test.h"
#include "eckit/types/Types.h"

#include "atlas/library/Library.h"
#include "atlas/runtime/Log.h"
#include "atlas/runtime/Trace.h"
#include "atlas/runtime/trace/StopWatch.h"
#include "atlas/util/Config.h"

namespace atlas {
namespace test {

using eckit::types::is_approximately_equal;

//----------------------------------------------------------------------------------------------------------------------

#ifdef MAYBE_UNUSED
#elif defined( __GNUC__ )
#define MAYBE_UNUSED __attribute__( ( unused ) )
#elif
#define MAYBE_UNUSED
#endif

// Redefine macro's defined in "eckit/testing/Test.h" to include trace
// information
#undef CASE
#define CASE( description )                                                                                          \
    void UNIQUE_NAME2( test_, __LINE__ )( std::string&, int&, int );                                                 \
    static eckit::testing::TestRegister UNIQUE_NAME2( test_registration_, __LINE__ )(                                \
        description, &UNIQUE_NAME2( test_, __LINE__ ) );                                                             \
    void UNIQUE_NAME2( traced_test_, __LINE__ )( std::string & _test_subsection, int& _num_subsections,              \
                                                 int _subsection );                                                  \
    void UNIQUE_NAME2( test_, __LINE__ )( std::string & _test_subsection, int& _num_subsections, int _subsection ) { \
        ATLAS_TRACE( description );                                                                                  \
        UNIQUE_NAME2( traced_test_, __LINE__ )                                                                       \
        ( _test_subsection, _num_subsections, _subsection );                                                         \
        if ( atlas::test::barrier_timeout( atlas::test::ATLAS_MPI_BARRIER_TIMEOUT() ) ) {                            \
            atlas::Log::warning() << "\nWARNING: Test \"" << description                                             \
                                  << "\" failed with MPI deadlock.  (${ATLAS_MPI_BARRIER_TIMEOUT}="                  \
                                  << atlas::test::ATLAS_MPI_BARRIER_TIMEOUT() << ").\nCalling MPI_Abort..."          \
                                  << std::endl;                                                                      \
            eckit::mpi::comm().abort();                                                                              \
        }                                                                                                            \
    }                                                                                                                \
    void UNIQUE_NAME2( traced_test_, __LINE__ )( MAYBE_UNUSED std::string & _test_subsection,                        \
                                                 MAYBE_UNUSED int& _num_subsections, MAYBE_UNUSED int _subsection )

#undef SECTION
#define SECTION( name )                                                                          \
    _num_subsections += 1;                                                                       \
    _test_subsection = ( name );                                                                 \
    if ( ( _num_subsections - 1 ) == _subsection ) {                                             \
        atlas::Log::info() << "Running section \"" << _test_subsection << "\" ..." << std::endl; \
    }                                                                                            \
    if ( ( _num_subsections - 1 ) == _subsection ) ATLAS_TRACE_SCOPE( _test_subsection )

#ifndef SETUP
#define SETUP( name )
#endif

//----------------------------------------------------------------------------------------------------------------------

static double ATLAS_MPI_BARRIER_TIMEOUT() {
    static int v = eckit::Resource<double>( "${ATLAS_MPI_BARRIER_TIMEOUT", 3. );
    return v;
}

static int barrier_timeout( double seconds ) {
    auto req = eckit::mpi::comm().iBarrier();
    runtime::trace::StopWatch watch;
    while ( not req.test() ) {
        watch.start();
        std::this_thread::sleep_for( std::chrono::milliseconds( 100 ) );
        watch.stop();
        if ( watch.elapsed() > seconds ) { return 1; }
    }
    return 0;
}

struct AtlasTestEnvironment {
    using Config = util::Config;

    AtlasTestEnvironment( int argc, char* argv[] ) {
        eckit::Main::initialise( argc, argv );
        eckit::Main::instance().taskID( eckit::mpi::comm( "world" ).rank() );
        if ( eckit::mpi::comm( "world" ).size() != 1 ) {
            long logtask = eckit::Resource<long>( "$ATLAS_LOG_TASK", 0 );
            if ( eckit::Main::instance().taskID() != logtask ) {
                eckit::Log::info().reset();
                eckit::Log::warning().reset();
                eckit::Log::debug().reset();
            }
            eckit::Log::error().setTarget(
                new eckit::PrefixTarget( "[" + std::to_string( eckit::mpi::comm().rank() ) + "]" ) );
        }
        eckit::LibEcKit::instance().setAbortHandler( [] {
            Log::error() << "Calling MPI_Abort" << std::endl;
            eckit::mpi::comm().abort();
        } );
        Library::instance().initialise();
        eckit::mpi::comm().barrier();
    }

    ~AtlasTestEnvironment() { Library::instance().finalise(); }
};

//----------------------------------------------------------------------------------------------------------------------

template <typename Environment>
int run( int argc, char* argv[] ) {
    Environment env( argc, argv );
    int errors = eckit::testing::run_tests( argc, argv, false );
    if ( eckit::mpi::comm().size() > 1 ) {
        if ( barrier_timeout( ATLAS_MPI_BARRIER_TIMEOUT() ) ) {
            eckit::Log::warning() << "\nWARNING: Tests failed with MPI deadlock "
                                     "(${ATLAS_MPI_BARRIER_TIMEOUT}="
                                  << ATLAS_MPI_BARRIER_TIMEOUT() << ").\nCalling MPI_Abort..." << std::endl;
            eckit::mpi::comm().abort();
        }
        eckit::mpi::comm().allReduceInPlace( errors, eckit::mpi::max() );
    }
    return errors;
}

int run( int argc, char* argv[] ) {
    return run<atlas::test::AtlasTestEnvironment>( argc, argv );
}

//----------------------------------------------------------------------------------------------------------------------

}  // end namespace test
}  // end namespace atlas
