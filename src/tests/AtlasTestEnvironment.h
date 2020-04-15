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

namespace atlas {
namespace test {

using eckit::types::is_approximately_equal;

//----------------------------------------------------------------------------------------------------------------------

#ifdef MAYBE_UNUSED
#elif defined( __GNUC__ )
#define MAYBE_UNUSED __attribute__( ( unused ) )
#else
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
    if ( ( _num_subsections - 1 ) == _subsection )                                               \
    ATLAS_TRACE_SCOPE( _test_subsection )

#ifdef EXPECT_EQ
#undef EXPECT_EQ
#endif
#define EXPECT_EQ( lhs, rhs )                                                                                \
    do {                                                                                                     \
        if ( !( lhs == rhs ) ) {                                                                             \
            throw eckit::testing::TestException( "EXPECT condition failed: " #lhs " == " #rhs                \
                                                 "\n"                                                        \
                                                 " --> " +                                                   \
                                                     std::to_string( lhs ) + " != " + std::to_string( rhs ), \
                                                 Here() );                                                   \
        }                                                                                                    \
    } while ( false )

//----------------------------------------------------------------------------------------------------------------------

static double ATLAS_MPI_BARRIER_TIMEOUT() {
    static double v = eckit::Resource<double>( "${ATLAS_MPI_BARRIER_TIMEOUT", 3. );
    return v;
}

static int barrier_timeout( double seconds ) {
    auto req = eckit::mpi::comm().iBarrier();
    runtime::trace::StopWatch watch;
    while ( not req.test() ) {
        watch.start();
        std::this_thread::sleep_for( std::chrono::milliseconds( 100 ) );
        watch.stop();
        if ( watch.elapsed() > seconds ) {
            return 1;
        }
    }
    return 0;
}


namespace {

int digits( int number ) {
    int d = 0;
    while ( number ) {
        number /= 10;
        d++;
    }
    return d;
}

static std::string debug_prefix( const std::string& libname ) {
    std::string s = libname;
    std::transform( s.begin(), s.end(), s.begin(), ::toupper );
    s += "_DEBUG";
    return s;
}

void debug_addTarget( eckit::LogTarget* target ) {
    for ( std::string libname : eckit::system::Library::list() ) {
        const eckit::system::Library& lib = eckit::system::Library::lookup( libname );
        if ( lib.debug() ) {
            lib.debugChannel().addTarget( new eckit::PrefixTarget( debug_prefix( libname ), target ) );
        }
    }
    if ( eckit::Log::debug() )
        eckit::Log::debug().addTarget( target );
}

void debug_setTarget( eckit::LogTarget* target ) {
    for ( std::string libname : eckit::system::Library::list() ) {
        const eckit::system::Library& lib = eckit::system::Library::lookup( libname );
        if ( lib.debug() ) {
            lib.debugChannel().setTarget( new eckit::PrefixTarget( debug_prefix( libname ), target ) );
        }
    }
    if ( eckit::Log::debug() )
        eckit::Log::debug().setTarget( target );
}

void debug_reset() {
    for ( std::string libname : eckit::system::Library::list() ) {
        const eckit::system::Library& lib = eckit::system::Library::lookup( libname );
        if ( lib.debug() ) {
            lib.debugChannel().reset();
        }
    }
    if ( eckit::Log::debug() )
        eckit::Log::debug().reset();
}

bool getEnv( const std::string& env, bool default_value ) {
    if ( ::getenv( env.c_str() ) ) {
        return eckit::Translator<std::string, bool>()( ::getenv( env.c_str() ) );
    }
    return default_value;
}

int getEnv( const std::string& env, int default_value ) {
    if ( ::getenv( env.c_str() ) ) {
        return eckit::Translator<std::string, int>()( ::getenv( env.c_str() ) );
    }
    return default_value;
}

}  // namespace

struct AtlasTestEnvironment {
    using Config = util::Config;

    AtlasTestEnvironment( int argc, char* argv[] ) {
        eckit::Main::initialise( argc, argv );
        eckit::Main::instance().taskID( eckit::mpi::comm( "world" ).rank() );

        long log_rank    = getEnv( "ATLAS_LOG_RANK", 0 );
        bool use_logfile = getEnv( "ATLAS_LOG_FILE", false );

        auto rank_str = []() {
            int d           = digits( eckit::mpi::comm().size() );
            std::string str = std::to_string( eckit::Main::instance().taskID() );
            for ( int i = str.size(); i < d; ++i )
                str = "0" + str;
            return str;
        };

        if ( use_logfile ) {
            eckit::LogTarget* logfile =
                new eckit::FileTarget( eckit::Main::instance().displayName() + ".log.p" + rank_str() );

            if ( eckit::Main::instance().taskID() == log_rank ) {
                eckit::Log::info().addTarget( logfile );
                eckit::Log::warning().addTarget( logfile );
                eckit::Log::error().setTarget( logfile );
                debug_addTarget( logfile );
            }
            else {
                eckit::Log::info().setTarget( logfile );
                eckit::Log::warning().setTarget( logfile );
                eckit::Log::error().setTarget( logfile );
                debug_setTarget( logfile );
            }
        }
        else {
            if ( eckit::Main::instance().taskID() != log_rank ) {
                eckit::Log::info().reset();
                eckit::Log::warning().reset();
                debug_reset();
            }
            eckit::Log::error().reset();
        }
        Log::error().addTarget( new eckit::PrefixTarget( "[" + std::to_string( eckit::mpi::comm().rank() ) + "]" ) );

        eckit::LibEcKit::instance().setAbortHandler( [] {
            Log::error() << "Calling MPI_Abort" << std::endl;
            eckit::mpi::comm().abort();
        } );
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
