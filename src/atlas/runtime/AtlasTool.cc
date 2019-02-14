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

#include "eckit/config/LibEcKit.h"
#include "eckit/log/FileTarget.h"
#include "eckit/log/PrefixTarget.h"

#include "atlas/library/config.h"
#include "atlas/runtime/AtlasTool.h"

namespace atlas {
static void usage( const std::string& name ) {
    Log::info() << "dummy usage" << std::endl;
}
}  // namespace atlas

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
    if ( eckit::Log::debug() ) eckit::Log::debug().addTarget( target );
}

void debug_setTarget( eckit::LogTarget* target ) {
    for ( std::string libname : eckit::system::Library::list() ) {
        const eckit::system::Library& lib = eckit::system::Library::lookup( libname );
        if ( lib.debug() ) {
            lib.debugChannel().setTarget( new eckit::PrefixTarget( debug_prefix( libname ), target ) );
        }
    }
    if ( eckit::Log::debug() ) eckit::Log::debug().setTarget( target );
}

void debug_reset() {
    for ( std::string libname : eckit::system::Library::list() ) {
        const eckit::system::Library& lib = eckit::system::Library::lookup( libname );
        if ( lib.debug() ) { lib.debugChannel().reset(); }
    }
    if ( eckit::Log::debug() ) eckit::Log::debug().reset();
}

bool getEnv( const std::string& env, bool default_value ) {
    if (::getenv( env.c_str() ) ) { return eckit::Translator<std::string, bool>()(::getenv( env.c_str() ) ); }
    return default_value;
}

int getEnv( const std::string& env, int default_value ) {
    if (::getenv( env.c_str() ) ) { return eckit::Translator<std::string, int>()(::getenv( env.c_str() ) ); }
    return default_value;
}

}  // namespace

void atlas::AtlasTool::add_option( eckit::option::Option* option ) {
    options_.push_back( option );
}

void atlas::AtlasTool::help( std::ostream& out ) {
    out << "NAME\n" << indent() << name();
    std::string brief = briefDescription();
    if ( brief.size() ) out << " - " << brief << '\n';

    std::string usg = usage();
    if ( usg.size() ) {
        out << '\n';
        out << "SYNOPSIS\n" << indent() << usg << '\n';
    }
    std::string desc = longDescription();
    if ( desc.size() ) {
        out << '\n';
        out << "DESCRIPTION\n" << indent() << desc << '\n';
    }
    out << '\n';
    out << "OPTIONS\n";
    for ( Options::const_iterator it = options_.begin(); it != options_.end(); ++it ) {
        out << indent() << ( **it ) << "\n\n";
    }
    out << std::flush;
}

bool atlas::AtlasTool::handle_help() {
    for ( int i = 1; i < argc(); ++i ) {
        if ( argv( i ) == "--help" || argv( i ) == "-h" ) {
            if ( taskID() == 0 ) help( std::cout );
            return true;
        }
    }
    return false;
}

atlas::AtlasTool::AtlasTool( int argc, char** argv ) : eckit::Tool( argc, argv ) {
    eckit::LibEcKit::instance().setAbortHandler( [] {
        Log::error() << "[" << atlas::mpi::comm().rank() << "] "
                     << "calling MPI_Abort" << std::endl;
        atlas::mpi::comm().abort( 1 );
    } );

    add_option( new SimpleOption<bool>( "help", "Print this help" ) );
    add_option( new SimpleOption<long>( "debug", "Debug level" ) );
    taskID( eckit::mpi::comm( "world" ).rank() );
}

int atlas::AtlasTool::start() {
    int status = 0;
    try {
        run();
    }
    catch ( eckit::Exception& e ) {
        status = 1;
        Log::error() << "** " << e.what() << e.location() << std::endl;
        Log::error() << "** Backtrace:\n" << e.callStack() << '\n';
        Log::error() << "** Exception  caught in " << Here() << " terminates " << name() << std::endl;
    }
    catch ( std::exception& e ) {
        status = 1;
        Log::error() << "** " << e.what() << " caught in " << Here() << '\n';
        Log::error() << "** Exception terminates " << name() << std::endl;
    }
    catch ( ... ) {
        status = 1;
        Log::error() << "** Exception caught in " << Here() << '\n';
        Log::error() << "** Exception terminates " << name() << std::endl;
    }
    if ( status ) {
        Log::error() << std::flush;
        eckit::LibEcKit::instance().abort();
    }

    return status;
}

void atlas::AtlasTool::run() {
    if ( handle_help() ) return;

    if ( argc() - 1 < minimumPositionalArguments() ) {
        if ( taskID() == 0 ) std::cout << "Usage: " << usage() << std::endl;
        return;
    }
    Options opts = options_;
    Args args( &atlas::usage, opts, numberOfPositionalArguments(), minimumPositionalArguments() );

    atlas::Library::instance().initialise();
    setupLogging();
    execute( args );
    atlas::Library::instance().finalise();
}

void atlas::AtlasTool::setupLogging() {
    int log_rank     = getEnv( "ATLAS_LOG_RANK", 0 );
    bool use_logfile = getEnv( "ATLAS_LOG_FILE", false );

    if ( use_logfile ) {
        int d               = digits( mpi::comm().size() );
        std::string rankstr = std::to_string( taskID() );
        for ( int i = rankstr.size(); i < d; ++i )
            rankstr = "0" + rankstr;

        eckit::LogTarget* logfile = new eckit::FileTarget( displayName() + ".log.p" + rankstr );

        if ( int( mpi::comm().rank() ) == log_rank ) {
            if ( Log::info() ) Log::info().addTarget( logfile );
            if ( Log::warning() ) Log::warning().addTarget( logfile );
            if ( Log::error() ) Log::error().addTarget( logfile );
            debug_addTarget( logfile );
        }
        else {
            if ( Log::info() ) Log::info().setTarget( logfile );
            if ( Log::warning() ) Log::warning().setTarget( logfile );
            if ( Log::error() ) Log::error().setTarget( logfile );
            debug_setTarget( logfile );
        }
    }
    else {
        if ( int( mpi::comm().rank() ) != log_rank ) {
            if ( Log::info() ) Log::info().reset();
            if ( Log::warning() ) Log::warning().reset();
            if ( Log::error() ) Log::error().reset();
            debug_reset();
        }
    }
}
