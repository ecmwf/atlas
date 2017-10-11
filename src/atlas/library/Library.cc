/*
 * (C) Copyright 1996-2017 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */



#include "atlas/library/config.h"

// Temporary until ECKIT-166 is fixed
#ifdef BUG_ECKIT_166
#include <mpi.h>
#endif

#ifdef ATLAS_HAVE_TRANS
#include "transi/version.h"
#endif

#include "eckit/runtime/Main.h"
#include "eckit/log/Log.h"
#include "eckit/filesystem/PathName.h"
#include "eckit/filesystem/LocalPathName.h"
#include "eckit/utils/Translator.h"
#include "eckit/log/PrefixTarget.h"
#include "eckit/log/OStreamTarget.h"

#include "atlas/parallel/mpi/mpi.h"
#include "atlas/runtime/Log.h"
#include "atlas/util/Config.h"
#include "atlas/library/Library.h"
#include "atlas/library/version.h"
#include "atlas/library/git_sha1.h"
#include "atlas/runtime/Trace.h"

using eckit::PathName;
using eckit::Main;
using eckit::LocalPathName;

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
    if( lib.gitsha1() != "not available" ) {
      ss << "  git-sha1 " << lib.gitsha1(7);
    }
    return ss.str();
  }
}

//----------------------------------------------------------------------------------------------------------------------

static Library libatlas;


Library::Library() : eckit::system::Library( std::string("atlas") ) {}

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

void Library::initialise(int argc, char **argv) {
    if( not Main::ready() ) {
        Main::initialise(argc, argv);
        Main::instance().taskID( eckit::mpi::comm("world").rank() );
        if( Main::instance().taskID() != 0 ) Log::reset();
        Log::debug() << "Atlas initialised eckit::Main.\n";
        if( eckit::mpi::comm("world").size() > 1 )
            Log::debug() <<
              "--> Only MPI rank 0 is logging. Please initialise eckit::Main \n"
              "    before to avoid this behaviour.\n";
    }
    initialise();
}

namespace{
  bool getEnv( const std::string& env, bool default_value ) {
      if (::getenv( env.c_str() ) ) {
        return eckit::Translator<std::string, bool>()(::getenv( env.c_str() ));
      }
      return default_value;
  }
}
void Library::initialise(const eckit::Parametrisation& config) {

    if( not config.get("trace",info_) ) {
      info_ = getEnv("ATLAS_INFO",info_);
    }
    if( not config.get("trace",trace_) ) {
      trace_ = getEnv("ATLAS_TRACE",trace_);
    }
    if( not config.get("barriers",barriers_) ) {
      barriers_ = getEnv("ATLAS_BARRIERS",barriers_);
    }

    // Summary
    std::ostream& out = Log::debug();
    out << "Executable        [" << Main::instance().name() << "]\n";
    out << " \n";
    out << "  current dir     [" << PathName(LocalPathName::cwd()).fullName() << "]\n";
    out << " \n";
    out << "  MPI\n";
    out << "    communicator  [" << parallel::mpi::comm() << "] \n";
    out << "    size          [" << parallel::mpi::comm().size() << "] \n";
    out << "    rank          [" << parallel::mpi::comm().rank() << "] \n";
    out << "    barriers      [" << barriers() << "] \n";
    out << " \n";
    out << "  trace           [" << str(trace()) << "] \n";
    out << "  debug           [" << str(debug()) << "] \n";
    out << " \n";
    out << atlas::Library::instance().information();
    out << std::flush;
}

void Library::initialise() {
  initialise( util::NoConfig() );
}

void Library::finalise() {
// Temporary until ECKIT-166 is fixed
#ifdef BUG_ECKIT_166
    const bool using_mpi = (::getenv("OMPI_COMM_WORLD_SIZE") || ::getenv("ALPS_APP_PE"));
    if( using_mpi ) {
      int finalized = 1;
      MPI_Finalized(&finalized);
      if( not finalized ) {
        MPI_Finalize();
      }
    }
#endif

    // Make sure that these specialised channels that wrap Log::info() are 
    // destroyed before Log::info gets destroyed.
    // Just in case someone still tries to log, we reset to empty channels.
    trace_channel_.reset( new eckit::Channel() );

    Log::debug() << "Atlas finalised" << std::endl;

    Log::flush();
}

std::ostream& Library::traceChannel() const {
  if( trace_channel_ ) return *trace_channel_;
  if( trace_ ) {
    trace_channel_.reset( new eckit::Channel(
      new eckit::PrefixTarget("ATLAS_TRACE", new eckit::OStreamTarget(eckit::Log::info())))); 
  } else {
    trace_channel_.reset( new eckit::Channel() );
  }
  return *trace_channel_;
}

std::ostream& Library::infoChannel() const {
  if( info_channel_ ) return *info_channel_;
  if( info_ ) {
    return eckit::Log::info();
  } else {
    info_channel_.reset( new eckit::Channel() );
  }
  return *trace_channel_;
}


//----------------------------------------------------------------------------------------------------------------------

void Library::Information::print( std::ostream& out ) const {
    out << "atlas version (" << atlas::Library::instance().version() << "), "
        << "git-sha1 "<< atlas::Library::instance().gitsha1(7) << '\n';
    out << " \n";
    out << "  Build:" << std::endl;
    out << "    build type      : " << ATLAS_BUILD_TYPE << '\n'
        << "    timestamp       : " << ATLAS_BUILD_TIMESTAMP << '\n'
        << "    op. system      : " << ATLAS_OS_NAME << " (" << ATLAS_OS_STR << ")"  << '\n'
        << "    processor       : " << ATLAS_SYS_PROCESSOR  << std::endl
        << "    c compiler      : " << ATLAS_C_COMPILER_ID << " " << ATLAS_C_COMPILER_VERSION << '\n'
        << "      flags         : " << ATLAS_C_FLAGS << '\n'
        << "    c++ compiler    : " << ATLAS_CXX_COMPILER_ID << " " << ATLAS_CXX_COMPILER_VERSION << '\n'
        << "      flags         : " << ATLAS_CXX_FLAGS << '\n'
#ifndef EC_HAVE_FORTRAN
        << "    fortran         : NO " << '\n'
#else
       << "    fortran compiler: " << ATLAS_Fortran_COMPILER_ID << " " << ATLAS_Fortran_COMPILER_VERSION << '\n'
       << "      flags         : " << ATLAS_Fortran_FLAGS << '\n'
#endif
       << " \n";

    bool feature_fortran(false);
    bool feature_OpenMP(false);
    bool feature_Trans(false);
    bool feature_Tesselation(false);
    bool feature_BoundsChecking(false);
#ifdef ATLAS_HAVE_FORTRAN
      feature_fortran = true;
#endif
#ifdef ATLAS_HAVE_OMP
      feature_OpenMP = true;
#endif
#ifdef ATLAS_HAVE_TRANS
      feature_Trans = true;
#endif
#ifdef ATLAS_HAVE_TESSELATION
      feature_Tesselation = true;
#endif
#ifdef ATLAS_HAVE_BOUNDSCHECKING
      feature_BoundsChecking = true;
#endif
      std::string array_data_store = "Native";
#ifdef ATLAS_HAVE_GRIDTOOLS_STORAGE
#if GRIDTOOLS_STORAGE_BACKEND_CUDA
      array_data_store = "GridTools-CUDA";
#endif
#if ATLAS_GRIDTOOLS_STORAGE_BACKEND_HOST
      array_data_store = "GridTools-host";
#endif
#endif
      out << "  Features:" << '\n'
        << "    Fortran        : " << str(feature_fortran) << '\n'
        << "    OpenMP         : " << str(feature_OpenMP) << '\n'
        << "    BoundsChecking : " << str(feature_BoundsChecking) << '\n'
        << "    Trans          : " << str(feature_Trans) << '\n'
        << "    Tesselation    : " << str(feature_Tesselation) << '\n'
        << "    ArrayDataStore : " << array_data_store << '\n'
        << "    gidx_t         : " << ATLAS_BITS_GLOBAL << " bit integer" << '\n'
        << " \n";

    out << "  Dependencies: " << "\n";

    if( Library::exists("eckit") ) {
      out << "    " << str( Library::lookup("eckit") ) << '\n';
    }
    if( Library::exists("fckit") ) {
      out << "    " << str( Library::lookup("fckit") ) << '\n';
    }

#ifdef ATLAS_HAVE_TRANS
    out << "    transi version (" << transi_version() << "), "
                << "git-sha1 "<< transi_git_sha1_abbrev(7) << '\n';
    out << "    trans version (" << trans_version() << "), "
                << "git-sha1 "<< trans_git_sha1_abbrev(7) << '\n';
#endif

}

} // namespace atlas

