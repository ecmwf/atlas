#include <stdlib.h>

#include "eckit/runtime/Main.h"
#include "eckit/parser/StringTools.h"
#include "eckit/filesystem/PathName.h"
#include "eckit/filesystem/LocalPathName.h"
#include "eckit/mpi/Comm.h"

#include "atlas/atlas.h"
#include "atlas/grid/grids.h"
#include "atlas/mesh/generators/MeshGenerator.h"
#include "atlas/field/FieldCreator.h"
#include "atlas/runtime/Log.h"
#include "atlas/util/Config.h"

using eckit::PathName;
using eckit::Main;
using eckit::LocalPathName;

namespace atlas {

class Environment {
private:
  Environment()
  {
    eckit::mpi::init( Main::instance().argc(), Main::instance().argv() );
  }

  ~Environment()
  {
    eckit::mpi::finalize();
  }

public:
  static Environment& setupMPI()
  {
    static Environment env;
    return env;
  }

  static bool usingMPI() {
      return (::getenv("OMPI_COMM_WORLD_SIZE") || ::getenv("ALPS_APP_PE"));
  }

};

std::string rundir()
{
  static PathName cwd( LocalPathName::cwd() );
  return cwd;
}

//----------------------------------------------------------------------------------------------------------------------

void atlas_info( std::ostream& out )
{
  out << "Atlas initialised [" << Main::instance().name() << "]\n";

  out << "atlas version [" << atlas_version() << "]\n";
  out << "atlas git     [" << atlas_git_sha1()<< "]\n";
  out << "eckit version [" << eckit_version() << "]\n";
  out << "eckit git     [" << eckit_git_sha1()<< "]\n";
  out << "current dir   [" << PathName(rundir()).fullName() << "]\n";

  if(eckit::mpi::comm().size() > 1) {
    out << "MPI\n";
    out << "communicator  [" << eckit::mpi::comm() << "] \n";
    out << "size          [" << eckit::mpi::comm().size() << "] \n";
    out << "rank          [" << eckit::mpi::comm().rank() << "] \n";
  }
  else
  {
    out << "MPI           [OFF]\n";
  }
}

//----------------------------------------------------------------------------------------------------------------------

void atlas_init()
{
  return atlas_init(util::NoConfig());
}


void atlas_init(const eckit::Parametrisation&)
{
  if( eckit::mpi::comm().rank() > 0 ) {
    Log::info().reset();
    Log::debug().reset();
    Log::warning().reset();
  }

  atlas_info( Log::debug() );

  // Load factories for static linking
  atlas::grid::load();
}


// This is only to be used from Fortran or unit-tests
void atlas_init(int argc, char** argv)
{
  Main::initialise(argc, argv);
  atlas_init();
}


void atlas_finalize()
{
  Log::debug() << "Atlas finalized\n" << std::flush;
}

//----------------------------------------------------------------------------------------------------------------------

void atlas__atlas_init(int argc, char* argv[])
{
  atlas_init(argc,argv);
}


void atlas__atlas_init_noargs()
{
  static int argc = 1;
  static char const *argv[] = {"atlas_program"};
  atlas_init(argc,(char**)argv);
}


void atlas__atlas_finalize()
{
  atlas_finalize();
}


const char* atlas__eckit_version()
{
  return eckit_version();
}


const char* atlas__eckit_git_sha1()
{
  return eckit_git_sha1();
}


const char* atlas__eckit_git_sha1_abbrev(int length)
{
  static std::string git_sha1(eckit_git_sha1());
  if( git_sha1.empty() ) git_sha1 = "not available";
  else                   git_sha1 = git_sha1.substr(0,std::min(length,40));
  return git_sha1.c_str();
}


const char* atlas__atlas_version()
{
  return atlas_version();
}


const char* atlas__atlas_git_sha1()
{
  return atlas_git_sha1();
}


const char* atlas__atlas_git_sha1_abbrev(int length)
{
  return atlas_git_sha1_abbrev(length);
}


const char* atlas__run_name ()
{
  static std::string str( Main::instance().name() );
  return str.c_str();
}


const char* atlas__display_name ()
{
  static std::string str( Main::instance().name() );
  return str.c_str();
}


const char* atlas__rundir ()
{
  return rundir().c_str();
}


const char* atlas__workdir ()
{
  static LocalPathName workdir = LocalPathName::cwd().fullName();
  return workdir.c_str();
}

//----------------------------------------------------------------------------------------------------------------------

} // namespace atlas


