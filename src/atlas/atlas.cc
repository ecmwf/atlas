#include <stdlib.h>
#include "eckit/runtime/Context.h"
#include "eckit/parser/StringTools.h"
#include "eckit/filesystem/PathName.h"
#include "eckit/filesystem/LocalPathName.h"
#include "atlas/atlas.h"
#include "atlas/parallel/mpi/mpi.h"
//#include "atlas/grid/grids.h"
#include "atlas/mesh/generators/MeshGenerator.h"
#include "atlas/field/FieldCreator.h"
#include "atlas/runtime/Log.h"
#include "atlas/runtime/Behavior.h"
#include "atlas/util/Config.h"

using namespace eckit;

namespace atlas {

std::string rundir()
{
  static PathName cwd( LocalPathName::cwd() );
  return cwd;
}

// ----------------------------------------------------------------------------

void atlas_info( std::ostream& out )
{
  out << "Atlas initialised [" << Context::instance().displayName() << "]\n";
  out << runtime::indent();
  out << "atlas version [" << atlas_version() << "]\n";
  out << "atlas git     [" << atlas_git_sha1()<< "]\n";
  out << "eckit version [" << eckit_version() << "]\n";
  out << "eckit git     [" << eckit_git_sha1()<< "]\n";
  Context::instance().behavior().reconfigure();
  out << "current dir   [" << PathName(rundir()).fullName() << "]\n";
  if( eckit::mpi::initialized() ) {
    out << "MPI\n" << runtime::indent();
    out << "communicator  [" << eckit::mpi::comm().fortran_communicator() << "] \n";
    out << "size          [" << eckit::mpi::comm().size() << "] \n";
    out << "rank          [" << eckit::mpi::comm().rank() << "] \n";
    out << runtime::dedent();
  }
  else
  {
    out << "MPI           [OFF]\n";
  }
  out << runtime::dedent();
}

// ----------------------------------------------------------------------------

void atlas_init()
{
  return atlas_init(util::NoConfig());
}

// ----------------------------------------------------------------------------

void atlas_init(const eckit::Parametrisation&)
{
  atlas_info( Log::debug() );

  // Load factories for static linking
  //atlas::grid::load();
}

// ----------------------------------------------------------------------------

// This is only to be used from Fortran or unit-tests
void atlas_init(int argc, char** argv)
{
  if( Context::instance().argc() == 0 )
  {
    if( argc>0 )
      Context::instance().setup(argc, argv);
    if( Context::instance().argc() > 0 )
      Context::instance().runName( PathName(Context::instance().argv(0)).baseName(false) );

    long debug(0);
    const char* env_debug = ::getenv("DEBUG");
    if( env_debug ) debug = ::atol(env_debug);
    // args.get("debug",debug);

    Context::instance().behavior( new atlas::runtime::Behavior() );
    Context::instance().debug(debug);
  }
  atlas_init();
}

// ----------------------------------------------------------------------------

void atlas_finalize()
{
  Log::debug() << "Atlas finalized\n" << std::flush;
}

// ----------------------------------------------------------------------------

void atlas__atlas_init(int argc, char* argv[])
{
  atlas_init(argc,argv);
}

// ----------------------------------------------------------------------------

void atlas__atlas_init_noargs()
{
  static int argc = 1;
  static char const *argv[] = {"atlas_program"};
  atlas_init(argc,(char**)argv);
}

// ----------------------------------------------------------------------------

void atlas__atlas_finalize()
{
  atlas_finalize();
}

// ----------------------------------------------------------------------------

const char* atlas__eckit_version()
{
  return eckit_version();
}

// ----------------------------------------------------------------------------

const char* atlas__eckit_git_sha1()
{
  return eckit_git_sha1();
}

// ----------------------------------------------------------------------------

const char* atlas__eckit_git_sha1_abbrev(int length)
{
  static std::string git_sha1(eckit_git_sha1());
  if( git_sha1.empty() ) git_sha1 = "not available";
  else                   git_sha1 = git_sha1.substr(0,std::min(length,40));
  return git_sha1.c_str();
}

// ----------------------------------------------------------------------------

const char* atlas__atlas_version()
{
  return atlas_version();
}

// ----------------------------------------------------------------------------

const char* atlas__atlas_git_sha1()
{
  return atlas_git_sha1();
}

// ----------------------------------------------------------------------------

const char* atlas__atlas_git_sha1_abbrev(int length)
{
  return atlas_git_sha1_abbrev(length);
}

// ----------------------------------------------------------------------------

const char* atlas__run_name ()
{
  static std::string str( Context::instance().runName() );
  return str.c_str();
}

// ----------------------------------------------------------------------------

const char* atlas__display_name ()
{
  static std::string str( Context::instance().displayName() );
  return str.c_str();
}

// ----------------------------------------------------------------------------

const char* atlas__rundir ()
{
  return rundir().c_str();
}

// ----------------------------------------------------------------------------

const char* atlas__workdir ()
{
  static LocalPathName workdir = LocalPathName::cwd().fullName();
  return workdir.c_str();
}

// ----------------------------------------------------------------------------

} // namespace atlas


