#include <stdlib.h>

// Temporary until ECKIT-166 is fixed
#ifdef BUG_ECKIT_166
#include <mpi.h>
#endif

#include "atlas/internals/atlas_defines.h"
#include "eckit/runtime/Main.h"
#include "eckit/parser/StringTools.h"
#include "eckit/filesystem/PathName.h"
#include "eckit/filesystem/LocalPathName.h"
#include "atlas/parallel/mpi/mpi.h"

#include "atlas/atlas.h"
#include "atlas/mesh/generators/MeshGenerator.h"
#include "atlas/field/FieldCreator.h"
#include "atlas/runtime/Log.h"
#include "atlas/util/Config.h"

using eckit::PathName;
using eckit::Main;
using eckit::LocalPathName;

namespace atlas {

std::string rundir()
{
  static PathName cwd( LocalPathName::cwd() );
  return cwd;
}

//----------------------------------------------------------------------------------------------------------------------

void atlas_info( std::ostream& out )
{
  out << "Executable        [" << Main::instance().name() << "]\n";

  out << "  atlas version   [" << atlas_version() << "]\n";
  out << "  atlas git       [" << atlas_git_sha1_abbrev(7)<< "]\n";
  out << "  eckit version   [" << eckit_version() << "]\n";
  out << "  eckit git       [" << atlas__eckit_git_sha1_abbrev(7)<< "]\n";
  out << "  current dir     [" << PathName(rundir()).fullName() << "]\n";

  out << "  MPI\n";
  out << "    communicator  [" << parallel::mpi::comm() << "] \n";
  out << "    size          [" << parallel::mpi::comm().size() << "] \n";
  out << "    rank          [" << parallel::mpi::comm().rank() << "] \n";
  out << std::flush;
}

//----------------------------------------------------------------------------------------------------------------------


void atlas_init(const eckit::Parametrisation&)
{
  atlas_info( Log::debug() );
}
void atlas_init() { return atlas_init(util::NoConfig()); }


/// deprecated
void atlas_init(int argc, char** argv)
{
  Main::initialise(argc, argv);
  Main::instance().taskID(eckit::mpi::comm("world").rank());
  if( Main::instance().taskID() != 0 ) Log::reset();
  atlas_init();
}

void atlas_finalize()
{
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

  Log::debug() << "Atlas finalized\n";
  Log::flush();
}

//----------------------------------------------------------------------------------------------------------------------

void atlas__atlas_init_noargs()
{
  atlas_init();
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

void atlas__info(std::ostream* channel)
{
  atlas_info(*channel);
}


//----------------------------------------------------------------------------------------------------------------------

} // namespace atlas


