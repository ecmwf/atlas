#include "eckit/runtime/Context.h"
#include "eckit/parser/StringTools.h"
#include "eckit/config/Resource.h"
#include "eckit/config/ResourceMgr.h"
#include "eckit/filesystem/PathName.h"
#include "eckit/filesystem/LocalPathName.h"
#include "atlas/atlas.h"
#include "atlas/grid/grids.h"
#include "atlas/mesh/generators/MeshGenerator.h"
#include "atlas/field/FieldCreator.h"
#include "atlas/util/runtime/Log.h"
#include "atlas/util/runtime/Behavior.h"

using namespace eckit;

namespace atlas {

std::string rundir()
{
  static PathName cwd( LocalPathName::cwd() );
  return cwd;
}


void atlas_init(int argc, char** argv)
{
  if( argc > 0 )
    Context::instance().setup(argc, argv);

  if( argc > 0 )
    Context::instance().runName( PathName(Context::instance().argv(0)).baseName(false) );


  Context::instance().behavior( new atlas::util::runtime::Behavior() );
  Context::instance().behavior().debug( Resource<int>("$DEBUG;--debug",0) );

  Context::instance().displayName( Resource<std::string>("--name;$ATLAS_DISPLAYNAME","") );

  std::vector<PathName> config_files;

  PathName atlas_config_path ( Resource<std::string>("atlas.configfile;$ATLAS_CONFIGFILE;--atlas_conf","atlas.cfg") );
  if( ResourceMgr::instance().appendConfig( atlas_config_path ) )
    config_files.push_back( atlas_config_path );

  PathName runname_config_path (
    Resource<std::string>("$"+StringTools::upper(Context::instance().runName())+"_CONFIGFILE",
                          Context::instance().runName()+".cfg") );
  if( ResourceMgr::instance().appendConfig( runname_config_path ) )
    config_files.push_back( runname_config_path );

  PathName displayname_config_path(
    Resource<std::string>(
       "$"+StringTools::upper(Context::instance().displayName())+"_CONFIGFILE;--conf",
       Context::instance().displayName()+".cfg") );
  if( ResourceMgr::instance().appendConfig( displayname_config_path ) )
    config_files.push_back( displayname_config_path );

  Context::instance().behavior().debug( Resource<int>("debug;$DEBUG;--debug",0) );

  Log::debug() << "Atlas program [" << Context::instance().displayName() << "]\n";
  Log::debug() << util::runtime::indent();
  Log::debug() << "atlas version [" << atlas_version() << "]\n";
  Log::debug() << "atlas git     [" << atlas_git_sha1()<< "]\n";
  Log::debug() << "eckit version [" << eckit_version() << "]\n";
  Log::debug() << "eckit git     [" << eckit_git_sha1()<< "]\n";

  std::vector<PathName>::iterator path_it;
  for( path_it = config_files.begin(); path_it!=config_files.end(); ++path_it )
  {
    Log::debug() << "Read config file " << path_it->fullName() << "\n";
  }
  Log::debug() << "Configuration read from scripts:\n";
  Log::debug() << util::runtime::indent();
  ResourceMgr::instance().printScript( Log::debug() );
  Log::debug() << util::runtime::dedent();

  Log::debug() << "rundir  : " << PathName(rundir()).fullName() << "\n";
  Log::debug() << util::runtime::dedent();

  // Load factories for static linking
  atlas::grid::load();
}

void atlas_finalize()
{
  Log::debug() << "Atlas finalized\n" << std::flush;
}

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
  static std::string str( Context::instance().runName() );
  return str.c_str();
}

const char* atlas__display_name ()
{
  static std::string str( Context::instance().displayName() );
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

} // namespace atlas


