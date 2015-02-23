#include <typeinfo>  // std::bad_cast
#include <unistd.h>
#include "eckit/runtime/Context.h"
#include "eckit/log/Log.h"
#include "eckit/log/ChannelBuffer.h"
#include "eckit/log/Colour.h"
#include "eckit/log/Channel.h"
#include "eckit/log/ChannelBuffer.h"
#include "eckit/log/MultiChannel.h"
#include "eckit/log/FileChannel.h"
#include "eckit/log/FormatChannel.h"
#include "eckit/log/TimeStamp.h"
#include "eckit/parser/StringTools.h"
#include "eckit/config/Resource.h"
#include "eckit/config/ResourceMgr.h"
#include "eckit/filesystem/PathName.h"
#include "eckit/thread/ThreadSingleton.h"
#include "eckit/os/BackTrace.h"
#include "eckit/parser/StringTools.h"

#include "atlas/Behavior.h"
#include "atlas/atlas.h"
#include "atlas/grids/grids.h"
#include "atlas/mpi/mpi.h"

using namespace eckit;

namespace atlas {

std::string read_config(const LocalPathName& path, const int master_task = 0)
{
	std::string config;

	/// This function lets only one MPI task read the config file,
	/// and broadcasts to remaining tasks. This is beneficial at
	/// large scale.
	std::stringstream stream;
	char* buf;
	int buf_len(0);
	if( eckit::mpi::rank() == master_task )
	{
		if( path.exists() )
		{
			std::fstream file( path.c_str(), std::ios_base::in );
			if ( file )
			{
				stream << file.rdbuf();
				file.close();
			}
			std::string str = stream.str();
			buf = const_cast<char*>(str.c_str());
			buf_len = str.size();
			MPI_Bcast(&buf_len,1,eckit::mpi::datatype<int >(),master_task,eckit::mpi::comm());
			if (buf_len)
				MPI_Bcast(buf,buf_len,eckit::mpi::datatype<char>(),master_task,eckit::mpi::comm());
		}
		else
		{
			MPI_Bcast(&buf_len,1,eckit::mpi::datatype<int >(),master_task,eckit::mpi::comm());
		}
	}
	else
	{
		MPI_Bcast(&buf_len,1,eckit::mpi::datatype<int>(),master_task,eckit::mpi::comm());
		if( buf_len )
		{
			buf = new char[buf_len];
			MPI_Bcast(buf,buf_len,eckit::mpi::datatype<char>(),master_task,eckit::mpi::comm());
			stream.write(buf,buf_len);
			delete[] buf;
		}
	}
	if (buf_len)
	{
		ResourceMgr::instance().appendConfig(stream);
		config = stream.str();
	}
	return config;
}

std::string rundir()
{
	static LocalPathName cwd( LocalPathName::cwd() );
	return cwd;
}

class Environment
{
public:
  Environment(): finalize_mpi_(true) {}
  static Environment& instance()
  {
    static Environment atlas_state;
    return atlas_state;
  }
  void set_finalize_mpi( bool finalize ) { finalize_mpi_ = finalize; }
  bool finalize_mpi() const { return finalize_mpi_; }
private:
  bool finalize_mpi_;
};

void atlas_init(int argc, char** argv)
{
  if( argc > 0 )
    Context::instance().setup(argc, argv);

  if( eckit::mpi::initialized() ) Environment::instance().set_finalize_mpi(false);

  eckit::mpi::init( Context::instance().argc(), Context::instance().argvs() );

  if( Context::instance().argc() > 0 )
  {
    Context::instance().runName(
      PathName(Context::instance().argv(0)).baseName(false));
  }
  Context::instance().displayName( Resource<std::string>("--name;$ATLAS_DISPLAYNAME","") );

  LocalPathName atlas_config_path ( Resource<std::string>("atlas.configfile;$ATLAS_CONFIGFILE;--atlas_conf","atlas.cfg") );
  std::string atlas_config = read_config( atlas_config_path );

  LocalPathName runname_config_path (
    Resource<std::string>("$"+StringTools::upper(Context::instance().runName())+"_CONFIGFILE",
                          Context::instance().runName()+".cfg") );
  std::string runname_config = read_config( runname_config_path );

  std::string default_conf("");
  if( Context::instance().displayName() != Context::instance().runName() )
    default_conf = Context::instance().displayName()+".cfg";

  LocalPathName displayname_config_path(
    Resource<std::string>(
       "$"+StringTools::upper(Context::instance().displayName())+"_CONFIGFILE;--conf",
       default_conf) );
  std::string displayname_config = read_config( displayname_config_path );

  LocalPathName workdir ( Resource<std::string>("workdir;$ATLAS_WORKDIR;--workdir",rundir()) );
  if( workdir != rundir() )
  {
    if( eckit::mpi::rank() == 0 ) workdir.mkdir();
    eckit::mpi::barrier();
    Log::debug() << "Changing working directory to " << workdir << std::endl;
    chdir(workdir.c_str());
  }

  Context::instance().behavior( new atlas::Behavior() );
  Context::instance().behavior().debug( Resource<int>("debug;$DEBUG;--debug",0) );

  atlas::grids::load();

  Log::debug() << "Atlas program " << Context::instance().displayName() << " initialized" << std::endl;
  Log::debug() << "    Atlas version [" << atlas_version() << "]" << std::endl;
  Log::debug() << "    Atlas git     [" << atlas_git_sha1()<< "]" << std::endl;

  if( atlas_config.size() ) {
    Log::debug() << "    Read config file " << atlas_config_path.fullName() << std::endl;
	Log::debug() << atlas_config << std::endl;
  }

  if( runname_config.size() ) {
    Log::debug() << "    Read config file " << runname_config_path.fullName() << std::endl;
    Log::debug() << runname_config << std::endl;
  }
  if( displayname_config.size() ) {
  	Log::debug() << "    Read config file " << displayname_config_path.fullName() << std::endl;
    Log::debug() << displayname_config << std::endl;
  }

  Log::debug() << "    rundir  : " << LocalPathName(rundir()).fullName() << std::endl;
  Log::debug() << "    workdir : " << LocalPathName::cwd().fullName() << std::endl;
}

void atlas_finalize()
{
  Log::debug() << "Atlas finalized\n" << std::flush;
  // Log::debug() << "Some more\n";
  //Log::debug() << std::flush;
  // << std::endl;
  //if( Environment::instance().finalize_mpi() )
  //  eckit::mpi::finalize();
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


