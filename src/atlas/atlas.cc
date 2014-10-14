#include <typeinfo>  // std::bad_cast

#include <eckit/runtime/LibBehavior.h>
#include <eckit/runtime/Context.h>
#include <eckit/log/Log.h>
#include <eckit/log/ChannelBuffer.h>
#include <eckit/log/Colour.h>
#include <eckit/log/Channel.h>
#include <eckit/log/ChannelBuffer.h>
#include <eckit/log/MultiChannel.h>
#include <eckit/log/FileChannel.h>
#include <eckit/log/FormatChannel.h>
#include "eckit/log/ColorizeFormat.h"
#include <eckit/config/Resource.h>
#include <eckit/config/ResourceMgr.h>
#include <eckit/filesystem/PathName.h>
#include <eckit/thread/ThreadSingleton.h>
#include <eckit/thread/AutoLock.h>
#include <eckit/thread/Mutex.h>
#include <eckit/thread/Once.h>
#include <eckit/thread/ThreadSingleton.h>
#include <eckit/os/BackTrace.h>


#include "atlas/atlas.h"
#include "atlas/mpl/MPL.h"

using namespace eckit;

namespace atlas {


static Once<Mutex> local_mutex;

template< typename TYPE >
struct OutAlloc {
    OutAlloc() {}
    TYPE* operator() ()
    {
    	return new TYPE( new ChannelBuffer( std::cout ) );
    }
};

template< typename TYPE >
struct ErrAlloc {
    ErrAlloc() {}
    TYPE* operator() ()
    {
    	return new TYPE( new ChannelBuffer( std::cerr ) );
    }
};


Channel& standard_out()
{
  static ThreadSingleton<Channel,OutAlloc<Channel> > x;
  return x.instance();
}

Channel& standard_error()
{
  static ThreadSingleton<Channel,ErrAlloc<Channel> > x;
  return x.instance();
}

struct CreateLogFile
{
	LocalPathName file_path;
	CreateLogFile(const LocalPathName& path) : file_path(path) {}
	FileChannel* operator()()
	{
		char s[6];
		std::sprintf(s, "%05d",MPL::rank());
		FileChannel* ch = new FileChannel(file_path+".p"+std::string(s)) ;
		return ch;
	}
};


// nawd: What I would really like is that depending on different file_path,
// a new singleton is created, but in case of same file_path, the logfiles
// should match for different channels
Channel& logfile( const CreateLogFile& alloc)
{
	static ThreadSingleton<Channel,CreateLogFile> x( alloc );
	return x.instance();
}

struct ChannelConfig
{
	std::string logfile_path;
	int  console_rank;
	bool console_enabled;
	bool logfile_enabled;
	bool callback_enabled;
	ColorizeFormat* console_format;
	ColorizeFormat* logfile_format;


	ChannelConfig()
	{
		int logfile_rank = Resource<int>("logfile_task;$ATLAS_LOGFILE_TASK;-logfile_task",-1);
		logfile_path    = Resource<std::string>("logfile;$ATLAS_LOGFILE;-logfile","");
		logfile_enabled = !logfile_path.empty() && ( logfile_rank < 0 || logfile_rank == MPL::rank() );
		console_rank = 0;
		console_enabled = true;
		console_format = new ColorizeFormat();
		logfile_format = new ColorizeFormat();
	}

	~ChannelConfig()
	{
		delete console_format;
		delete logfile_format;
	}

	void apply(Channel& ch)
	{
		MultiChannel* mc;
		try
		{
			mc = dynamic_cast<MultiChannel*>(&ch);
		}
		catch (std::bad_cast& e)
		{
			throw BadCast("Cannot cast Channel to MultiChannel",Here());
		}

		if( logfile_enabled && !mc->has("logfile") )
			mc->add( "logfile", new FormatChannel(logfile(CreateLogFile(logfile_path)),logfile_format) );

		if( console_enabled && !mc->has("console") && (console_rank < 0 || console_rank == MPL::rank()) )
			mc->add( "console" , new FormatChannel(standard_out(),console_format) );

		if( mc->has("console") && (!console_enabled || (console_rank >= 0 && console_rank != MPL::rank() ) ) )
			mc->remove("console");

		if( !mc->has("callback") )
			mc->add( "callback" , new CallbackChannel() );
	}
};


struct CreateChannel
{
	Channel* operator()()
	{
		MultiChannel* mc = new MultiChannel();
		return mc;
	}
};


struct CreateDebugChannel : CreateChannel {};
struct CreateInfoChannel  : CreateChannel {};
struct CreateWarnChannel  : CreateChannel {};
struct CreateErrorChannel : CreateChannel {};
struct CreateStatsChannel : CreateChannel {};


class Behavior : public eckit::ContextBehavior {

public:

	/// Contructors
	Behavior() : ContextBehavior()
	{
		// Console format
		char p[6];
		std::sprintf(p, "%05d",MPL::rank());
		debug_ctxt.console_format->prefix("(P"+std::string(p)+" D) -- ");
		info_ctxt. console_format->prefix("(P"+std::string(p)+" I) -- ");
		warn_ctxt. console_format->prefix("(P"+std::string(p)+" W) -- ");
		error_ctxt.console_format->prefix("(P"+std::string(p)+" E) -- ");
		stats_ctxt.console_format->prefix("(P"+std::string(p)+" S) -- ");

		// Logfile format
		debug_ctxt.logfile_format->prefix("(D) -- ");
		info_ctxt. logfile_format->prefix("(I) -- ");
		warn_ctxt. logfile_format->prefix("(W) -- ");
		error_ctxt.logfile_format->prefix("(E) -- ");
		stats_ctxt.logfile_format->prefix("(S) -- ");

		// Debug configuration
		debug_ctxt.console_enabled = false;
		debug_ctxt.apply(debugChannel());

		// Info configuration
		info_ctxt.apply(infoChannel());

		// Warning configuration
		warn_ctxt.apply(warnChannel());

		// Error configuration
		//error_ctxt.console_rank = -1; // all ranks log errors to console
		error_ctxt.apply(errorChannel());

		// Stats configuration
		stats_ctxt.console_enabled = false;
		stats_ctxt.apply(statsChannel());
	}

	/// Destructor
	~Behavior() {}

	/// Info channel
	virtual Channel& infoChannel()
	{
		static ThreadSingleton<Channel,CreateInfoChannel> x;
		return x.instance();
	}

	/// Warning channel
	virtual Channel& warnChannel()
	{
		static ThreadSingleton<Channel,CreateWarnChannel> x;
		return x.instance();
	}

	/// Error channel
	virtual Channel& errorChannel()
	{
		static ThreadSingleton<Channel,CreateErrorChannel> x;
		return x.instance();
	}

	/// Debug channel
	virtual Channel& debugChannel()
	{
		static ThreadSingleton<Channel,CreateDebugChannel> x;
		return x.instance();
	}

	/// Stats channel
	virtual Channel& statsChannel()
	{
		static ThreadSingleton<Channel,CreateStatsChannel> x;
		return x.instance();
	}

  enum ChannelCategory { ERROR=0, WARN=1, INFO=2, DEBUG=3, STATS=4 };

	virtual Channel& channel(int cat)
	{
		switch( cat ) {
			case ERROR: return errorChannel();
			case WARN:  return warnChannel();
			case INFO:  return infoChannel();
			case DEBUG: return debugChannel();
			case STATS: return statsChannel();
		}
		throw Exception("Logging category "+Translator<int,std::string>()(cat)+" not known.",Here());
		return infoChannel();
	}

	/// Configuration
	virtual void reconfigure()
	{
	}

private:

	ChannelConfig debug_ctxt;
	ChannelConfig info_ctxt;
	ChannelConfig warn_ctxt;
	ChannelConfig error_ctxt;
	ChannelConfig stats_ctxt;
};

void read_config(const LocalPathName& path, const int master_task = 0)
{
	/// This function lets only one MPI task read the config file,
	/// and broadcasts to remaining tasks. This is beneficial at
	/// large scale.
	std::stringstream stream;
	char* buf;
	int buf_len(0);
	if( MPL::rank() == master_task )
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
			MPI_Bcast(&buf_len,1,MPL::TYPE<int >(),master_task,MPI_COMM_WORLD);
			if (buf_len)
				MPI_Bcast(buf,buf_len,MPL::TYPE<char>(),master_task,MPI_COMM_WORLD);
		}
	}
	else
	{
		MPI_Bcast(&buf_len,1,MPL::TYPE<int>(),master_task,MPI_COMM_WORLD);
		if( buf_len )
		{
			buf = new char[buf_len];
			MPI_Bcast(buf,buf_len,MPL::TYPE<char>(),master_task,MPI_COMM_WORLD);
			stream.write(buf,buf_len);
			delete[] buf;
		}
	}
	if (buf_len)
		ResourceMgr::instance().appendConfig(stream);
}


void atlas_init(int argc, char** argv)
{
	Context::instance().setup(argc, argv);
	MPL::init( Context::instance().argc(), Context::instance().argvs() );

	if( Context::instance().argc() > 0 )
	{
		Context::instance().runName(
			PathName(Context::instance().argv(0)).baseName(false));
	}
	Context::instance().displayName( Resource<std::string>("-name","") );

	read_config( LocalPathName( Resource<std::string>("$ATLAS_CONFIGFILE;-atlas_conf","atlas.cfg") ) );

	Context::instance().behavior( new atlas::Behavior() );
	Context::instance().behavior().debug( Resource<int>("debug;$DEBUG;-debug",0) );

	Log::debug() << "Atlas program " << Context::instance().runName() << " initialized" << std::endl;
	Log::debug() << "    Atlas version [" << atlas_version() << "]" << std::endl;
	Log::debug() << "    Atlas git     [" << atlas_git_sha1()<< "]" << std::endl;

#if 0 // DEBUG
	bool  triangulate = Resource< bool > ( &Context::instance(), "meshgen.triangulate;-triangulate", false );
	if( triangulate == false )
	{
		Log::error() << "triangulate =" << triangulate << std::endl;
	}

		Log::error() << "angle =" << Resource< double > ( &Context::instance(), "meshgen.angle", 15.0 ) << std::endl;
#endif

}

void atlas_finalize()
{
  Log::debug() << "Atlas finalized" << std::endl;
  MPL::finalize();
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

const char* atlas__atlas_version()
{
	return atlas_version();
}

const char* atlas__atlas_git_sha1()
{
	return atlas_git_sha1();
}

} // namespace atlas


