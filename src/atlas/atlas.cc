#include <eckit/runtime/LibBehavior.h>
#include <eckit/runtime/Context.h>
#include <eckit/log/Log.h>

#include <eckit/log/ChannelBuffer.h>
#include <eckit/thread/ThreadSingleton.h>

#include <eckit/log/Colour.h>
#include <eckit/log/Channel.h>
#include <eckit/log/ChannelBuffer.h>
#include <eckit/log/FormatChannel.h>
#include "eckit/log/ColorizeFormat.h"

#include <eckit/config/Resource.h>
#include <eckit/thread/ThreadSingleton.h>
#include <eckit/thread/AutoLock.h>
#include <eckit/thread/Mutex.h>
#include <eckit/thread/Once.h>

#include <eckit/os/BackTrace.h>


#include "atlas/atlas.h"
#include "atlas/mpl/MPL.h"

using namespace eckit;

namespace atlas {

class Behavior : public eckit::LibBehavior {

public:

	struct LogContext {
		int myproc;
		int logger_proc;
		std::string rank_prefix;
	};

public:

	Behavior();

	virtual ~Behavior();

	virtual void reconfigure();

private:
	LogContext debug_ctxt;
	LogContext info_ctxt;
	LogContext warn_ctxt;
	LogContext error_ctxt;
};

static void debug_callback( void* ctxt, const char* msg )
{
	Behavior::LogContext& log_ctxt = *reinterpret_cast<Behavior::LogContext*>(ctxt);
	if( log_ctxt.logger_proc < 0 || log_ctxt.logger_proc == log_ctxt.myproc )
		std::cout << log_ctxt.rank_prefix << "[DEBUG] -- " << msg ;
}

static void info_callback( void* ctxt, const char* msg )
{
	Behavior::LogContext& log_ctxt = *reinterpret_cast<Behavior::LogContext*>(ctxt);
	if( log_ctxt.logger_proc < 0 || log_ctxt.logger_proc == log_ctxt.myproc )
		std::cout << log_ctxt.rank_prefix << "[INFO]  -- " << msg ;
}

static void warn_callback( void* ctxt, const char* msg )
{
	Behavior::LogContext& log_ctxt = *reinterpret_cast<Behavior::LogContext*>(ctxt);
	if( log_ctxt.logger_proc < 0 || log_ctxt.logger_proc == log_ctxt.myproc )
		std::cout << log_ctxt.rank_prefix << "[WARN]  -- " << msg ;
}

static void error_callback( void* ctxt, const char* msg )
{
	Behavior::LogContext& log_ctxt = *reinterpret_cast<Behavior::LogContext*>(ctxt);
	if( log_ctxt.logger_proc < 0 || log_ctxt.logger_proc == log_ctxt.myproc )
		std::cout << log_ctxt.rank_prefix << "[ERROR] -- " << msg ;
}


Behavior::Behavior() : eckit::LibBehavior()
{
	std::stringstream stream; stream << "["<<MPL::rank()<<"] -- ";
	std::string rank_prefix = ( MPL::size() > 1 ) ? stream.str() : std::string();
	debug_ctxt.myproc = MPL::rank();
	debug_ctxt.logger_proc = 0;
	debug_ctxt.rank_prefix = rank_prefix;

	info_ctxt.myproc = MPL::rank();
	info_ctxt.logger_proc = 0;
	info_ctxt.rank_prefix = rank_prefix;

	warn_ctxt.myproc = MPL::rank();
	warn_ctxt.logger_proc = 0;
	warn_ctxt.rank_prefix = rank_prefix;

	error_ctxt.myproc = MPL::rank();
	error_ctxt.logger_proc = -1;
	error_ctxt.rank_prefix = rank_prefix;
}

void Behavior::reconfigure()
{
	CallbackChannel& debug = *dynamic_cast<CallbackChannel*>(&Context::instance().debugChannel());
	debug.register_callback( &debug_callback, &debug_ctxt );

	CallbackChannel& info = *dynamic_cast<CallbackChannel*>(&Context::instance().infoChannel());
	info.register_callback( &info_callback, &info_ctxt );

	CallbackChannel& warn = *dynamic_cast<CallbackChannel*>(&Context::instance().warnChannel());
	warn.register_callback( &warn_callback, &warn_ctxt );

	CallbackChannel& error = *dynamic_cast<CallbackChannel*>(&Context::instance().errorChannel());
	error.register_callback( &error_callback, &error_ctxt );
}

Behavior::~Behavior() {}


void atlas_init(int argc, char** argv)
{
	MPL::init(argc,argv);
	Context::instance().setup(argc, argv);

	Context::instance().behavior( new atlas::Behavior() );
	Context::instance().debug( Resource<int>(&Context::instance(),"debug;$DEBUG;-debug",0) );
	Context::instance().reconfigure();
	Log::info() << "Atlas initialized" << std::endl;
	Log::info() << "    version [" << atlas_version() << "]" << std::endl;
	Log::info() << "    git     [" << atlas_git_sha1()<< "]" << std::endl;
}

void atlas_initf(int argc, char** argv)
{
	atlas_init(argc,argv);
}

} // namespace atlas


