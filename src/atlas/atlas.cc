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

	struct CallBackContext {
		CallBackContext() { }
		std::string prefix;
		int my_proc;
		int log_proc;
		std::ostream* stream;
	};
	static void callback( void* ctxt, const char* msg )
	{
		CallBackContext& context = *reinterpret_cast<CallBackContext*>(ctxt);
		if( context.log_proc < 0 || context.log_proc == context.my_proc )
			*(context.stream) << context.prefix << msg ;
	}

	CallBackContext debug_ctxt;
	CallBackContext info_ctxt;
	CallBackContext warn_ctxt;
	CallBackContext error_ctxt;
};


Behavior::Behavior() : eckit::LibBehavior()
{
	std::stringstream stream; stream << "["<<MPL::rank()<<"] -- ";
	std::string rank_prefix = ( MPL::size() > 1 ) ? stream.str() : std::string();
	debug_ctxt.stream = &std::cout;
	debug_ctxt.my_proc = MPL::rank();
	debug_ctxt.log_proc = 0;
	debug_ctxt.prefix = rank_prefix + "[DEBUG] -- ";

	info_ctxt.stream = &std::cout;
	info_ctxt.my_proc = MPL::rank();
	info_ctxt.log_proc = 0;
	info_ctxt.prefix = rank_prefix + "[INFO]  -- ";

	warn_ctxt.stream = &std::cout;
	warn_ctxt.my_proc = MPL::rank();
	warn_ctxt.log_proc = 0;
	warn_ctxt.prefix = rank_prefix + "[WARN]  -- ";

	error_ctxt.stream = &std::cout;
	error_ctxt.my_proc = MPL::rank();
	error_ctxt.log_proc = -1;
	error_ctxt.prefix = rank_prefix + "[ERROR] -- ";
}

void Behavior::reconfigure()
{
	CallbackChannel& debug = *dynamic_cast<CallbackChannel*>(&Context::instance().debugChannel());
	debug.register_callback( &callback, &debug_ctxt );

	CallbackChannel& info = *dynamic_cast<CallbackChannel*>(&Context::instance().infoChannel());
	info.register_callback( &callback, &info_ctxt );

	CallbackChannel& warn = *dynamic_cast<CallbackChannel*>(&Context::instance().warnChannel());
	warn.register_callback( &callback, &warn_ctxt );

	CallbackChannel& error = *dynamic_cast<CallbackChannel*>(&Context::instance().errorChannel());
	error.register_callback( &callback, &error_ctxt );
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

void atlas_finalize()
{
  Log::info() << "Atlas finalized" << std::endl;
  MPL::finalize();
}

void atlas__atlas_init(int argc, char **argv)
{
	atlas_init(argc,argv);
}

void atlas__atlas_finalize()
{
	atlas_finalize();
}

} // namespace atlas


