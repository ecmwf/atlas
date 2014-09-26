#include "atlas/fortran/atlas_logging.h"
#include "eckit/log/Log.h"
#include "eckit/runtime/Context.h"

extern "C"
{

	void atlas__log_debug_set_level (int level)
	{
		eckit::Log::info() << "Set Log::debug() to " << level << std::endl;
		eckit::Context::instance().debug( level );
		eckit::Context::instance().reconfigure();
	}

	void atlas__log_debug(int lvl, char *msg, int endl, int flush)
	{
		eckit::Log::debug(lvl) << msg;
		if( endl )
			eckit::Log::debug(lvl) << std::endl;
		else if ( flush )
			eckit::Log::debug(lvl) << std::flush;
	}

	void atlas__log(int cat, int lvl, char *msg, int endl, int flush)
	{
		eckit::Log::channel(cat,lvl) << msg;
		if( endl )
			eckit::Log::channel(cat,lvl) << std::endl;
		else if ( flush )
			eckit::Log::channel(cat,lvl) << std::flush;
	}

}
