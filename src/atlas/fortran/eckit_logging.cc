#include "atlas/fortran/eckit_logging.h"
#include "eckit/log/Log.h"
#include "eckit/config/Resource.h"
#include "eckit/runtime/Context.h"

extern "C"
{

	void eckit__log_debug_set_level (int level)
	{
		eckit::Log::info() << "Set Log::debug() to " << level << std::endl;
		eckit::Context::instance().debug( level );
		eckit::Context::instance().reconfigure();
	}

	void eckit__log_debug (int level, char* name)
	{
		eckit::Log::debug(level) << name;
	}

	void eckit__log_info (char* name)
	{
		eckit::Log::info() << name << std::flush;
	}

	void eckit__log_warning (char* name)
	{
		eckit::Log::warning() << name;
	}

	void eckit__log_error (char* name)
	{
		eckit::Log::error() << name;
	}

	void eckit__log_panic (char* name)
	{
		eckit::Log::panic() << name;
	}

	void eckit__log_debug_endl (int level, char* name)
	{
		eckit::Log::debug(level) << name << std::endl;
	}

	void eckit__log_info_endl (char* name)
	{
		eckit::Log::info() << name <<std::endl;
	}

	void eckit__log_warning_endl (char* name)
	{
		eckit::Log::warning() << name << std::endl;
	}

	void eckit__log_error_endl (char* name)
	{
		eckit::Log::error() << name << std::endl;
	}

	void eckit__log_panic_endl (char* name)
	{
		eckit::Log::panic() << name << std::endl;
	}

}
