#include "atlas/fortran/atlas_logging.h"
#include "eckit/log/Log.h"
#include "eckit/runtime/Context.h"
#include "eckit/runtime/ContextBehavior.h"
#include "eckit/log/MultiChannel.h"
#include "eckit/log/CallbackChannel.h"

using namespace eckit;

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

extern "C" { void atlas_write_to_fortran_unit(int unit, const char* msg); }

void write_to_fortran_unit( void* ctxt, const char* msg )
{
  atlas_write_to_fortran_unit( *static_cast<int*>(ctxt), msg );
}

class FortranUnitChannel: public CallbackChannel {
public:

  FortranUnitChannel(int unit) : CallbackChannel()
  {
    unit_ = unit;
    register_callback(&write_to_fortran_unit, &unit_);
  }

private:
  int unit_;
};


void atlas__cat__connect_fortran_unit (int cat, int unit)
{
  Channel& ch = eckit::Log::channel(cat);
  MultiChannel* mc = dynamic_cast< MultiChannel* > (&ch);
  if( !mc )
    throw BadCast("Channel is not a MultiChannel,"
                  "so I cannot connect fortran unit to it",Here());

  FortranUnitChannel* fortran_unit = new FortranUnitChannel(unit);

  std::stringstream channel_name; channel_name << "fortran_unit_"<<unit;
  mc->add( channel_name.str() , fortran_unit ); // pass ownership
}

