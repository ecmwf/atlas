#include "atlas/fortran/atlas_logging.h"
#include "eckit/log/Log.h"
#include "eckit/runtime/Context.h"
#include "eckit/runtime/ContextBehavior.h"
#include "eckit/log/MultiChannel.h"
#include "eckit/log/FormatChannel.h"
#include "eckit/log/CallbackChannel.h"
#include "atlas/LogFormat.h"
#include "atlas/Behavior.h"

using namespace eckit;
using namespace atlas;

namespace atlas {

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

  int unit() const { return unit_; }
private:
  int unit_;
};

class FormattedFortranUnitChannel: public FormatChannel
{
public:
  FormattedFortranUnitChannel( FortranUnitChannel* fortran_unit, LogFormat* format ):
    FormatChannel( fortran_unit, format )
  {
    fortran_unit_ = fortran_unit;
    format_ = format;
  }
  virtual ~FormattedFortranUnitChannel()
  {
  }

  void set_prefix( const std::string& prefix )
  {
    format_->setPrefix(prefix);
  }

  int unit() const
  {
    return fortran_unit_->unit();
  }

private:
  FortranUnitChannel* fortran_unit_;
  LogFormat* format_;
};

} // namespace atlas


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

// ----------------------------------------------------------------------------

MultiChannel& atlas__get_log_channel(int cat)
{
  Channel& ch = eckit::Log::channel(cat);
  MultiChannel* mc = dynamic_cast< MultiChannel* > (&ch);
  if( !mc )
    throw BadCast("Channel is not a MultiChannel,"
                  "so I cannot connect fortran unit to it",Here());
  return *mc;
}

// ----------------------------------------------------------------------------

MultiChannel* atlas__log_channel(int cat)
{
  return &atlas__get_log_channel(cat);
}

// ----------------------------------------------------------------------------

void atlas__Channel__connect_stdout( MultiChannel* ch )
{
  if( !ch->has("console") )
    ch->add( "console" , new FormatChannel(standard_out(), new LogFormat() ) );
}

void atlas__cat__connect_stdout (int cat)
{
  atlas__Channel__connect_stdout( atlas__log_channel(cat) );
}

// ----------------------------------------------------------------------------

void atlas__Channel__disconnect_stdout (MultiChannel* ch)
{
  if( ch->has("console") )
  {
    ch->remove("console");
  }
}

void atlas__cat__disconnect_stdout (int cat)
{
  atlas__Channel__disconnect_stdout( atlas__log_channel(cat) );
}

// ----------------------------------------------------------------------------

void atlas__Channel__connect_stderr( MultiChannel* ch )
{
  if( !ch->has("stderr") )
    ch->add( "stderr" , new FormatChannel(standard_error(), new LogFormat() ) );
}

void atlas__cat__connect_stderr (int cat)
{
  atlas__Channel__connect_stderr( atlas__log_channel(cat) );
}

// ----------------------------------------------------------------------------


void atlas__Channel__disconnect_stderr (MultiChannel* ch)
{
  if( ch->has("stderr") )
  {
    ch->remove("stderr");
  }
}

void atlas__cat__disconnect_stderr (int cat)
{
  atlas__Channel__disconnect_stderr( atlas__log_channel(cat) );
}

// ----------------------------------------------------------------------------

void atlas__Channel__connect_fortran_unit (MultiChannel* ch, int unit)
{
  FormattedFortranUnitChannel* formatted_fortran_unit =
      new FormattedFortranUnitChannel( new FortranUnitChannel(unit), new LogFormat() );

  std::stringstream channel_name; channel_name << "fortran_unit_"<<unit;
  ch->add( channel_name.str() , formatted_fortran_unit ); // pass ownership
}

void atlas__cat__connect_fortran_unit (int cat, int unit)
{
  atlas__Channel__connect_fortran_unit( atlas__log_channel(cat), unit );
}

// ----------------------------------------------------------------------------

void atlas__Channel__disconnect_fortran_unit ( MultiChannel* ch,  int unit )
{
  std::stringstream channel_name; channel_name << "fortran_unit_"<<unit;
  if( ch->has(channel_name.str()))
    ch->remove(channel_name.str());
}

void atlas__cat__disconnect_fortran_unit (int cat, int unit)
{
  atlas__Channel__disconnect_fortran_unit( atlas__log_channel(cat), unit );
}

// ----------------------------------------------------------------------------

void atlas__Channel__set_prefix_stdout (AtlasChannel* ch, char* prefix)
{
  if( ch->has("console") )
  {
    FormattedChannel* formatted_ch = dynamic_cast<atlas::FormattedChannel*>(&ch->get("console"));
    if( !formatted_ch )
      throw BadCast("Cannot cast channel to atlas::FormattedChannel");
    formatted_ch->set_prefix(std::string(prefix));
  }
}


void atlas__cat__set_prefix_stdout(int cat, char *prefix)
{
  atlas__Channel__set_prefix_stdout( atlas__log_channel(cat), prefix );
}

// ----------------------------------------------------------------------------

void atlas__Channel__set_prefix_stderr (AtlasChannel* ch, char* prefix)
{
  if( ch->has("stderr") )
  {
    FormattedChannel* formatted_ch = dynamic_cast<atlas::FormattedChannel*>(&ch->get("stderr"));
    if( !formatted_ch )
      throw BadCast("Cannot cast channel to atlas::FormattedChannel");
    formatted_ch->set_prefix(std::string(prefix));
  }
}

void atlas__cat__set_prefix_stderr(int cat, char *prefix)
{
  atlas__Channel__set_prefix_stderr( atlas__log_channel(cat), prefix );
}

// ----------------------------------------------------------------------------


void atlas__Channel__set_prefix_fortran_unit ( MultiChannel* ch,  int unit, char* prefix )
{
  std::stringstream channel_name; channel_name << "fortran_unit_"<<unit;
  if( !ch->has(channel_name.str()) )
  {
    std::stringstream msg;
    msg << "Channel does not have fortran unit " << unit << " connected. I cannot change the prefix.";
    throw BadParameter( msg.str(), Here() );
  }

  FormattedFortranUnitChannel* formatted_fortran_unit =
      dynamic_cast< FormattedFortranUnitChannel* >( &ch->get(channel_name.str()) );

  if( !formatted_fortran_unit )
    throw BadCast("Cannot cast channel to a FormattedFortranUnitChannel");

  formatted_fortran_unit->set_prefix(prefix);
}

void atlas__cat__set_prefix_fortran_unit (int cat, int unit, char* prefix )
{
  atlas__Channel__set_prefix_fortran_unit( atlas__log_channel(cat), unit, prefix );
}

// ----------------------------------------------------------------------------



