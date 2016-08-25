#include "atlas_f/internals/atlas_logging.h"
#include "eckit/runtime/Context.h"
#include "eckit/runtime/ContextBehavior.h"
#include "eckit/log/MultiChannel.h"
#include "eckit/log/CallbackChannel.h"
#include "atlas/runtime/LogFormat.h"
#include "atlas/runtime/Behavior.h"
#include "atlas/runtime/Log.h"

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

class FormattedFortranUnitChannel: public runtime::FormattedChannel
{
public:
  FormattedFortranUnitChannel( FortranUnitChannel* fortran_unit, runtime::LogFormat* format ):
    runtime::FormattedChannel( fortran_unit, format )
  {
    fortran_unit_ = fortran_unit;
  }
  virtual ~FormattedFortranUnitChannel()
  {
  }

  int unit() const
  {
    return fortran_unit_->unit();
  }

private:
  FortranUnitChannel* fortran_unit_;
};

} // namespace atlas


void atlas__log_set_debug (int level)
{
  eckit::Context::instance().debug( level );
}

void atlas__log_debug(int lvl, char *msg, int endl, int flush)
{
  atlas::Log::debug() << msg;
  if( endl )
    atlas::Log::debug() << std::endl;
  else if ( flush )
    atlas::Log::debug() << std::flush;
}

void atlas__log_cat(int cat, int lvl, char *msg, int endl, int flush)
{
  atlas::Log::channel(cat) << msg;
  if( endl )
    atlas::Log::channel(cat) << std::endl;
  else if ( flush )
  atlas::Log::channel(cat) << std::flush;
}

// ----------------------------------------------------------------------------

MultiChannel& atlas__get_log_channel(int cat)
{
  Channel& ch = atlas::Log::channel(cat);
  MultiChannel* mc = dynamic_cast< MultiChannel* > (&ch);
  if( !mc )
    throw BadCast("Channel is not a MultiChannel,"
                  "so I cannot connect fortran unit to it. "
                  "Did you forget to call atlas_init()?",Here());
  return *mc;
}

// ----------------------------------------------------------------------------

MultiChannel* atlas__LogChannel_cat(int cat)
{
  return &atlas__get_log_channel(cat);
}

// ----------------------------------------------------------------------------

void atlas__LogChannel__connect_stdout( MultiChannel* ch )
{
  if( !ch->has("console") )
    ch->add( "console" , new FormatChannel(runtime::standard_out(), new runtime::LogFormat() ) );
}

void atlas__logcat__connect_stdout (int cat)
{
  atlas__LogChannel__connect_stdout( atlas__LogChannel_cat(cat) );
}

// ----------------------------------------------------------------------------

void atlas__LogChannel__disconnect_stdout (MultiChannel* ch)
{
  if( ch->has("console") )
  {
    ch->remove("console");
  }
}

void atlas__logcat__disconnect_stdout (int cat)
{
  atlas__LogChannel__disconnect_stdout( atlas__LogChannel_cat(cat) );
}

// ----------------------------------------------------------------------------

void atlas__LogChannel__connect_stderr( MultiChannel* ch )
{
  if( !ch->has("stderr") )
    ch->add( "stderr" , new FormatChannel(runtime::standard_error(), new runtime::LogFormat() ) );
}

void atlas__logcat__connect_stderr (int cat)
{
  atlas__LogChannel__connect_stderr( atlas__LogChannel_cat(cat) );
}

// ----------------------------------------------------------------------------


void atlas__LogChannel__disconnect_stderr (MultiChannel* ch)
{
  if( ch->has("stderr") )
  {
    ch->remove("stderr");
  }
}

void atlas__logcat__disconnect_stderr (int cat)
{
  atlas__LogChannel__disconnect_stderr( atlas__LogChannel_cat(cat) );
}

// ----------------------------------------------------------------------------

void atlas__LogChannel__connect_fortran_unit (MultiChannel* ch, int unit)
{
  FormattedFortranUnitChannel* formatted_fortran_unit =
      new FormattedFortranUnitChannel( new FortranUnitChannel(unit), new runtime::LogFormat() );

  std::stringstream channel_name; channel_name << "fortran_unit_"<<unit;
  ch->add( channel_name.str() , formatted_fortran_unit ); // pass ownership
}

void atlas__logcat__connect_fortran_unit (int cat, int unit)
{
  atlas__LogChannel__connect_fortran_unit( atlas__LogChannel_cat(cat), unit );
}

// ----------------------------------------------------------------------------

void atlas__LogChannel__disconnect_fortran_unit ( MultiChannel* ch,  int unit )
{
  std::stringstream channel_name; channel_name << "fortran_unit_"<<unit;
  if( ch->has(channel_name.str()))
    ch->remove(channel_name.str());
}

void atlas__logcat__disconnect_fortran_unit (int cat, int unit)
{
  atlas__LogChannel__disconnect_fortran_unit( atlas__LogChannel_cat(cat), unit );
}

// ----------------------------------------------------------------------------

void atlas__LogChannel__set_prefix (MultiChannel* ch, char* prefix)
{
  MultiChannel::iterator it;
  for( it=ch->begin(); it!=ch->end(); ++it )
  {
     runtime::FormattedChannel* formatted_ch = dynamic_cast<runtime::FormattedChannel*>(it->second.get());
     if( formatted_ch )
       formatted_ch->format().set_prefix(std::string(prefix));
  }
}

// ----------------------------------------------------------------------------

void atlas__LogChannel__set_prefix_stdout (MultiChannel* ch, char* prefix)
{
  if( ch->has("console") )
  {
    runtime::FormattedChannel* formatted_ch = dynamic_cast<runtime::FormattedChannel*>(&ch->get("console"));
    if( !formatted_ch )
      throw BadCast("Cannot cast channel to atlas::FormattedChannel");
    formatted_ch->format().set_prefix(std::string(prefix));
  }
}

void atlas__logcat__set_prefix_stdout(int cat, char *prefix)
{
  atlas__LogChannel__set_prefix_stdout( atlas__LogChannel_cat(cat), prefix );
}

// ----------------------------------------------------------------------------

void atlas__LogChannel__set_prefix_stderr (MultiChannel* ch, char* prefix)
{
  if( ch->has("stderr") )
  {
    runtime::FormattedChannel* formatted_ch = dynamic_cast<runtime::FormattedChannel*>(&ch->get("stderr"));
    if( !formatted_ch )
      throw BadCast("Cannot cast channel to runtime::FormattedChannel");
    formatted_ch->format().set_prefix(std::string(prefix));
  }
}

void atlas__logcat__set_prefix_stderr(int cat, char *prefix)
{
  atlas__LogChannel__set_prefix_stderr( atlas__LogChannel_cat(cat), prefix );
}

// ----------------------------------------------------------------------------

void atlas__LogChannel__set_prefix_fortran_unit ( MultiChannel* ch,  int unit, char* prefix )
{
  std::stringstream channel_name; channel_name << "fortran_unit_"<<unit;
  if( !ch->has(channel_name.str()) )
  {
    std::stringstream msg;
    msg << "Channel does not have fortran unit " << unit << " connected. I cannot set the prefix.";
    throw BadParameter( msg.str(), Here() );
  }

  runtime::FormattedChannel* formatted_ch = dynamic_cast<runtime::FormattedChannel*>(&ch->get(channel_name.str()));
  if( !formatted_ch )
    throw BadCast("Cannot cast channel to runtime::FormattedChannel");
  formatted_ch->format().set_prefix(prefix);
}

void atlas__logcat__set_prefix_fortran_unit (int cat, int unit, char* prefix )
{
  atlas__LogChannel__set_prefix_fortran_unit( atlas__LogChannel_cat(cat), unit, prefix );
}

// ----------------------------------------------------------------------------

void atlas__LogChannel__indent (MultiChannel* ch, char* indent)
{
  MultiChannel::iterator it;
  for( it=ch->begin(); it!=ch->end(); ++it )
  {
     runtime::FormattedChannel* formatted_ch = dynamic_cast<runtime::FormattedChannel*>(it->second.get());
     if( formatted_ch )
       formatted_ch->format().indent( std::string(indent) );
  }
}

void atlas__logcat__indent (int cat, char* indent)
{
  atlas__LogChannel__indent( atlas__LogChannel_cat(cat), indent );
}


// ----------------------------------------------------------------------------

void atlas__LogChannel__indent_stdout (MultiChannel* ch, char* indent)
{
  if( ch->has("console") )
  {
    runtime::FormattedChannel* formatted_ch = dynamic_cast<runtime::FormattedChannel*>(&ch->get("console"));
    if( !formatted_ch )
      throw BadCast("Cannot cast channel to runtime::FormattedChannel");
    formatted_ch->format().indent( std::string(indent) );
  }
}

// ----------------------------------------------------------------------------

void atlas__LogChannel__indent_stderr (MultiChannel* ch, char* indent)
{
  if( ch->has("stderr") )
  {
    runtime::FormattedChannel* formatted_ch = dynamic_cast<runtime::FormattedChannel*>(&ch->get("stderr"));
    if( !formatted_ch )
      throw BadCast("Cannot cast channel to runtime::FormattedChannel");
    formatted_ch->format().indent( std::string(indent) );
  }
}

// ----------------------------------------------------------------------------

void atlas__LogChannel__indent_fortran_unit (MultiChannel* ch, int unit, char* indent)
{
  std::stringstream channel_name; channel_name << "fortran_unit_"<<unit;
  if( !ch->has(channel_name.str()) )
  {
    std::stringstream msg;
    msg << "Channel does not have fortran unit " << unit << " connected. I cannot indent.";
    throw BadParameter( msg.str(), Here() );
  }

  runtime::FormattedChannel* formatted_ch = dynamic_cast<runtime::FormattedChannel*>(&ch->get(channel_name.str()));
  if( !formatted_ch )
    throw BadCast("Cannot cast channel to runtime::FormattedChannel");
  formatted_ch->format().indent( std::string(indent) );
}

// ----------------------------------------------------------------------------

void atlas__LogChannel__dedent (MultiChannel* ch)
{
  MultiChannel::iterator it;
  for( it=ch->begin(); it!=ch->end(); ++it )
  {
     runtime::FormattedChannel* formatted_ch = dynamic_cast<runtime::FormattedChannel*>(it->second.get());
     if( formatted_ch )
       formatted_ch->format().dedent();
  }
}

void atlas__logcat__dedent (int cat)
{
  atlas__LogChannel__dedent( atlas__LogChannel_cat(cat) );
}

// ----------------------------------------------------------------------------

void atlas__LogChannel__dedent_stdout (MultiChannel* ch)
{
  if( ch->has("console") )
  {
    runtime::FormattedChannel* formatted_ch = dynamic_cast<runtime::FormattedChannel*>(&ch->get("console"));
    if( !formatted_ch )
      throw BadCast("Cannot cast channel to atlas::FormattedChannel");
    formatted_ch->format().dedent();
  }
}

// ----------------------------------------------------------------------------

void atlas__LogChannel__dedent_stderr (MultiChannel* ch)
{
  if( ch->has("stderr") )
  {
    runtime::FormattedChannel* formatted_ch = dynamic_cast<runtime::FormattedChannel*>(&ch->get("stderr"));
    if( !formatted_ch )
      throw BadCast("Cannot cast channel to runtime::FormattedChannel");
    formatted_ch->format().dedent();
  }
}

// ----------------------------------------------------------------------------

void atlas__LogChannel__dedent_fortran_unit (MultiChannel* ch, int unit)
{
  std::stringstream channel_name; channel_name << "fortran_unit_"<<unit;
  if( !ch->has(channel_name.str()) )
  {
    std::stringstream msg;
    msg << "Channel does not have fortran unit " << unit << " connected. I cannot dedent.";
    throw BadParameter( msg.str(), Here() );
  }

  runtime::FormattedChannel* formatted_ch = dynamic_cast<runtime::FormattedChannel*>(&ch->get(channel_name.str()));
  if( !formatted_ch )
    throw BadCast("Cannot cast channel to runtime::FormattedChannel");
  formatted_ch->format().dedent();
}

// ----------------------------------------------------------------------------

void atlas__LogChannel__clear_indentation (MultiChannel* ch)
{
  MultiChannel::iterator it;
  for( it=ch->begin(); it!=ch->end(); ++it )
  {
     runtime::FormattedChannel* formatted_ch = dynamic_cast<runtime::FormattedChannel*>(it->second.get());
     if( formatted_ch )
       formatted_ch->format().clear_indentation();
  }
}

void atlas__logcat__clear_indentation (int cat)
{
  atlas__LogChannel__clear_indentation( atlas__LogChannel_cat(cat) );
}

// ----------------------------------------------------------------------------

void atlas__LogChannel__clear_indentation_stdout (MultiChannel* ch)
{
  if( ch->has("console") )
  {
    runtime::FormattedChannel* formatted_ch = dynamic_cast<runtime::FormattedChannel*>(&ch->get("console"));
    if( !formatted_ch )
      throw BadCast("Cannot cast channel to runtime::FormattedChannel");
    formatted_ch->format().clear_indentation();
  }
}

// ----------------------------------------------------------------------------

void atlas__LogChannel__clear_indentation_stderr (MultiChannel* ch)
{
  if( ch->has("stderr") )
  {
    runtime::FormattedChannel* formatted_ch = dynamic_cast<runtime::FormattedChannel*>(&ch->get("stderr"));
    if( !formatted_ch )
      throw BadCast("Cannot cast channel to runtime::FormattedChannel");
    formatted_ch->format().clear_indentation();
  }
}

// ----------------------------------------------------------------------------

void atlas__LogChannel__clear_indentation_fortran_unit (MultiChannel* ch, int unit)
{
  std::stringstream channel_name; channel_name << "fortran_unit_"<<unit;
  if( !ch->has(channel_name.str()) )
  {
    std::stringstream msg;
    msg << "Channel does not have fortran unit " << unit << " connected. I cannot clear indentation.";
    throw BadParameter( msg.str(), Here() );
  }

  runtime::FormattedChannel* formatted_ch = dynamic_cast<runtime::FormattedChannel*>(&ch->get(channel_name.str()));
  if( !formatted_ch )
    throw BadCast("Cannot cast channel to runtime::FormattedChannel");
  formatted_ch->format().clear_indentation();
}

// ----------------------------------------------------------------------------



