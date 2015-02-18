#include <typeinfo>  // std::bad_cast
#include <unistd.h>
#include "eckit/log/Log.h"
#include "eckit/log/ChannelBuffer.h"
#include "eckit/log/CallbackChannel.h"
#include "eckit/log/MultiChannel.h"
#include "eckit/log/FileChannel.h"
#include "eckit/config/Resource.h"
#include "eckit/config/ResourceMgr.h"
#include "eckit/filesystem/PathName.h"
#include "eckit/thread/ThreadSingleton.h"
#include "eckit/thread/AutoLock.h"
#include "eckit/thread/Mutex.h"
#include "eckit/thread/Once.h"
#include "eckit/thread/ThreadSingleton.h"

#include "atlas/Behavior.h"
#include "atlas/mpi/mpi.h"
#include "atlas/LogFormat.h"

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
    std::sprintf(s, "%05zu",eckit::mpi::rank());
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

ChannelConfig::ChannelConfig()
{
  int logfile_rank = Resource<int>("atlas.logfile_task;$ATLAS_LOGFILE_TASK;--logfile_task",-1);
  logfile_path    = Resource<std::string>("atlas.logfile;$ATLAS_LOGFILE;--logfile","");
  logfile_enabled = !logfile_path.empty() && ( logfile_rank < 0 || logfile_rank == eckit::mpi::rank() );
  console_rank = Resource<int>("atlas.console_task;$ATLAS_CONSOLE_TASK;--console_task",0);
  console_enabled = true;
  console_format = new LogFormat();
  logfile_format = new LogFormat();
}

ChannelConfig::~ChannelConfig()
{
}

void ChannelConfig::apply(Channel& ch)
{
  MultiChannel* mc;
  mc = dynamic_cast<MultiChannel*>(&ch);
  if( !mc )
    throw BadCast("Cannot cast Channel to MultiChannel",Here());

  if( logfile_enabled && !mc->has("logfile") )
   mc->add( "logfile", new FormattedChannel(logfile(CreateLogFile(logfile_path)),logfile_format) );

  if( console_enabled && !mc->has("console") && (console_rank < 0 || console_rank == eckit::mpi::rank()) )
    mc->add( "console" , new FormattedChannel(standard_out(),console_format) );

  if( mc->has("console") && (!console_enabled || (console_rank >= 0 && console_rank != eckit::mpi::rank() ) ) )
    mc->remove("console");

  if( !mc->has("callback") )
    mc->add( "callback" , new CallbackChannel() );
}

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

Behavior::Behavior() : ContextBehavior()
{
  // Console format
  char p[6];
  std::sprintf(p, "%05zu",eckit::mpi::rank());
  debug_ctxt.console_format->setPrefix("[%p] (%Y-%m-%d T %H:%M:%S) (D) -- ");
  info_ctxt. console_format->setPrefix("[%p] (%Y-%m-%d T %H:%M:%S) (I) -- ");
  warn_ctxt. console_format->setPrefix("[%p] (%Y-%m-%d T %H:%M:%S) (W) -- ");
  error_ctxt.console_format->setPrefix("[%p] (%Y-%m-%d T %H:%M:%S) (E) -- ");
  stats_ctxt.console_format->setPrefix("[%p] (%Y-%m-%d T %H:%M:%S) (S) -- ");

  // Logfile format
  debug_ctxt.logfile_format->setPrefix("(%Y-%m-%d T %H:%M:%S) (D) -- ");
  info_ctxt. logfile_format->setPrefix("(%Y-%m-%d T %H:%M:%S) (I) -- ");
  warn_ctxt. logfile_format->setPrefix("(%Y-%m-%d T %H:%M:%S) (W) -- ");
  error_ctxt.logfile_format->setPrefix("(%Y-%m-%d T %H:%M:%S) (E) -- ");
  stats_ctxt.logfile_format->setPrefix("(%Y-%m-%d T %H:%M:%S) (S) -- ");

  // Debug configuration
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

/// Info channel
Channel& Behavior::infoChannel()
{
  static ThreadSingleton<Channel,CreateInfoChannel> x;
  return x.instance();
}

/// Warning channel
Channel& Behavior::warnChannel()
{
  static ThreadSingleton<Channel,CreateWarnChannel> x;
  return x.instance();
}

/// Error channel
Channel& Behavior::errorChannel()
{
  static ThreadSingleton<Channel,CreateErrorChannel> x;
  return x.instance();
}

  /// Debug channel
Channel& Behavior::debugChannel()
{
  static ThreadSingleton<Channel,CreateDebugChannel> x;
  return x.instance();
}

/// Stats channel
Channel& Behavior::statsChannel()
{
  static ThreadSingleton<Channel,CreateStatsChannel> x;
  return x.instance();
}

Channel& Behavior::channel(int cat)
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
void Behavior::reconfigure()
{
}

} // namespace atlas


