#include <typeinfo>
#include <unistd.h>

#include "eckit/log/ChannelBuffer.h"
#include "eckit/log/CallbackChannel.h"
#include "eckit/log/MultiChannel.h"
#include "eckit/log/FileChannel.h"
#include "eckit/filesystem/LocalPathName.h"
#include "eckit/filesystem/PathName.h"
#include "eckit/thread/ThreadSingleton.h"
#include "eckit/thread/AutoLock.h"
#include "eckit/thread/Mutex.h"
#include "eckit/thread/Once.h"
#include "eckit/thread/ThreadSingleton.h"
#include "eckit/utils/Translator.h"
#include "atlas/runtime/Log.h"
#include "atlas/runtime/Behavior.h"
#include "atlas/runtime/LogFormat.h"
#include "eckit/mpi/Comm.h"

using namespace eckit;

namespace atlas {
namespace runtime {

namespace {
enum { DIRECT_FILE_READ_POLICY=0, PARALLEL_FILE_READ_POLICY=1 };
}

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
  PathName file_path;
  CreateLogFile(const PathName& path) : file_path(path) {}
  FileChannel* operator()()
  {
    char s[6];
    std::sprintf(s, "%05zu",eckit::mpi::comm().rank());
    FileChannel* ch = new FileChannel(LocalPathName(file_path+".p"+std::string(s))) ;
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

namespace { // anonymous
template <typename T>
static bool get_env(const std::string& var,T& val){
  const char* env = ::getenv(var.c_str());
  if( env ) {
    std::string val_str = env;
    val = eckit::Translator<std::string,T>()(val_str);
    return true;
  }
  return false;
}
} // namespace anonymous

ChannelConfig::ChannelConfig()
{
  int logfile_rank = -1;
  get_env("ATLAS_LOGFILE_TASK",logfile_rank);
  
  logfile_path = "";
  get_env("ATLAS_LOGFILE",logfile_path);
  
  logfile_enabled = !logfile_path.empty() && ( logfile_rank < 0 || size_t(logfile_rank) == eckit::mpi::comm().rank() );

  console_rank = 0;
  get_env("ATLAS_CONSOLE_TASK",console_rank);

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

  if( console_enabled && !mc->has("console") && (console_rank < 0 || size_t(console_rank) == eckit::mpi::comm().rank()) )
    mc->add( "console" , new FormattedChannel(standard_out(),console_format) );

  if( mc->has("console") && (!console_enabled || (console_rank >= 0 && size_t(console_rank) != eckit::mpi::comm().rank() ) ) )
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

Behavior::Behavior() /* : ParallelContextBehavior() */
{
  // Read Policy
  const char* read_policy = ::getenv("ATLAS_FILE_READ_POLICY");
  read_policy_ = PARALLEL_FILE_READ_POLICY;
  if( read_policy )
  {
    if(
        std::string(read_policy) == "direct" ||
        std::string(read_policy) == "DIRECT" ||
        std::string(read_policy) == "0" ) {
      read_policy_ = DIRECT_FILE_READ_POLICY;
    }
    else if(
        std::string(read_policy) == "parallel" ||
        std::string(read_policy) == "PARALLEL" ||
        std::string(read_policy) == "1" ) {
      read_policy_ = PARALLEL_FILE_READ_POLICY;
    }
  }

  // Console format
  char p[6];
  std::sprintf(p, "%05zu",eckit::mpi::comm().rank());
  debug_ctxt.console_format->set_prefix("[%p] (%Y-%m-%d T %H:%M:%S) (D) -- ");
  info_ctxt. console_format->set_prefix("[%p] (%Y-%m-%d T %H:%M:%S) (I) -- ");
  warn_ctxt. console_format->set_prefix("[%p] (%Y-%m-%d T %H:%M:%S) (W) -- ");
  error_ctxt.console_format->set_prefix("[%p] (%Y-%m-%d T %H:%M:%S) (E) -- ");
  stats_ctxt.console_format->set_prefix("[%p] (S) -- ");

  // Logfile format
  debug_ctxt.logfile_format->set_prefix("(%Y-%m-%d T %H:%M:%S) (D) -- ");
  info_ctxt. logfile_format->set_prefix("(%Y-%m-%d T %H:%M:%S) (I) -- ");
  warn_ctxt. logfile_format->set_prefix("(%Y-%m-%d T %H:%M:%S) (W) -- ");
  error_ctxt.logfile_format->set_prefix("(%Y-%m-%d T %H:%M:%S) (E) -- ");
  stats_ctxt.logfile_format->set_prefix("(S) -- ");

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
  Log::debug() << "File read policy: " << read_policy_ << std::endl;
}

eckit::FileReadPolicy Behavior::fileReadPolicy()
{
  switch( read_policy_ )
  {
    case PARALLEL_FILE_READ_POLICY:
//      return eckit::mpi::ParallelContextBehavior::fileReadPolicy();
      return eckit::StandardBehavior::fileReadPolicy();
    case DIRECT_FILE_READ_POLICY:
      return eckit::StandardBehavior::fileReadPolicy();
    default:
      throw eckit::BadParameter("Unrecognised read policy", Here());
  }
  return eckit::StandardBehavior::fileReadPolicy();
}

} // namespace runtime
} // namespace atlas


