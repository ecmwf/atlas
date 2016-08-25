#ifndef atlas_Behavior_h
#define atlas_Behavior_h

#include "eckit/log/Channel.h"
#include "eckit/runtime/StandardBehavior.h"

#include "atlas/runtime/LogFormat.h"

namespace atlas {
namespace runtime {

//----------------------------------------------------------------------------------------------------------------------

eckit::Channel& standard_out();

eckit::Channel& standard_error();

struct ChannelConfig
{
  std::string logfile_path;
  int  console_rank;
  bool console_enabled;
  bool logfile_enabled;
  bool callback_enabled;
  LogFormat* console_format;
  LogFormat* logfile_format;

  ChannelConfig();
  ~ChannelConfig();

  void apply(eckit::Channel& ch);
};

class Behavior : public eckit::StandardBehavior {

public:

  /// Contructors
  Behavior();

  /// Destructor
  ~Behavior() {}

  /// Info channel
  virtual eckit::Channel& infoChannel();

  /// Warning channel
  virtual eckit::Channel& warnChannel();

  /// Error channel
  virtual eckit::Channel& errorChannel();

  /// Debug channel
  virtual eckit::Channel& debugChannel();

  /// Stats channel
  virtual eckit::Channel& statsChannel();

  enum ChannelCategory { ERROR=0, WARN=1, INFO=2, DEBUG=3, STATS=4 };

  virtual eckit::Channel& channel(int cat);

  /// Configuration
  virtual void reconfigure();

  virtual eckit::FileReadPolicy fileReadPolicy();

private:

  ChannelConfig debug_ctxt;
  ChannelConfig info_ctxt;
  ChannelConfig warn_ctxt;
  ChannelConfig error_ctxt;
  ChannelConfig stats_ctxt;

  int read_policy_;
};

//----------------------------------------------------------------------------------------------------------------------

} // namespace runtime
} // namespace atlas

#endif // atlas_Behavior_h

