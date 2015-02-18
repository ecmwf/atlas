#ifndef atlas_LogFormat_h
#define atlas_LogFormat_h

#include <string>
#include "eckit/log/FormatBuffer.h"
#include "eckit/log/FormatChannel.h"

namespace atlas {

class LogFormat : public eckit::FormatBuffer {
public:

    LogFormat( std::size_t size = 1024 );

    virtual ~LogFormat(){ pubsync(); }

    virtual void beginLine();
    virtual void endLine();

    std::string parsedPrefix() const;

    const std::string& prefix() const;

    void setPrefix( const std::string& );

protected:

    std::string prefix_;

    mutable std::map<std::string,std::string> subst_;

};

class FormattedChannel : public eckit::FormatChannel {
public:

  FormattedChannel( std::ostream* channel, LogFormat* format );

  FormattedChannel( std::ostream& channel, LogFormat* format );

  virtual ~FormattedChannel();

  void set_prefix( const std::string& prefix );

  const std::string& prefix() const;

private:
  std::ostream* channel_;
  LogFormat* format_;
};


} // namespace atlas

#endif // atlas_LogFormat_h
