#ifndef atlas_LogFormat_h
#define atlas_LogFormat_h

#include <string>
#include <stack>
#include "eckit/log/FormatBuffer.h"
#include "eckit/log/FormatChannel.h"

namespace atlas {
namespace runtime {

class LogFormat : public eckit::FormatBuffer {
public:

    LogFormat( std::size_t size = 1024 );

    virtual ~LogFormat(){ pubsync(); }

    const std::string& prefix() const;

    void set_prefix( const std::string& );

    void indent( const std::string& = std::string("  ") );

    void dedent();

    void clear_indentation();

private:

    std::string parsed_prefix() const;
    virtual void beginLine();
    virtual void endLine();

private:

    std::vector<std::string> indent_stack_;
    std::string indent_;
    std::string prefix_;

    mutable std::map<std::string,std::string> subst_;

};

class FormattedChannel : public eckit::FormatChannel {
public:

  FormattedChannel( std::ostream* channel, LogFormat* format );

  FormattedChannel( std::ostream& channel, LogFormat* format );

  virtual ~FormattedChannel();

  const LogFormat& format() const { return *format_; };
        LogFormat& format()       { return *format_; };

private:
  // std::ostream* channel_;
  LogFormat* format_;
};

class indent
{
public:
  indent(const std::string& s = std::string("  ")) :
    indent_(s)
  {}
  operator std::string() const { return indent_; }
private:
  std::string indent_;
};

class dedent
{
public:
  dedent() {}
};

std::ostream& operator<< ( std::ostream& stream, const indent& );
std::ostream& operator<< ( std::ostream& stream, const dedent& );

} // namespace runtime
} // namespace atlas

#endif // atlas_LogFormat_h
