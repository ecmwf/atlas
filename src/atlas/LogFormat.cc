#include <cmath>

#include "eckit/log/TimeStamp.h"
#include "atlas/LogFormat.h"
#include "atlas/mpi/mpi.h"
using namespace eckit;

namespace atlas {

namespace {

bool replace(std::string& str, const std::string& from, const std::string& to)
{
  size_t start_pos = str.find(from);
  if(start_pos == std::string::npos)
    return false;
  str.replace(start_pos, from.length(), to);
  return true;
}

}

LogFormat::LogFormat(std::size_t size)
 : FormatBuffer(0,size)
{
  for( int j=1; j<10; ++j )
  {
    char p[j+1];
    char P[j+1];
    std::stringstream fmt; fmt << "%0"<<j<<"zu";
    std::sprintf(p, fmt.str().c_str(),eckit::mpi::rank());
    std::sprintf(P, fmt.str().c_str(),eckit::mpi::rank()+1);
    std::stringstream sub; sub << "%"<<j<<"p";
    subst_[sub.str()] = std::string(p);
    std::stringstream sub_plus_1; sub_plus_1 << "%"<<j<<"P";
    subst_[sub_plus_1.str()] = std::string(P);
  }
  int digits = std::floor(std::log10(double(eckit::mpi::size())))+1;
  std::stringstream p; p<<"%"<<digits<<"p";
  std::stringstream P; P<<"%"<<digits<<"P";
  subst_["%p"] = subst_[p.str()];
  subst_["%P"] = subst_[P.str()];

}

void LogFormat::beginLine()
{
  *target() << parsedPrefix();
}

void LogFormat::endLine()
{
}

const std::string& LogFormat::prefix() const
{
  return prefix_;
}

void LogFormat::setPrefix( const std::string& p )
{
  prefix_ = p;
}

std::string LogFormat::parsedPrefix() const
{
  subst_["%Y"] = TimeStamp("%Y");
  subst_["%m"] = TimeStamp("%m");
  subst_["%d"] = TimeStamp("%d");
  subst_["%H"] = TimeStamp("%H");
  subst_["%M"] = TimeStamp("%M");
  subst_["%S"] = TimeStamp("%S");
  std::string pref = prefix();
  std::map<std::string,std::string>::const_iterator it;
  for( it = subst_.begin(); it!=subst_.end(); ++it )
  {
    replace(pref, it->first, it->second);
  }
  return pref;
}


FormattedChannel::FormattedChannel( std::ostream* channel, LogFormat* format ) :
  FormatChannel( channel, format ),
  channel_(channel),
  format_(format)
{
}

FormattedChannel::FormattedChannel( std::ostream& channel, LogFormat* format ) :
  FormatChannel( channel, format ),
  channel_(&channel),
  format_(format)
{
}

FormattedChannel::~FormattedChannel()
{
}

void FormattedChannel::set_prefix( const std::string& p )
{
  format_->setPrefix(p);
}

const std::string& FormattedChannel::prefix() const
{
  return format_->prefix();
}








} // namespace atlas


