/*
 * (C) Copyright 1996-2017 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include <sstream>
#include <regex>
#include <string>

#include "Timings.h"
#include "eckit/config/Configuration.h"
#include "atlas/util/detail/CallStack.h"
#include "atlas/util/Config.h"

//-----------------------------------------------------------------------------------------------------------

namespace atlas {
namespace runtime {
namespace timer {

class TimingsRegistry {

  using CallStack = Timings::CallStack;

private:

    std::vector<long>                 counts_;
    std::vector<double>               tot_timings_;
    std::vector<double>               min_timings_;
    std::vector<double>               max_timings_;
    std::vector<std::string>          titles_;
    std::vector<eckit::CodeLocation>  locations_;
    std::vector<long>                 nest_;
    std::vector<CallStack>            stack_;
    std::map<std::string,size_t> index_;

    TimingsRegistry() {
    }

public:

    static TimingsRegistry& instance() {
        static TimingsRegistry registry;
        return registry;
    }

    size_t add( const CallStack& stack, const std::string& title );

    void update( size_t idx, double seconds );

    size_t size() const;

    void report( std::ostream& out, const eckit::Configuration& config );
private:
    std::string prefix( long indent, long nest ) const;
    std::string filter_filepath( const std::string& filepath ) const;


};

size_t TimingsRegistry::add( const CallStack& stack, const std::string& title ) {
    std::stringstream ss; ss << stack;
    std::string key = ss.str();
    auto it = index_.find( key );
    if( it == index_.end() ) {
      size_t idx = size();
      index_[key] = idx;
      counts_.emplace_back( 0 );
      tot_timings_.emplace_back( 0 );
      min_timings_.emplace_back( 0 );
      max_timings_.emplace_back( 0 );
      titles_.emplace_back( title );
      locations_.emplace_back( stack.front() );
      nest_.emplace_back( stack.size() );
      stack_.emplace_back( stack );
      return idx;
    }
    else {
      return it->second;
    }
}

void TimingsRegistry::update( size_t idx, double seconds ) {
    counts_[idx] += 1;
    tot_timings_[idx] += seconds;
    min_timings_[idx] = counts_[idx] == 1 ? seconds : std::min( seconds, min_timings_[idx] );
    max_timings_[idx] = counts_[idx] == 1 ? seconds : std::max( seconds, max_timings_[idx] );
}

size_t TimingsRegistry::size() const { return counts_.size(); }

void TimingsRegistry::report( std::ostream& out, const eckit::Configuration& config ) {
    long indent   = config.getLong("indent",4);
    long max_nest = config.getLong("depth",0);
    long decimals = config.getLong("decimals",5);

    auto digits_before_decimal = [](double x) -> int {
      return std::floor(std::log10(std::trunc( std::max(1.,x)) ) )+1;
    };
    auto digits = [](long x) -> long {
      return std::floor(std::log10( std::max(1l,x) ) )+1l;
    };

    size_t max_title_length(0);
    long max_count(0);
    double max_seconds(0);
    for( size_t j=0; j<size(); ++ j ) {
      size_t nest = nest_[j];
      if( (not max_nest) or (nest <= max_nest) ) {
        max_title_length = std::max(max_title_length,titles_[j].size()+nest_[j]*indent);
        max_count        = std::max(max_count,counts_[j]);
        max_seconds      = std::max(max_seconds,tot_timings_[j]);
      }
    }
    size_t max_count_length          = digits(max_count);
    size_t max_digits_before_decimal = digits_before_decimal(max_seconds);

    auto print_time = [max_digits_before_decimal,decimals](double x) -> std::string {
      std::stringstream out;
      char unit = 's';
      if( std::floor(x) >= 60 ) {
        x/=60.;
        unit = 'm';
      }
      out << std::right
          << std::fixed
          << std::setprecision(decimals)
          << std::setw(max_digits_before_decimal+decimals+1)
          << x << unit;
      return out.str();
    };


    for( size_t j=0; j<size(); ++ j ) {

      const auto& tot    = tot_timings_[j];
      const auto& min    = min_timings_[j];
      const auto& max    = max_timings_[j];
      const auto& count  = counts_[j];
      const auto& title  = titles_[j];
      const auto& loc    = locations_[j];
      const auto& nest   = nest_[j];

      if( (not max_nest) or (nest <= max_nest) ) {

        std::string sep(" │ ");
        out << prefix(indent,nest)
            << std::left << std::setw(max_title_length-nest*indent) << title
            << sep << "count: " << std::left  << std::setw(max_count_length) << count
            << sep << "tot: "   << print_time(tot)
            << sep << "avg: "   << print_time(tot/double(count))
            << sep << "min: "   << print_time(min)
            << sep << "max: "   << print_time(max)
            << sep << filter_filepath(loc.file()) << " +"<<loc.line()
            << std::endl;

      }

    }
}

std::string TimingsRegistry::prefix( long indent, long nest ) const {
    std::stringstream out;
    std::string first;
    std::string fill;
    for( long i=0; i<nest; ++i ) {
      if( i==nest-1 ) {
        first = "├";
        fill  = "─";
      }
      else {
        first = "│";
        fill  = " ";
      }
      out << first;
      for( long j=1; j<indent; ++j )
        out << fill;
    }
    return out.str();
}

std::string TimingsRegistry::filter_filepath( const std::string& filepath ) const {
    std::regex filepath_re("(.*)?/atlas/src/(.*)");
    std::smatch matches;
    std::string filtered(filepath);
    if( std::regex_search( filepath, matches, filepath_re ) ) {
        filtered = matches[2];
    }
    return filtered;
}

Timings::Identifier Timings::add( const CallStack& stack, const std::string& title ) {
    return TimingsRegistry::instance().add( stack, title );
}

void Timings::update( const Identifier& id, double seconds ) {
    TimingsRegistry::instance().update( id, seconds );
}

std::string Timings::report() {
  return report( util::NoConfig() );
}

std::string Timings::report( const eckit::Configuration& config ) {
    std::stringstream out;
    TimingsRegistry::instance().report( out, config );
    return out.str();
}


} // namespace timer
} // namespace runtime
} // namespace atlas

