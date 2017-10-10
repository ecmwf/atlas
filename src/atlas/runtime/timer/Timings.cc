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
#include <cmath>
#include <limits>

#include "Timings.h"
#include "eckit/config/Configuration.h"
#include "eckit/filesystem/PathName.h"
#include "atlas/util/detail/CallStack.h"
#include "atlas/util/Config.h"
#include "atlas/runtime/Debug.h"
#include "atlas/parallel/mpi/mpi.h"

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
    std::vector<double>               var_timings_;
    std::vector<std::string>          titles_;
    std::vector<eckit::CodeLocation>  locations_;
    std::vector<long>                 nest_;
    std::vector<CallStack>            stack_;
    std::map<size_t,size_t> index_;

    std::map<std::string, std::vector<size_t> > labels_;

    TimingsRegistry() {
    }

public:

    static TimingsRegistry& instance() {
        static TimingsRegistry registry;
        return registry;
    }

    size_t add( const CallStack& stack, const std::string& title, const Timings::Labels& );

    void update( size_t idx, double seconds );

    size_t size() const;

    void report( std::ostream& out, const eckit::Configuration& config );
private:

    std::string filter_filepath( const std::string& filepath ) const;


};

size_t TimingsRegistry::add( const CallStack& stack, const std::string& title, const Timings::Labels& labels ) {

    size_t key = stack.hash();
    auto it = index_.find( key );
    if( it == index_.end() ) {
      size_t idx = size();
      index_[key] = idx;
      counts_.emplace_back( 0 );
      tot_timings_.emplace_back( 0 );
      min_timings_.emplace_back( std::numeric_limits<double>::max() );
      max_timings_.emplace_back( 0 );
      var_timings_.emplace_back( 0 );
      titles_.emplace_back( title );
      locations_.emplace_back( stack.loc() );
      nest_.emplace_back( stack.size() );
      stack_.emplace_back( stack );

      for( const auto& label: labels ) {
        labels_[label].emplace_back(idx);
      }

      return idx;
    }
    else {
      return it->second;
    }
}

void TimingsRegistry::update( size_t idx, double seconds ) {

    auto sqr = [](double x) { return x*x; };
    double n = counts_[idx]+1;
    double avg_nm1 = tot_timings_[idx] / std::max(n,1.);
    double var_nm1 = var_timings_[idx];
    var_timings_[idx] = counts_[idx] == 0 ? 0. : (n-2.)/(n-1.) * var_nm1 + 1./n * sqr(seconds-avg_nm1);
    min_timings_[idx] = std::min( seconds, min_timings_[idx] );
    max_timings_[idx] = std::max( seconds, max_timings_[idx] );
    tot_timings_[idx] += seconds;
    counts_[idx]      += 1;
}

size_t TimingsRegistry::size() const { return counts_.size(); }

void TimingsRegistry::report( std::ostream& out, const eckit::Configuration& config ) {
    long indent   = config.getLong("indent",2);
    long depth = config.getLong("depth",0);
    long decimals = config.getLong("decimals",5);
    bool header   = config.getBool("header",true);
    std::vector<std::string> excluded_labels_vector = config.getStringVector("exclude",std::vector<std::string>());
    std::vector<std::string> include_back;

    for( auto& label : excluded_labels_vector ) {
      size_t found = label.find("/*");
      if (found!=std::string::npos) {
        label.erase(found,2);
        include_back.push_back(label);
      }
    }

    std::set<std::string> excluded_labels(excluded_labels_vector.begin(),excluded_labels_vector.end());

    auto digits_before_decimal = [](double x) -> int {
      return std::floor(std::log10(std::trunc( std::max(1.,x)) ) )+1;
    };
    auto digits = [](long x) -> long {
      return std::floor(std::log10( std::max(1l,x) ) )+1l;
    };

    std::vector<size_t> excluded_timers_vector;
    for( auto label: labels_ ) {
      auto name = label.first;
      if( excluded_labels.count(name) ) {
        auto timers = label.second;
        for( size_t j : timers ) {
          excluded_timers_vector.push_back(j);
        }
      }
    }
    std::set<size_t> excluded_timers(excluded_timers_vector.begin(),excluded_timers_vector.end());


    auto excluded = [&](size_t i) -> bool {
      if( depth and nest_[i] > depth )
        return true;
      return excluded_timers.count(i);
    };

    std::vector<long> excluded_nest_stored(size());
    long excluded_nest=size();
    for( size_t j=0; j<size(); ++j ) {
      if( nest_[j] > excluded_nest ) {
        excluded_timers.insert(j);
      }
      if( not excluded(j) ) {
        excluded_nest = nest_[j]+1;
      } else {
        excluded_nest = std::min(excluded_nest,nest_[j]);
      }
      excluded_nest_stored[j] = excluded_nest;
    }
    for( auto& label: include_back ) {
      auto timers = labels_[label];
      for( size_t j : timers ) {
        if( nest_[j] == excluded_nest_stored[j] )
          excluded_timers.erase(j);
      }
    }




    size_t max_title_length(0);
    size_t max_location_length(0);
    size_t max_nest(0);
    long max_count(0);
    double max_seconds(0);
    for( size_t j=0; j<size(); ++ j ) {
      size_t nest = nest_[j];
      max_nest = std::max(max_nest,nest);
      if( not excluded(j) ) {
        const auto& loc    = locations_[j];
        max_title_length = std::max(max_title_length,titles_[j].size()+nest_[j]*indent);
        max_count        = std::max(max_count,counts_[j]);
        max_seconds      = std::max(max_seconds,tot_timings_[j]);
        size_t location_length = filter_filepath(loc.file()).size() + 2 + digits(loc.line());
        max_location_length = std::max(max_location_length,location_length);
      }
    }
    size_t max_count_length = digits(max_count);
    if( header ) {
      max_count_length = std::max(std::string("cnt").size(),max_count_length);
    }
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

    auto print_line = [](size_t length) -> std::string {
      std::stringstream ss;
      for( size_t i=0; i<length; ++i ) {
        ss << "─";
      }
      return ss.str();
    };
    auto print_horizontal = [&](const std::string& sep) -> std::string {
      std::stringstream ss;
      ss  << print_line(max_title_length  + digits(size()) + 3)
          << sep << print_line(max_count_length)
          << sep << print_line(max_digits_before_decimal+decimals+2)
          << sep << print_line(max_digits_before_decimal+decimals+2)
          << sep << print_line(max_digits_before_decimal+decimals+2)
          << sep << print_line(max_digits_before_decimal+decimals+2)
          << sep << print_line(max_digits_before_decimal+decimals+2)
          << sep << print_line(max_location_length);
      return ss.str();
    };

    std::string sept("─┬─");
    std::string seph("─┼─");
    std::string sep (" │ ");
    std::string sepf("─┴─");


    out << print_horizontal( sept ) << std::endl;
    out << std::left << std::setw(max_title_length + digits(size()) + 3) << "Timers"
        << sep << std::setw(max_count_length) << "cnt"
        << sep << std::setw(max_digits_before_decimal+decimals+2ul) << "tot"
        << sep << std::setw(max_digits_before_decimal+decimals+2ul) << "avg"
        << sep << std::setw(max_digits_before_decimal+decimals+2ul) << "std"
        << sep << std::setw(max_digits_before_decimal+decimals+2ul) << "min"
        << sep << std::setw(max_digits_before_decimal+decimals+2ul) << "max"
        << sep << "location"
        << std::endl;
    out << print_horizontal( seph ) << std::endl;


    std::vector<std::string> prefix_(size());
    if( indent ) {
      std::vector<bool> active(max_nest,false);
      for( long j=long(size())-1; j>=0; --j ) {
        const auto& nest   = nest_[j];

        const CallStack& this_stack = stack_[j];
        const CallStack& next_stack = (j==size()-1) ? this_stack : stack_[j+1];

        auto this_it = this_stack.rbegin();
        auto next_it = next_stack.rbegin();
        for( size_t i=0; this_it!=this_stack.rend() && next_it!=next_stack.rend(); ++i, ++this_it, ++next_it ) {
          if( *this_it == *next_it ) {
            active[i] = active[i] or false;
          } else {
            active[i] = true;
          }
        }
        for( size_t i=nest; i<active.size(); ++i ) {
          active[i] = false;
        }

        std::stringstream out;
        for( long i=0; i<nest-1; ++i ) {
          if(active[i])
            out << "│";
          else
            out << " ";
          for( size_t j=1; j<indent; ++j )
            out << " ";
        }
        if( active[nest-1] )
          out << "├";
        else
          out << "└";
        for( size_t j=1; j<indent; ++j )
          out << "─";


        prefix_[j] = out.str();
      }
    }

    for( size_t j=0; j<size(); ++j ) {

      auto& tot    = tot_timings_[j];
      auto& min    = min_timings_[j];
      auto& max    = max_timings_[j];
      auto& count  = counts_[j];
      auto& title  = titles_[j];
      auto& loc    = locations_[j];
      auto& nest   = nest_[j];
      auto  std    = std::sqrt( var_timings_[j] );
      auto  avg    = tot/double(count);

      //parallel::mpi::comm().allReduceInPlace(min,eckit::mpi::min());
      //parallel::mpi::comm().allReduceInPlace(max,eckit::mpi::max());

      if( not excluded(j) ) {

        out << std::setw(digits(size())) << j << " : " << prefix_[j] //prefix(indent,nest,next_nest)
            << std::left << std::setw(max_title_length-nest*indent) << title
            << sep << std::string(header ? "" : "count: ") << std::left  << std::setw(max_count_length) << count
            << sep << std::string(header ? "" : "tot: "  ) << print_time(tot)
            << sep << std::string(header ? "" : "avg: "  ) << print_time(avg)
            << sep << std::string(header ? "" : "std: "  ) << print_time(std)
            << sep << std::string(header ? "" : "min: "  ) << print_time(min)
            << sep << std::string(header ? "" : "max: "  ) << print_time(max)
            << sep << filter_filepath(loc.file()) << " +"<<loc.line()
            << std::endl;

      }

    }

    out << print_horizontal(sepf) << std::endl;

    std::string sepc = "───";
    // std::string sep (" │ ");

    out << "\n";
    out << print_horizontal(sepc) << std::endl;
    out << "Timers accumulated by label" << std::endl;
    out << print_horizontal(sepc) << std::endl;
    for( auto label : labels_ ) {
      auto name = label.first;
      auto timers = label.second;
      double tot(0);
      double count(0);
      for( size_t j : timers ) {
        tot += tot_timings_[j];
        count += counts_[j];
      }
      out << std::left << std::setw(40) << name
          << sep << print_time(tot)
          << sep << std::string(header ? "" : "count: ") << std::left  << std::setw(max_count_length) << count
          << std::endl;
    }
    out << print_horizontal(sepc) << std::endl;
}


std::string TimingsRegistry::filter_filepath( const std::string& filepath ) const {
    std::regex filepath_re("(.*)?/atlas/src/(.*)");
    std::smatch matches;
    std::string filtered("");
    bool is_atlas = false;
    if( std::regex_search( filepath, matches, filepath_re ) ) {
        // filtered = matches[2];
        filtered = "[atlas] ";
    }
    filtered += eckit::PathName(filepath).baseName();
    return filtered;
    //
    // return filtered;
  // return filepath;
}

Timings::Identifier Timings::add( const CallStack& stack, const std::string& title, const Labels& labels ) {
    return TimingsRegistry::instance().add( stack, title, labels );
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

