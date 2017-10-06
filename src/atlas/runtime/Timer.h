/*
 * (C) Copyright 1996-2017 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#pragma once

#include <unordered_map>
#include <regex>
#include <stack>
#include "eckit/log/Timer.h"
#include "eckit/log/Statistics.h"
#include "atlas/library/Library.h"
#include "atlas/library/config.h"
#include "atlas/runtime/Log.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/util/Config.h"

//-----------------------------------------------------------------------------------------------------------

#define __ATLAS__REVERSE 5, 4, 3, 2, 1, 0
#define __ATLAS__ARGN(_1, _2, _3, _4, _5, N, ...) N
#define __ATLAS__NARG__(dummy, ...)  __ATLAS__ARGN(__VA_ARGS__)
#define __ATLAS__NARG_(...)          __ATLAS__NARG__(dummy, ##__VA_ARGS__, __ATLAS__REVERSE)
#define __ATLAS__SPLICE(a,b)        __ATLAS__SPLICE_1(a,b)
#define __ATLAS__SPLICE_1(a,b)      __ATLAS__SPLICE_2(a,b)
#define __ATLAS__SPLICE_2(a,b)      a##b


#define __ATLAS__ARG16(_0, _1, _2, _3, _4, _5, _6, _7, _8, _9, _10, _11, _12, _13, _14, _15, ...) _15
#define __ATLAS__HAS_COMMA(...) __ATLAS__ARG16(__VA_ARGS__, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0)
#define __ATLAS__TRIGGER_PARENTHESIS(...) ,
#define __ATLAS__ISEMPTY(...)                                                                      \
__ATLAS__ISEMPTY_(                                                                                 \
  /* test if there is just one argument, eventually an empty one */                                \
  __ATLAS__HAS_COMMA(__VA_ARGS__),                                                                 \
  /* test if _TRIGGER_PARENTHESIS_ together with the argument adds a comma */                      \
  __ATLAS__HAS_COMMA(__ATLAS__TRIGGER_PARENTHESIS __VA_ARGS__),                                    \
  /* test if the argument together with a parenthesis adds a comma */                              \
  __ATLAS__HAS_COMMA(__VA_ARGS__ (/*empty*/)),                                                     \
  /* test if placing it between __ATLAS__TRIGGER_PARENTHESIS and the parenthesis adds a comma */   \
  __ATLAS__HAS_COMMA(__ATLAS__TRIGGER_PARENTHESIS __VA_ARGS__ (/*empty*/))                         \
)

#define __ATLAS__PASTE5(_0, _1, _2, _3, _4) _0 ## _1 ## _2 ## _3 ## _4
#define __ATLAS__ISEMPTY_(_0, _1, _2, _3) __ATLAS__HAS_COMMA(__ATLAS__PASTE5(__ATLAS__IS_EMPTY_CASE_, _0, _1, _2, _3))
#define __ATLAS__IS_EMPTY_CASE_0001 ,

#define __ATLAS__CALL_NARG_1(...) 0
#define __ATLAS__CALL_NARG_0 __ATLAS__NARG_
#define __ATLAS__NARG(...) __ATLAS__SPLICE( __ATLAS__CALL_NARG_, __ATLAS__ISEMPTY( __VA_ARGS__ ) ) (__VA_ARGS__)

//-----------------------------------------------------------------------------------------------------------

#if ATLAS_HAVE_TIMINGS

#define ATLAS_TIME(...)       ATLAS_TIME_( __ATLAS__NARG(__VA_ARGS__), ##__VA_ARGS__ )
#define ATLAS_TIME_SCOPE(...) ATLAS_TIME_SCOPE_( __ATLAS__NARG(__VA_ARGS__), ##__VA_ARGS__ )

#define ATLAS_TIME_SCOPE_(N, ...) __ATLAS__SPLICE( ATLAS_TIME_SCOPE_, N)(__VA_ARGS__)
#define ATLAS_TIME_SCOPE_0() \
  for( ::atlas::Timer __ATLAS__SPLICE( timer, __LINE__ ) (Here());\
    __ATLAS__SPLICE( timer, __LINE__ ) .running(); \
    __ATLAS__SPLICE( timer, __LINE__ ) .stop() )
#define ATLAS_TIME_SCOPE_1(title) \
  for( ::atlas::Timer __ATLAS__SPLICE( timer, __LINE__ ) (Here(),title);\
    __ATLAS__SPLICE( timer, __LINE__ ) .running(); \
    __ATLAS__SPLICE( timer, __LINE__ ) .stop() )

#define ATLAS_TIME_(N, ...)  __ATLAS__SPLICE( ATLAS_TIME_, N)(__VA_ARGS__)
#define ATLAS_TIME_0()       ::atlas::Timer __ATLAS__SPLICE( timer, __LINE__ ) (Here());
#define ATLAS_TIME_1(title)  ::atlas::Timer __ATLAS__SPLICE( timer, __LINE__ ) (Here(),title);

#else

#define ATLAS_TIME(...)
#define ATLAS_TIME_SCOPE(...)

#endif


namespace atlas {

class CallStack: public std::list< eckit::CodeLocation > {
private:
  using Base = std::list< eckit::CodeLocation >;
public:
  using Base::Base;

  void print(std::ostream& out) const {
    const_iterator it = begin();
    out << "\nstack:\n";
    for( size_t i=0; i<size(); ++i, ++it ) {
      out << i << ": " << it->file() << " +"<< it->line() << '\n';
    }
    out << "\n";
    out << std::flush;
  }

  friend std::ostream& operator<< ( std::ostream& out, const CallStack& p ) {
    p.print(out);
    return out;
  }

};

template< typename TimerTraits >
class TimerT {
public:
    using Barrier      = typename TimerTraits::Barrier;
    using Log          = typename TimerTraits::Log;
    using Report       = typename TimerTraits::Report;
    using Nest         = typename TimerTraits::Nest;
    using CodeLocation = eckit::CodeLocation;
    using ReportIdentifier = typename Report::Identifier;
public:

    TimerT( const CodeLocation&, const std::string& msg, std::ostream& out = Log::channel() );
    TimerT( const CodeLocation&, std::ostream& out = Log::channel() );

    ~TimerT();

    bool running() const;

    void start();

    void stop();

    double elapsed() const;

private:

    void barrier() const;

    void updateTimings() const;

    void registerTimer();

private:
    mutable eckit::Timer timer_;
    CodeLocation loc_;
    std::ostream& out_;
    std::string msg_;
    ReportIdentifier id_;
    Nest nest_;
};

// Definitions

template< typename TimerTraits >
inline TimerT<TimerTraits>::TimerT( const CodeLocation& loc, const std::string& msg, std::ostream& out ) :
  loc_(loc),
  out_(out),
  msg_(msg),
  nest_(loc) {
  start();
}

template< typename TimerTraits >
inline TimerT<TimerTraits>::TimerT( const CodeLocation& loc, std::ostream& out ) :
  loc_(loc),
  out_(out),
  msg_( loc_ ? loc_.func() : "" ),
  nest_(loc_) {
  start();
}

template< typename TimerTraits >
inline TimerT<TimerTraits>::~TimerT() {
    stop();
}

template< typename TimerTraits >
inline void TimerT<TimerTraits>::barrier() const {
    Barrier::execute();
}

template< typename TimerTraits >
inline void TimerT<TimerTraits>::registerTimer() {
    id_ = Report::add( nest_, msg_ );
}

template< typename TimerTraits >
inline void TimerT<TimerTraits>::updateTimings() const {
    Report::update( id_, timer_ );
}

template< typename TimerTraits >
inline bool TimerT<TimerTraits>::running() const {
    return timer_.running();
}

template< typename TimerTraits >
inline void TimerT<TimerTraits>::start() {
    timer_.stop();
    registerTimer();
    out_ << msg_ << " ..." << std::endl;
    barrier();
    timer_.start();
}

template< typename TimerTraits >
inline void TimerT<TimerTraits>::stop() {
    if( running() ) {
        barrier();
        timer_.stop();
        updateTimings();
        out_ << msg_ << " ... done : " << timer_.elapsed() << " seconds" << std::endl;
    }
}

template< typename TimerTraits >
inline double TimerT<TimerTraits>::elapsed() const {
    return timer_.elapsed();
}


class TimerBarrier {
private:
    class State {
    private:
        State() {
            barriers_ = atlas::Library::instance().timer().barriers();
        }
        bool barriers_;
    public:
        State(State const&)           = delete;
        void operator=(State const&)  = delete;
        static State& instance() {
            static State state;
            return state;
        }
        operator bool() const {
            return barriers_;
        }
        void set( bool state ) {
            barriers_ = state;
        }
    };

    bool previous_state_;

public:
    TimerBarrier(bool state) :
        previous_state_( State::instance() ) {
        State::instance().set(state);
    }
    ~TimerBarrier() {
        restore();
    }
    void restore() {
        State::instance().set( previous_state_ );
    }
    static bool state() { return State::instance(); }
    static void execute() {
      if( state() )
          parallel::mpi::comm().barrier();
    }
};



class TimerNest {
private:
    class State {
    private:
        State() {}
        CallStack stack_;
    public:
        State(State const&)           = delete;
        void operator=(State const&)  = delete;
        static State& instance() {
            static State state;
            return state;
        }
        operator CallStack() const {
            return stack_;
        }
        CallStack& push( const eckit::CodeLocation& loc ) {
            stack_.push_front(loc);
            return stack_;
        }
        CallStack& pop() {
            stack_.pop_front();
            return stack_;
        }
    };

    long depth_;
    CallStack stack_;

public:
    TimerNest( const eckit::CodeLocation& loc ) :
        stack_( State::instance().push( loc ) ) {
    }
    ~TimerNest() {
        State::instance().pop();
    }
    operator long() const { return stack_.size(); }
    operator CallStack() const { return stack_; }
};

class TimerLog {
private:
    class State {
    private:
        std::ostream* channel_;

        State() {
            channel_ = &atlas::Library::instance().timer().channel();
        }

    public:

        static eckit::Channel& empty_channel() {
            static eckit::Channel channel;
            return channel;
        }

        static State& instance() {
            static State channel;
            return channel;
        }

        operator std::ostream&() { return *channel_; }
        operator std::ostream*() { return channel_; }

        void set( std::ostream& channel ) { channel_ = &channel; }
        void set( bool state ) { if( state == false ) channel_ = &empty_channel(); }
    };

    std::ostream* previous_state_;

public:
    TimerLog( bool state ) :
        previous_state_( State::instance() ) {
        State::instance().set( state );
    }
    TimerLog( std::ostream& channel ) :
        previous_state_( State::instance() ) {
        State::instance().set( channel );
    }
    ~TimerLog() {
          State::instance().set( *previous_state_ );
    }
    static std::ostream& channel() { return State::instance(); }
};

class TimerReport {
  using Timing = eckit::Timing;
private:
    class Registry {
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

        Registry() {
        }

    public:

        static Registry& instance() {
            static Registry registry;
            return registry;
        }

        size_t add( const CallStack& stack, const std::string& title ) {
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

        void update( size_t idx, eckit::Timer& timer ) {
            double time = timer.elapsed();
            counts_[idx] += 1;
            tot_timings_[idx] += time;
            min_timings_[idx] = counts_[idx] == 1 ? time : std::min( time, min_timings_[idx] );
            max_timings_[idx] = counts_[idx] == 1 ? time : std::max( time, max_timings_[idx] );
        }

        size_t size() const { return counts_.size(); }

        void report( std::ostream& out, const eckit::Configuration& config ) {
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

                std::string sep("  |  ");
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
    private:
        std::string prefix( long indent, long nest ) {
            std::stringstream out;
            switch( indent ) {
              case 2:
                for( long i=0; i<nest; ++i ) {
                  out << ( i==nest-1 ? "|-" : "| " );
                }
                break;
              case 4:
              default:
                for( long i=0; i<nest; ++i ) {
                  out << ( i==nest-1 ? "  |-" : "  | " );
                }
                break;
            }
            return out.str();
        }

        std::string filter_filepath( const std::string& filepath ) {
          std::regex filepath_re("(.*)?/atlas/src/(.*)");
          std::smatch matches;
          std::string filtered(filepath);
          if( std::regex_search( filepath, matches, filepath_re ) ) {
            filtered = matches[2];
          }
          return filtered;
        }

    };


public:
    using Identifier = size_t;

    static Identifier add( const CallStack& stack, const std::string& title ) {
      return Registry::instance().add( stack, title );
    }
    static void update( const Identifier& id, eckit::Timer& timer ) {
      Registry::instance().update( id, timer );
    }
    static void report( std::ostream& out, const eckit::Configuration& config = util::NoConfig() ) {
      Registry::instance().report( out, config );
    }

    static void report( const eckit::Configuration& config = util::NoConfig() ) {
      report( Log::info(), config );
    }

};

struct TimerTraits {
    using Barrier = TimerBarrier;
    using Log     = TimerLog;
    using Report  = TimerReport;
    using Nest    = TimerNest;
};

class Timer : public TimerT< TimerTraits > {
public:
    using TimerT::TimerT;
};

} // namespace atlas
