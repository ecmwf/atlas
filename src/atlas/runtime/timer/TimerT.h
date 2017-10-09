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

#include <iosfwd>
#include <string>
#include <vector>
#include <string>
#include "atlas/runtime/timer/StopWatch.h"

//-----------------------------------------------------------------------------------------------------------

namespace eckit {
  class CodeLocation;
  class Configuration;
}

//-----------------------------------------------------------------------------------------------------------

namespace atlas {
namespace runtime {
namespace timer {

//-----------------------------------------------------------------------------------------------------------

template< typename TimerTraits >
class TimerT {
public:
    using Barriers         = typename TimerTraits::Barriers;
    using Logging          = typename TimerTraits::Logging;
    using Labels           = std::vector<std::string>;

public: // static methods

    static std::string report();
    static std::string report( const eckit::Configuration& config );

public:

    TimerT( const eckit::CodeLocation&, const std::string& title, std::ostream& out = Logging::channel() );
    TimerT( const eckit::CodeLocation&, std::ostream& out = Logging::channel() );

    TimerT( const eckit::CodeLocation&, const std::string& title, const Labels&, std::ostream& out = Logging::channel() );
    TimerT( const eckit::CodeLocation&, const Labels&, std::ostream& out = Logging::channel() );

    ~TimerT();

    bool running() const;

    void start();

    void stop();

    void pause();

    void resume();

    double elapsed() const;

private: // types

    using Nesting    = typename TimerTraits::Nesting;
    using Timings    = typename TimerTraits::Timings;
    using Identifier = typename Timings::Identifier;

private: // member functions

    void barrier() const;

    void updateTimings() const;

    void registerTimer();

private: // member data

    bool running_{true};
    StopWatch stopwatch_;
    eckit::CodeLocation loc_;
    std::ostream& out_;
    std::string title_;
    Identifier id_;
    Nesting nesting_;
    Labels labels_;
};

//-----------------------------------------------------------------------------------------------------------
// Definitions

template< typename TimerTraits >
inline TimerT<TimerTraits>::TimerT( const eckit::CodeLocation& loc, const std::string& title, std::ostream& out ) :
  loc_(loc),
  out_(out),
  title_(title),
  nesting_(loc) {
  start();
}

template< typename TimerTraits >
inline TimerT<TimerTraits>::TimerT( const eckit::CodeLocation& loc, std::ostream& out ) :
  loc_(loc),
  out_(out),
  title_( loc_ ? loc_.func() : ""),
  nesting_(loc_) {
  start();
}

template< typename TimerTraits >
inline TimerT<TimerTraits>::TimerT( const eckit::CodeLocation& loc, const std::string& title, const Labels& labels, std::ostream& out ) :
  loc_(loc),
  out_(out),
  title_(title),
  nesting_(loc),
  labels_(labels) {
  start();
}

template< typename TimerTraits >
inline TimerT<TimerTraits>::TimerT( const eckit::CodeLocation& loc, const Labels& labels, std::ostream& out ) :
  loc_(loc),
  out_(out),
  title_( loc_ ? loc_.func() : "" ),
  nesting_(loc_),
  labels_(labels) {
  start();
}

template< typename TimerTraits >
inline TimerT<TimerTraits>::~TimerT() {
    stop();
}

template< typename TimerTraits >
inline void TimerT<TimerTraits>::barrier() const {
    Barriers::execute();
}

template< typename TimerTraits >
inline void TimerT<TimerTraits>::registerTimer() {
    id_ = Timings::add( nesting_, title_ + (Barriers::state() ? " [b]" : ""), labels_ );
}

template< typename TimerTraits >
inline void TimerT<TimerTraits>::updateTimings() const {
  Timings::update( id_, stopwatch_.elapsed() );
}

template< typename TimerTraits >
inline bool TimerT<TimerTraits>::running() const {
    return running_;
}

template< typename TimerTraits >
inline void TimerT<TimerTraits>::start() {
    registerTimer();
    out_ << title_ << " ..." << std::endl;
    barrier();
    stopwatch_.start();
}

template< typename TimerTraits >
inline void TimerT<TimerTraits>::stop() {
    if( running() ) {
        barrier();
        stopwatch_.stop();
        nesting_.stop();
        updateTimings();
        out_ << title_ << " ... done : " << stopwatch_.elapsed() << "s" << std::endl;
        running_ = false;
    }
}

template< typename TimerTraits >
inline void TimerT<TimerTraits>::pause() {
    if( running() ) {
      barrier();
      stopwatch_.stop();
      nesting_.stop();
    }
}

template< typename TimerTraits >
inline void TimerT<TimerTraits>::resume() {
    barrier();
    nesting_.start();
    stopwatch_.start();
}

template< typename TimerTraits >
inline double TimerT<TimerTraits>::elapsed() const {
    return stopwatch_.elapsed();
}

template< typename TimerTraits >
inline std::string TimerT<TimerTraits>::report() {
    return Timings::report() + Barriers::report();
}

template< typename TimerTraits >
inline std::string TimerT<TimerTraits>::report( const eckit::Configuration& config ) {
  return Timings::report(config) + Barriers::report();
}

//-----------------------------------------------------------------------------------------------------------

} // timer
} // runtime
} // atlas
