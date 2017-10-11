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
#include "atlas/runtime/timer/Timings.h"
#include "atlas/runtime/timer/TimerNesting.h"

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
    using Tracing          = typename TimerTraits::Tracing;
    using Labels           = std::vector<std::string>;

public: // static methods

    static std::string report();
    static std::string report( const eckit::Configuration& config );

public:

    TimerT( const eckit::CodeLocation& );
    TimerT( const eckit::CodeLocation&, const std::string& title );
    TimerT( const eckit::CodeLocation&, const std::string& title, const Labels& );

    ~TimerT();

    bool running() const;

    void start();

    void stop();

    void pause();

    void resume();

    double elapsed() const;

private: // types

    using Nesting    = timer::TimerNesting;
    using Timings    = timer::Timings;
    using Identifier = timer::Timings::Identifier;

private: // member functions

    void barrier() const;

    void updateTimings() const;

    void registerTimer();

private: // member data

    bool running_{true};
    StopWatch stopwatch_;
    eckit::CodeLocation loc_;
    std::string title_;
    Identifier id_;
    Nesting nesting_;
    Labels labels_;
};

//-----------------------------------------------------------------------------------------------------------
// Definitions

template< typename TimerTraits >
inline TimerT<TimerTraits>::TimerT( const eckit::CodeLocation& loc, const std::string& title ) :
  loc_(loc),
  title_(title),
  nesting_(loc) {
  start();
}

template< typename TimerTraits >
inline TimerT<TimerTraits>::TimerT( const eckit::CodeLocation& loc ) :
  loc_(loc),
  title_( loc_ ? loc_.func() : ""),
  nesting_(loc_) {
  start();
}

template< typename TimerTraits >
inline TimerT<TimerTraits>::TimerT( const eckit::CodeLocation& loc, const std::string& title, const Labels& labels ) :
  loc_(loc),
  title_(title),
  nesting_(loc),
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
    id_ = Timings::add( loc_, nesting_, title_ + (Barriers::state() ? " [b]" : ""), labels_ );
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
    Tracing::start( title_ );
    barrier();
    stopwatch_.start();
}

template< typename TimerTraits >
inline void TimerT<TimerTraits>::stop() {
    if( running_ ) {
        barrier();
        stopwatch_.stop();
        nesting_.stop();
        updateTimings();
        Tracing::stop( title_, stopwatch_.elapsed() );
        running_ = false;
    }
}

template< typename TimerTraits >
inline void TimerT<TimerTraits>::pause() {
    if( running_ ) {
      barrier();
      stopwatch_.stop();
      nesting_.stop();
    }
}

template< typename TimerTraits >
inline void TimerT<TimerTraits>::resume() {
    if( running_ ) {
      barrier();
      nesting_.start();
      stopwatch_.start();
    }
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
