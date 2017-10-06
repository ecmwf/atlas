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
#include "eckit/log/Timer.h"

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

public: // static methods

    static std::string report();
    static std::string report( const eckit::Configuration& config );

public:

    TimerT( const eckit::CodeLocation&, const std::string& msg, std::ostream& out = Logging::channel() );
    TimerT( const eckit::CodeLocation&, std::ostream& out = Logging::channel() );

    ~TimerT();

    bool running() const;

    void start();

    void stop();

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

    mutable eckit::Timer timer_;
    eckit::CodeLocation loc_;
    std::ostream& out_;
    std::string msg_;
    Identifier id_;
    Nesting nesting_;
};

//-----------------------------------------------------------------------------------------------------------
// Definitions

template< typename TimerTraits >
inline TimerT<TimerTraits>::TimerT( const eckit::CodeLocation& loc, const std::string& msg, std::ostream& out ) :
  loc_(loc),
  out_(out),
  msg_(msg),
  nesting_(loc) {
  start();
}

template< typename TimerTraits >
inline TimerT<TimerTraits>::TimerT( const eckit::CodeLocation& loc, std::ostream& out ) :
  loc_(loc),
  out_(out),
  msg_( loc_ ? loc_.func() : "" ),
  nesting_(loc_) {
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
    id_ = Timings::add( nesting_, msg_ );
}

template< typename TimerTraits >
inline void TimerT<TimerTraits>::updateTimings() const {
    Timings::update( id_, timer_.elapsed() );
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
        out_ << msg_ << " ... done : " << timer_.elapsed() << "s" << std::endl;
    }
}

template< typename TimerTraits >
inline double TimerT<TimerTraits>::elapsed() const {
    return timer_.elapsed();
}

template< typename TimerTraits >
inline std::string TimerT<TimerTraits>::report() {
    return Timings::report();
}

template< typename TimerTraits >
inline std::string TimerT<TimerTraits>::report( const eckit::Configuration& config ) {
    return Timings::report(config);
}

//-----------------------------------------------------------------------------------------------------------

} // timer
} // runtime
} // atlas
