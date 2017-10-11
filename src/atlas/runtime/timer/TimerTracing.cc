/*
 * (C) Copyright 1996-2017 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "TimerTracing.h"

#include <iostream>
#include "eckit/log/Channel.h"
#include "atlas/library/Library.h"

//-----------------------------------------------------------------------------------------------------------

namespace atlas {
namespace runtime {
namespace timer {

//-----------------------------------------------------------------------------------------------------------

class TimerTracingState {
private:
    std::ostream* channel_;

    TimerTracingState() {
        channel_ = &atlas::Library::instance().traceChannel();
    }

public:

    static eckit::Channel& empty_channel() {
        static eckit::Channel channel;
        return channel;
    }

    static TimerTracingState& instance() {
        static TimerTracingState channel;
        return channel;
    }

    operator std::ostream&() { return *channel_; }
    operator std::ostream*() { return channel_; }

    void set( std::ostream& channel ) { channel_ = &channel; }
    void set( bool state ) { if( state == false ) channel_ = &empty_channel(); }
};

//-----------------------------------------------------------------------------------------------------------

TimerTracing::TimerTracing( bool state ) :
    previous_state_( TimerTracingState::instance() ) {
    TimerTracingState::instance().set( state );
}

TimerTracing::TimerTracing( std::ostream& channel ) :
    previous_state_( TimerTracingState::instance() ) {
    TimerTracingState::instance().set( channel );
}

TimerTracing::~TimerTracing() {
      TimerTracingState::instance().set( *previous_state_ );
}

std::ostream& TimerTracing::channel() {
  return TimerTracingState::instance();
}

bool TimerTracing::enabled() {
  return TimerTracingState::instance();
}
void TimerTracing::start( const std::string& title ) {
  if( enabled() )
    channel() << title << " ..." << std::endl;
}

void TimerTracing::stop( const std::string& title, double seconds ) {
  if( enabled() )
    channel() << title << " ... done : " << seconds << "s" << std::endl;
}
//-----------------------------------------------------------------------------------------------------------

std::ostream& TimerTracingNone::channel() {
  return TimerTracingState::empty_channel();
}

//-----------------------------------------------------------------------------------------------------------

void TimerTracingResult::stop( const std::string& title, double seconds ) {
  if( enabled() )
    channel() << title << " : " << seconds << "s" << std::endl;
}

//-----------------------------------------------------------------------------------------------------------

} // namespace timer
} // namespace runtime
} // namespace atlas

