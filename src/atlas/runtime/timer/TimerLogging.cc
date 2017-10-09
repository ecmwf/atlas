/*
 * (C) Copyright 1996-2017 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "TimerLogging.h"

#include <iostream>
#include "eckit/log/Channel.h"
#include "atlas/library/Library.h"

//-----------------------------------------------------------------------------------------------------------

namespace atlas {
namespace runtime {
namespace timer {

class TimerLoggingState {
private:
    std::ostream* channel_;

    TimerLoggingState() {
        channel_ = &atlas::Library::instance().timer().channel();
    }

public:

    static eckit::Channel& empty_channel() {
        static eckit::Channel channel;
        return channel;
    }

    static TimerLoggingState& instance() {
        static TimerLoggingState channel;
        return channel;
    }

    operator std::ostream&() { return *channel_; }
    operator std::ostream*() { return channel_; }

    void set( std::ostream& channel ) { channel_ = &channel; }
    void set( bool state ) { if( state == false ) channel_ = &empty_channel(); }
};

TimerLogging::TimerLogging( bool state ) :
    previous_state_( TimerLoggingState::instance() ) {
    TimerLoggingState::instance().set( state );
}

TimerLogging::TimerLogging( std::ostream& channel ) :
    previous_state_( TimerLoggingState::instance() ) {
    TimerLoggingState::instance().set( channel );
}

TimerLogging::~TimerLogging() {
      TimerLoggingState::instance().set( *previous_state_ );
}

std::ostream& TimerLogging::channel() {
  return TimerLoggingState::instance();
}

std::ostream& TimerNoLogging::channel() {
  return TimerLoggingState::empty_channel();
}


} // namespace timer
} // namespace runtime
} // namespace atlas

