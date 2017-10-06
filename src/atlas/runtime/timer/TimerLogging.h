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

#include <iostream>
#include "eckit/log/Channel.h"
#include "atlas/library/Library.h"

//-----------------------------------------------------------------------------------------------------------

namespace atlas {
namespace runtime {
namespace timer {


class TimerLogging {
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
    TimerLogging( bool state ) :
        previous_state_( State::instance() ) {
        State::instance().set( state );
    }
    TimerLogging( std::ostream& channel ) :
        previous_state_( State::instance() ) {
        State::instance().set( channel );
    }
    ~TimerLogging() {
          State::instance().set( *previous_state_ );
    }
    static std::ostream& channel() { return State::instance(); }
};

} // namespace timer
} // namespace runtime
} // namespace atlas

