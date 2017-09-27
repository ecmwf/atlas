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

#include "eckit/log/Timer.h"
#include "atlas/library/Library.h"
#include "atlas/runtime/Log.h"
#include "atlas/parallel/mpi/mpi.h"

namespace atlas {

template< typename TimerTraits >
class TimerT {
public:
    using Barrier = typename TimerTraits::Barrier;
    using Log     = typename TimerTraits::Log;

public:

    TimerT(const std::string& msg, std::ostream& out = Log::channel() );

    TimerT();

    ~TimerT();

    bool running() const;

    void start();

    void stop();

    double elapsed() const;

private:

    void barrier() const;

    static eckit::Channel& empty_channel();

private:
    mutable eckit::Timer timer_;
    std::ostream& out_;
    std::string msg_;
    bool barrier_;
};

// Definitions

template< typename TimerTraits >
inline TimerT<TimerTraits>::TimerT(const std::string& msg, std::ostream& out ) :
  out_(out),
  msg_(msg),
  barrier_(Barrier::state()) {
  start();
}

template< typename TimerTraits >
inline TimerT<TimerTraits>::TimerT() :
    TimerT(std::string(), empty_channel() ) {
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
inline bool TimerT<TimerTraits>::running() const {
    return timer_.running();
}

template< typename TimerTraits >
inline void TimerT<TimerTraits>::start() {
    timer_.stop();
    out_ << msg_ << " ..." << std::endl;
    barrier();
    timer_.start();
}

template< typename TimerTraits >
inline void TimerT<TimerTraits>::stop() {
    if( running() ) {
        barrier();
        timer_.stop();
        out_ << msg_ << " ... done : " << timer_.elapsed() << " seconds" << std::endl;
    }
}

template< typename TimerTraits >
inline double TimerT<TimerTraits>::elapsed() const {
    return timer_.elapsed();
}

template< typename TimerTraits >
inline eckit::Channel& TimerT<TimerTraits>::empty_channel() {
    static eckit::Channel channel;
    return channel;
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

struct TimerTraits {
    using Barrier = TimerBarrier;
    using Log     = TimerLog;
};

class Timer : public TimerT< TimerTraits > {
public:
    using TimerT::TimerT;
};

} // namespace atlas
