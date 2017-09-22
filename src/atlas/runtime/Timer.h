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
	using Barriers = typename TimerTraits::Barriers;
	using Log      = typename TimerTraits::Log;

private:

    TimerT(const std::string& msg, std::ostream* _out) : 
        out_(_out),
        finished_(false),
        msg_(msg) {
        timer_.stop();
        out() << msg << " ..." << std::endl;
        barrier();
        timer_.start();
    }

public:

    TimerT(const std::string& msg, std::ostream& out) : 
        TimerT(msg,&out){
    }

    TimerT(const std::string& msg) : 
    	TimerT(msg,nullptr) {
    }
    
    TimerT() : 
        TimerT(std::string(),&empty_channel()) {
    }

    ~TimerT() {
        finish();
    }

    static eckit::Channel& empty_channel() {
        static eckit::Channel channel;
        return channel;
    }

    std::ostream& out() const {
    	if( out_ )
    		return *out_;
    	else
    		return Log::channel();
    }

    void barrier() const {
        Barriers::execute();
    }
    
    bool finished() const { return finished_; }
    
    void finish() { 
        if( not finished_ ) {
            barrier();
            timer_.stop();
            out() << msg_ << " ... done : " << timer_.elapsed() << " seconds" << std::endl; 
            finished_ = true;
        }
    }
    
    double elapsed() const { return timer_.elapsed(); }

private:
  mutable eckit::Timer timer_;
  std::ostream* out_{nullptr};
  bool finished_;
  std::string msg_;
};


namespace {

class TimerBarriers {
private:
    class State {
    private:
	    State() {
            barriers_ = atlas::Library::instance().timer().barriers();
	    }
	    bool barriers_;
    public:
	    static State& instance() {
		    static State state;
		    return state;
        }
	    operator bool() const { return barriers_; }
	    void set(bool state) { barriers_ = state; }
    };

    bool previous_state_;

public:
    TimerBarriers(bool state) :
        previous_state_( State::instance() ) {
  	    State::instance().set(state);
    }
    ~TimerBarriers() {
  	    State::instance().set(previous_state_);
    }
    static void execute() {
        if( State::instance() ) {
            parallel::mpi::comm().barrier();
        }
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
	    
	    void set(std::ostream& channel ) { channel_ = &channel; }
    };

    std::ostream* previous_state_;

public:
    TimerLog(bool state) :
        previous_state_( State::instance() ) {
        if( state == false )
  	        State::instance().set( State::empty_channel() );
    }
    TimerLog(std::ostream& channel) :
        previous_state_( State::instance() ) {
        State::instance().set( channel );
    }
    ~TimerLog() {
  	    State::instance().set( *previous_state_ );
    }
    static std::ostream& channel() { return State::instance(); }
};

struct TimerTraits {
    using Barriers = TimerBarriers;
    using Log      = TimerLog;
};

}

using Timer = TimerT< TimerTraits >;

} // namespace atlas
