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

#include "atlas/library/Library.h"
#include "atlas/parallel/mpi/mpi.h"

//-----------------------------------------------------------------------------------------------------------

namespace atlas {
namespace runtime {
namespace timer {

class TimerBarriers {
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
    TimerBarriers(bool state) :
        previous_state_( State::instance() ) {
        State::instance().set(state);
    }
    ~TimerBarriers() {
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

} // namespace timer
} // namespace runtime
} // namespace atlas

