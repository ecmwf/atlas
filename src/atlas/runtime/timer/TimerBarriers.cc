/*
 * (C) Copyright 1996-2017 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include <sstream>
#include "TimerBarriers.h"
#include "atlas/library/Library.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/runtime/timer/StopWatch.h"


//-----------------------------------------------------------------------------------------------------------

namespace atlas {
namespace runtime {
namespace timer {

class TimerBarriersState {
private:
    TimerBarriersState() {
        barriers_ = atlas::Library::instance().barriers();
    }
    bool barriers_;
    StopWatch stopwatch_;
public:
    TimerBarriersState(TimerBarriersState const&)           = delete;
    void operator=(TimerBarriersState const&)  = delete;
    static TimerBarriersState& instance() {
        static TimerBarriersState state;
        return state;
    }
    operator bool() const {
        return barriers_;
    }
    void set( bool state ) {
        barriers_ = state;
    }
    StopWatch& stopwatch() { return stopwatch_; }
    
    std::string report() const {
      std::stringstream out;
      double time = stopwatch_.elapsed();
      if( time ) {
        out << "Total time spent in mpi barriers due to load imbalance : " << time << "s" << std::endl;
      }
      return out.str();
    }
};

TimerBarriers::TimerBarriers(bool state) :
  previous_state_( TimerBarriersState::instance() ) {
  TimerBarriersState::instance().set(state);
}

TimerBarriers::~TimerBarriers() {
  restore();
}

void TimerBarriers::restore() {
  TimerBarriersState::instance().set( previous_state_ );
}

bool TimerBarriers::state() {
  return TimerBarriersState::instance();
}

void TimerBarriers::execute() {
  if( state() ) {
    TimerBarriersState::instance().stopwatch().start();
    parallel::mpi::comm().barrier();
    TimerBarriersState::instance().stopwatch().stop();
  }
}

double TimerBarriers::time() {
  return TimerBarriersState::instance().stopwatch().elapsed();
}

double TimerBarriersNone::time() {
  return TimerBarriersState::instance().stopwatch().elapsed();
}

std::string TimerBarriers::report() {
  return TimerBarriersState::instance().report();
}

std::string TimerBarriersNone::report() {
  return TimerBarriersState::instance().report();
}
  
} // namespace timer
} // namespace runtime
} // namespace atlas

