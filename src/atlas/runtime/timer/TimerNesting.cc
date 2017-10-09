/*
 * (C) Copyright 1996-2017 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "TimerNesting.h"

//-----------------------------------------------------------------------------------------------------------

namespace atlas {
namespace runtime {
namespace timer {
  

class TimerNestingState {
  using CallStack = util::detail::CallStack;
private:
  TimerNestingState() {}
  CallStack stack_;
public:
  TimerNestingState(TimerNestingState const&)           = delete;
  void operator=(TimerNestingState const&)  = delete;
  static TimerNestingState& instance() {
    static TimerNestingState state;
    return state;
  }
  operator CallStack() const {
    return stack_;
  }
  CallStack& push( const eckit::CodeLocation& loc ) {
    stack_.push_front(loc);
    return stack_;
  }
  void pop() {
    stack_.pop_front();
  }
};

TimerNesting::TimerNesting( const eckit::CodeLocation& loc ) :
  loc_(loc),
  stack_( TimerNestingState::instance().push( loc ) ) {
}

TimerNesting::~TimerNesting() {
  if( running_ )
    stop();
}

void TimerNesting::stop() {
  if( running_ )
    TimerNestingState::instance().pop();
  running_ = false;
}

void TimerNesting::start() {
  if( not running_ )
    TimerNestingState::instance().push( loc_ );
  running_ = true;
}

} // namespace timer
} // namespace runtime
} // namespace atlas

