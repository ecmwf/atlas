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

#include "atlas/util/detail/CallStack.h"

//-----------------------------------------------------------------------------------------------------------

namespace atlas {
namespace runtime {
namespace timer {


class TimerNesting {
private:
    using CallStack = util::detail::CallStack;
    class State {
    private:
        State() {}
        CallStack stack_;
    public:
        State(State const&)           = delete;
        void operator=(State const&)  = delete;
        static State& instance() {
            static State state;
            return state;
        }
        operator CallStack() const {
            return stack_;
        }
        CallStack& push( const eckit::CodeLocation& loc ) {
            stack_.push_front(loc);
            return stack_;
        }
        CallStack& pop() {
            stack_.pop_front();
            return stack_;
        }
    };

    long depth_;
    CallStack stack_;

public:
    TimerNesting( const eckit::CodeLocation& loc ) :
        stack_( State::instance().push( loc ) ) {
    }
    ~TimerNesting() {
        State::instance().pop();
    }
    operator long() const { return stack_.size(); }
    operator CallStack() const { return stack_; }
};


} // namespace timer
} // namespace runtime
} // namespace atlas

