/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor
 * does it submit to any jurisdiction.
 */

#include "Nesting.h"

//-----------------------------------------------------------------------------------------------------------

namespace atlas {
namespace runtime {
namespace trace {

class NestingState {
private:
    NestingState() {}
    CallStack stack_;

public:
    NestingState( NestingState const& ) = delete;
    void operator=( NestingState const& ) = delete;
    static NestingState& instance() {
        static NestingState state;
        return state;
    }
    operator CallStack() const { return stack_; }
    CallStack& push( const eckit::CodeLocation& loc ) {
        stack_.push_front( loc );
        return stack_;
    }
    void pop() { stack_.pop_front(); }
};

Nesting::Nesting( const eckit::CodeLocation& loc ) : loc_( loc ), stack_( NestingState::instance().push( loc ) ) {}

Nesting::~Nesting() {
    stop();
}

void Nesting::stop() {
    if ( running_ ) {
        NestingState::instance().pop();
        running_ = false;
    }
}

void Nesting::start() {
    if ( not running_ ) {
        NestingState::instance().push( loc_ );
        running_ = true;
    }
}

}  // namespace trace
}  // namespace runtime
}  // namespace atlas
