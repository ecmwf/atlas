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
public:
    using CallStack = util::detail::CallStack;
    TimerNesting( const eckit::CodeLocation& );
    ~TimerNesting();
    operator long() const { return stack_.size(); }
    operator CallStack() const { return stack_; }
    void stop();
    void start();
private:
    CallStack stack_;
    eckit::CodeLocation loc_;
    bool running_{true};
};


} // namespace timer
} // namespace runtime
} // namespace atlas

