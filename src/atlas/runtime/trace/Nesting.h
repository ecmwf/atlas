/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#pragma once

#include "atlas/runtime/trace/CallStack.h"
#include "atlas/runtime/trace/CodeLocation.h"
#include "atlas/runtime/trace/Logging.h"

//-----------------------------------------------------------------------------------------------------------

namespace atlas {
namespace runtime {
namespace trace {

class CurrentCallStack {
private:
    CurrentCallStack() {}
    CallStack stack_;

public:
    CurrentCallStack(CurrentCallStack const&) = delete;
    CurrentCallStack& operator=(CurrentCallStack const&) = delete;
    static CurrentCallStack& instance() {
        static CurrentCallStack state;
        return state;
    }
    operator CallStack() const { return stack_; }
    CallStack& push(const CodeLocation& loc, const std::string& id) {
        if (Control::enabled())
            stack_.push(loc, id);
        return stack_;
    }
    void pop() {
        if (Control::enabled())
            stack_.pop();
    }
};

}  // namespace trace
}  // namespace runtime
}  // namespace atlas
