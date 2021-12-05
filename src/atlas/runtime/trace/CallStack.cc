/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "CallStack.h"

#include <functional>

#include "atlas/runtime/trace/CodeLocation.h"

namespace atlas {
namespace runtime {
namespace trace {

void CallStack::push_front(const CodeLocation& loc, const std::string& id) {
    stack_.push_front(std::hash<std::string>{}(loc.asString() + id));
}

void CallStack::pop_front() {
    stack_.pop_front();
}

size_t CallStack::hash() const {
    if (hash_) {
        return hash_;
    }
    for (auto h : stack_) {
        hash_ ^= (h << 1);
    }
    return hash_;
}

}  // namespace trace
}  // namespace runtime
}  // namespace atlas
