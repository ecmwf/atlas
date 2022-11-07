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

void CallStack::push(const CodeLocation& loc, const std::string& id) {
    if (stack_.size() == size_) {
        stack_.resize(2 * size_);
    }
    stack_[size_++] = std::hash<std::string>{}(loc.asString() + id);
}

void CallStack::pop() {
    --size_;
}

size_t CallStack::hash() const {
    if (hash_) {
        return hash_;
    }
    for (long i = size_ - 1; i >= 0; --i) {
        auto h = stack_[i];
        hash_ ^= (h << 1);
    }
    return hash_;
}

}  // namespace trace
}  // namespace runtime
}  // namespace atlas
