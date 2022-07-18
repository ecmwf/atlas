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

#include <cstddef>
#include <string>
#include <vector>

namespace atlas {
class CodeLocation;
}

namespace atlas {
namespace runtime {
namespace trace {

/// @class CallStack
/// Instances of CallStack can keep track of nested CodeLocations
class CallStack {
public:
    using const_iterator = std::vector<size_t>::const_iterator;

public:
    void push(const CodeLocation&, const std::string& id = "");
    void pop();

    const_iterator begin() const { return stack_.begin(); }
    const_iterator end() const { return stack_.begin() + size_; }

    size_t hash() const;
    size_t size() const { return size_; }

    operator bool() const { return size_ > 0; }

public:
    CallStack(): stack_(64){};
    CallStack(const CallStack& other): stack_(other.stack_), size_(other.size_) {}
    CallStack& operator=(const CallStack& other) {
        stack_ = other.stack_;
        size_  = other.size_;
        hash_  = 0;
        return *this;
    }

private:
    std::vector<size_t> stack_;
    size_t size_{0};
    mutable size_t hash_{0};
};

}  // namespace trace
}  // namespace runtime
}  // namespace atlas
