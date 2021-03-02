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
#include <list>
#include <string>

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
    using const_iterator         = std::list<size_t>::const_iterator;
    using const_reverse_iterator = std::list<size_t>::const_reverse_iterator;

public:
    void push_front( const CodeLocation&, const std::string& id = "" );
    void pop_front();

    const_iterator begin() const { return stack_.begin(); }
    const_iterator end() const { return stack_.end(); }

    const_reverse_iterator rbegin() const { return stack_.rbegin(); }
    const_reverse_iterator rend() const { return stack_.rend(); }

    size_t hash() const;
    size_t size() const { return stack_.size(); }

    operator bool() const { return not stack_.empty(); }

public:
    CallStack() = default;
    CallStack( const CallStack& other ) : stack_( other.stack_ ) {}
    CallStack& operator=( const CallStack& other ) {
        stack_ = other.stack_;
        hash_  = 0;
        return *this;
    }

private:
    std::list<size_t> stack_;
    mutable size_t hash_{0};
};

}  // namespace trace
}  // namespace runtime
}  // namespace atlas
