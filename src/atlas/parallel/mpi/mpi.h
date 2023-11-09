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

#include <string_view>

#include "eckit/mpi/Comm.h"
#include "atlas/parallel/mpi/Statistics.h"

namespace atlas {
namespace mpi {

using Comm = eckit::mpi::Comm;

inline const Comm& comm() {
    return eckit::mpi::comm();
}

inline const Comm& comm(std::string_view name) {
    return eckit::mpi::comm(name.data());
}

inline idx_t rank() {
    return static_cast<idx_t>(comm().rank());
}

inline int size() {
    return static_cast<idx_t>(comm().size());
}

void finalize();
void finalise();

class CommStack {
public:
    using const_iterator = std::vector<std::string>::const_iterator;

public:
    void push(std::string_view name);
    void pop();
    void pop(std::string_view name); // verifies name matches compared to version above
    const std::string& name() const;
    const mpi::Comm& comm() const;

    const_iterator begin() const { return stack_.begin(); }
    const_iterator end() const { return stack_.begin() + size_; }

    size_t size() const { return size_; }

    static CommStack& instance() {
        static CommStack instance;
        return instance;
    }

private:
    CommStack();
private:
    std::vector<std::string> stack_;
    size_t size_{0};
};
void push(std::string_view name);
void pop(std::string_view name);
void pop();
struct Scope {
    Scope(std::string_view name) : name_(name) {
        push(name_);
    }
    ~Scope() {
        pop(name_);
    }
    std::string name_;
};

}  // namespace mpi
}  // namespace atlas
