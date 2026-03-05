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

namespace atlas::mpi {

using Comm = eckit::mpi::Comm;

inline const Comm& comm() {
    return eckit::mpi::comm();
}

inline const Comm& comm(std::string_view name) {
    return eckit::mpi::comm(name.data());
}

const Comm& comm(int communicator);

bool has_comm(std::string_view name);

const Comm& register_comm(std::string_view name, int communicator);
void unregister_comm(std::string_view name);

inline idx_t rank() {
    return static_cast<idx_t>(comm().rank());
}

inline int size() {
    return static_cast<idx_t>(comm().size());
}

void finalize();
void finalise();

namespace scope {
    /// @brief Push a new scope for MPI communicators. The current default communicator is restored when the scope is destructed using pop()
    void push();
    /// @brief Push a new scope for MPI communicators. The given communicator is set as the default communicator for the duration of the scope. The previous default communicator is restored when the scope is destructed using pop()
    void push(std::string_view name);
    /// @brief Push a new scope for MPI communicators. The given communicator is set as the default communicator for the duration of the scope. The previous default communicator is restored when the scope is destructed using pop()
    void push(const Comm& comm);
    /// @brief Push a new scope for MPI communicators. The given communicator is set as the default communicator for the duration of the scope. The previous default communicator is restored when the scope is destructed using pop()
    /// If the given communicator is not already registered, it will be registered with a generated name "int.<communicator>" for the duration of the scope, and unregistered when the scope is destructed.
    void push(int communicator);
    /// @brief Pop the current MPI communicator scope, restoring the previous default communicator
    void pop();
}

/// @brief RAII helper class to manage MPI communicator scopes. The constructor pushes a new scope, and the destructor pops the scope, ensuring that the previous default communicator is restored even in case of exceptions.
struct Scope {
    Scope() { scope::push(); }
    Scope(std::string_view name) { scope::push(name); }
    Scope(const Comm& comm) { scope::push(comm); };
    /// @brief Constructor using integer value. If the given communicator is not already registered, it will be registered with a generated name "int.<communicator>" for the duration of the scope, and unregistered when the scope is destructed.
    /// @param communicator The integer value of the MPI communicator to use for the scope
    Scope(int communicator) { scope::push(communicator); }

    Scope(const Scope&) = delete;
    Scope(Scope&&) = delete;
    Scope& operator=(const Scope&) = delete;
    Scope& operator=(Scope&&) = delete;
    ~Scope() { scope::pop(); }
};

[[deprecated("Use atlas::mpi::scope::push instead")]] inline void push() { scope::push(); }
[[deprecated("Use atlas::mpi::scope::pop instead")]]  inline void pop() { scope::pop(); }
[[deprecated("Use atlas::mpi::scope::push instead")]] inline void push(std::string_view name) { scope::push(name); }
[[deprecated("Use atlas::mpi::scope::pop instead")]]  inline void pop(std::string_view /*name*/) { scope::pop(); }

}  // namespace atlas::mpi
