/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "atlas/parallel/mpi/mpi.h"

#include <functional>
#include <memory>
#include <optional>
#include <sstream>
#include <stack>
#include <vector>

#include "atlas/runtime/Log.h"
#include "atlas/runtime/Exception.h"

namespace atlas::mpi {

void finalize() {
    finalise();
}
void finalise() {
    Log::debug() << "atlas::mpi::finalize() --> Finalizing MPI" << std::endl;
    eckit::mpi::finaliseAllComms();
}

namespace {
const Comm& comm_world() {
    static const Comm& world = atlas::mpi::comm("world");
    return world;
}
const Comm& comm_self() {
    static const Comm& self = atlas::mpi::comm("self");
    return self;
}
int comm_world_int() {
    static int world = comm_world().communicator();
    return world;
}
int comm_self_int() {
    static int self = comm_self().communicator();
    return self;
}

std::string comm_name_from_communicator( int communicator ) {
    return "int." + std::to_string(communicator);
}

std::optional<std::reference_wrapper<const Comm>> lookup_comm( int communicator ) {
    // First check if the communicator corresponds to the well-known comm_world or comm_self, or the current default comm
    // Then fall back to a more expensive lookup in the registered communicators, if any
    if ( communicator == comm_world_int() ) {
        return comm_world();
    }
    else if ( communicator == comm_self_int() ) {
        return comm_self();
    }
    else {
        const auto& default_comm = comm();
        if ( default_comm.communicator() == communicator ) {
            return default_comm;
        }
        else {
            for( std::string name : eckit::mpi::listComms() ) {
                const auto& comm = atlas::mpi::comm(name);
                if ( comm.communicator() == communicator ) {
                    return comm;
                }
            }
        }
    }
    return std::nullopt;
}

class CommScope {
public:
    CommScope() {
        Log::debug() << "atlas::mpi::scope::push()" << std::endl;
        previous_comm_name_ = mpi::comm().name();
    }
    CommScope(std::string_view name) {
        Log::debug() << "atlas::mpi::scope::push(" << name << ")" << std::endl;
        previous_comm_name_ = mpi::comm().name();
        eckit::mpi::setCommDefault(name.data());
    }
    CommScope(const Comm& comm) {
        Log::debug() << "atlas::mpi::scope::push(" << comm.name() << ")" << std::endl;
        previous_comm_name_ = mpi::comm().name();
        eckit::mpi::setCommDefault(comm.name().data());
    }
    CommScope(int communicator) {
        previous_comm_name_ = mpi::comm().name();
        auto comm_opt = lookup_comm(communicator);
        if (comm_opt) {
            Log::debug() << "atlas::mpi::scope::push(" << communicator << ") : found existing comm " << comm_opt->get().name() << std::endl;
            eckit::mpi::setCommDefault(comm_opt->get().name().data());
        }
        else {
            new_registered_comm_name_ = comm_name_from_communicator(communicator);
            Log::debug() << "atlas::mpi::scope::push(" << communicator << ") : registered new comm " << new_registered_comm_name_ << std::endl;
            register_comm(new_registered_comm_name_, communicator);
            eckit::mpi::setCommDefault(new_registered_comm_name_.c_str());
        }
    }
    ~CommScope() {
        Log::debug() << "atlas::mpi::scope::pop() : restored to " << previous_comm_name_;
        eckit::mpi::setCommDefault(previous_comm_name_.c_str());
        if (not new_registered_comm_name_.empty()) {
            Log::debug() << " and unregistered comm " << new_registered_comm_name_;
            unregister_comm(new_registered_comm_name_);
        }
        Log::debug() << '\n';
    }
private:
    std::string previous_comm_name_;
    std::string new_registered_comm_name_;
};

class ScopeStack {
    using value_type = std::unique_ptr<CommScope>;
    using container_type = std::vector<value_type>;
    using stack_type = std::stack<value_type, container_type>;
public:
    void push() {
        stack_.push(std::make_unique<CommScope>());
    }
    void push(std::string_view name) {
        stack_.push(std::make_unique<CommScope>(name));
    }
    void push(const Comm& comm) {
        stack_.push(std::make_unique<CommScope>(comm));
    }
    void push(int communicator) {
        stack_.push(std::make_unique<CommScope>(communicator));
    }
    void pop() {
        if( !stack_.empty() ) {
            stack_.pop();
        }
    }
    static ScopeStack& instance() {
        static ScopeStack instance;
        return instance;
    }

private:
    ScopeStack() {
        container_type underlying;
        underlying.reserve(64);
        stack_ = stack_type(std::move(underlying));
    }
private:
    stack_type stack_;
};
}

void scope::push() {
    ScopeStack::instance().push();
}

void scope::push(std::string_view name) {
    ScopeStack::instance().push(name);
}

void scope::push(const Comm& comm) {
    ScopeStack::instance().push(comm);
}

void scope::push(int communicator) {
    ScopeStack::instance().push(communicator);
}

void scope::pop() {
    ScopeStack::instance().pop();
}

bool has_comm( std::string_view name ) {
    return eckit::mpi::hasComm(name.data());
}

const Comm& comm(int communicator) {
    auto comm_opt = lookup_comm( communicator );
    if ( comm_opt ) {
        return comm_opt->get();
    }
    throw_Exception("No MPI communicator found with integer value " + std::to_string(communicator), Here());
}

const Comm& register_comm(std::string_view name, int communicator) {
    eckit::mpi::addComm(name.data(), communicator);
    return comm(name);
}

void unregister_comm(std::string_view name) {
#if ATLAS_ECKIT_VERSION_AT_LEAST(2, 0, 0)
    eckit::mpi::unregisterComm(name.data());
#else
    Log::warning() << "atlas::mpi::unregister_comm is a no-op with eckit versions prior to 2.0.0" << std::endl;
#endif
}


}  // namespace atlas::mpi
