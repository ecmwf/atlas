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
#include "atlas/runtime/Log.h"
#include "atlas/runtime/Exception.h"

namespace atlas {
namespace mpi {

void finalize() {
    finalise();
}
void finalise() {
    Log::debug() << "atlas::mpi::finalize() --> Finalizing MPI" << std::endl;
    eckit::mpi::finaliseAllComms();
}

void CommStack::push(std::string_view name) {
    if (stack_.size() == size_) {
        stack_.resize(2 * size_);
    }
    stack_[size_++] = name;
    eckit::mpi::setCommDefault(name.data());
}

void CommStack::pop(std::string_view _name) {
    ATLAS_ASSERT(_name == name());
    pop();
}

void CommStack::pop() {
    --size_;
    eckit::mpi::setCommDefault(name().c_str());
}

const std::string& CommStack::name() const {
    return stack_[size_-1];
}

const mpi::Comm& CommStack::comm() const {
    return mpi::comm(name());
}

CommStack::CommStack(): stack_(64){
        stack_[size_++] = mpi::comm().name();
    };

void push(std::string_view name) {
    Log::debug() << "atlas::mpi::push("<<name<<")" << std::endl;
    CommStack::instance().push(name);
}

void pop(std::string_view name) {
    Log::debug() << "atlas::mpi::pop("<<mpi::comm().name()<<")" << std::endl;
    CommStack::instance().pop(name);
}

void pop() {
    Log::debug() << "atlas::mpi::pop("<<mpi::comm().name()<<")" << std::endl;
    CommStack::instance().pop();
}

}  // namespace mpi
}  // namespace atlas
