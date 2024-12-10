/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <iostream>

#include "atlas/runtime/Exception.h"
#include "atlas/runtime/Log.h"
#include "atlas/util/Factory.h"

// #define DEBUG_FACTORY_REGISTRATION

using lock_guard = std::lock_guard<std::mutex>;

namespace atlas {
namespace util {

bool FactoryRegistry::has(const std::string& builder) const {
    lock_guard lock(mutex_);
    return (factories_.find(builder) != factories_.end());
}

FactoryBase* FactoryRegistry::get(const std::string& builder) const {
    lock_guard lock(mutex_);
    auto iterator = factories_.find(builder);

    if (iterator == factories_.end()) {
        Log::error() << "No " << factory_ << " for [" << builder << "]" << std::endl;
        Log::error() << "Factories are:" << std::endl;
        for (const auto& map_pair : factories_) {
            Log::error() << "   " << map_pair.first << std::endl;
        }
        throw_Exception(std::string("No ") + factory_ + std::string(" called ") + builder);
    }
    else {
        return iterator->second;
    }
}

void FactoryRegistry::add(const std::string& builder, FactoryBase* factory) {
    lock_guard lock(mutex_);
    ATLAS_ASSERT(factories_.find(builder) == factories_.end(), "Cannot find builder in factory");
    factories_[builder] = factory;
#ifdef DEBUG_FACTORY_REGISTRATION
    std::cout << "Registered " << builder << " (" << factory << ") in " << factory_ << std::endl << std::flush;
#endif
}

void FactoryRegistry::remove(const std::string& builder) {
    lock_guard lock(mutex_);
    ATLAS_ASSERT(factories_.find(builder) != factories_.end());
#ifdef DEBUG_FACTORY_REGISTRATION
    FactoryBase* factory = factories_[builder];
    std::cout << "Unregistered " << builder << " (" << factory << ") from " << factory_ << std::endl << std::flush;
#endif
    factories_.erase(builder);
}

FactoryRegistry::FactoryRegistry(const std::string& factory): factory_(factory) {
#ifdef DEBUG_FACTORY_REGISTRATION
    std::cout << "Created " << factory << std::endl;
#endif
}

FactoryRegistry::~FactoryRegistry() {
    ATLAS_ASSERT(factories_.empty(), "Registry should only be destroyed once all builders have been unregistrered");
#ifdef DEBUG_FACTORY_REGISTRATION
    std::cout << "Destroyed " << factory_ << std::endl << std::flush;
#endif
}

std::vector<std::string> FactoryRegistry::keys() const {
    lock_guard lock(mutex_);
    std::vector<std::string> _keys;
    _keys.reserve(factories_.size());
    for (const auto& key_value : factories_) {
        _keys.emplace_back(key_value.first);
    }
    return _keys;
}

void FactoryRegistry::list(std::ostream& out) const {
    lock_guard lock(mutex_);
    const char* sep = "";
    for (const auto& map_pair : factories_) {
        out << sep << map_pair.first;
        sep = ", ";
    }
}


//----------------------------------------------------------------------------------------------------------------------

FactoryBase::FactoryBase(FactoryRegistry& registry, const std::string& builder):
    registry_(registry), builder_(builder) {
    if (not builder_.empty()) {
        registry_.add(builder, this);
    }
}

FactoryBase::~FactoryBase() {
    if (not builder_.empty()) {
        registry_.remove(builder_);
    }
}

//----------------------------------------------------------------------------------------------------------------------

std::shared_ptr<FactoryRegistry> FactoryRegistry::instance(const std::string& factory) {
    static std::map<std::string,std::shared_ptr<FactoryRegistry>> registries;
    if (registries.find(factory) == registries.end()) {
        auto [it_pair, inserted] = registries.emplace(factory, new FactoryRegistry(factory));
        return it_pair->second;
    }
    return registries.at(factory);
}

}  // namespace util
}  // namespace atlas
