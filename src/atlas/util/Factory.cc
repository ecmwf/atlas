/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "atlas/util/Factory.h"

#include <iostream>
#include <cstdlib>

#include "atlas/runtime/Exception.h"
#include "atlas/runtime/Log.h"

// #define DEBUG_FACTORY_REGISTRATION

using lock_guard = std::lock_guard<std::mutex>;

namespace atlas {
namespace util {

static bool ATLAS_DEPRECATION_WARNINGS() {
    const char* val = std::getenv("ATLAS_DEPRECATION_WARNINGS");
    if (val != nullptr) {
        return std::atoi(val);
    }
    return true;
}

static bool ATLAS_DEPRECATION_ERRORS() {
    const char* val = std::getenv("ATLAS_DEPRECATION_ERRORS");
    if (val != nullptr) {
        return std::atoi(val);
    }
    return false;
}

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
        auto* factory = iterator->second;
        if (factory->deprecated()) {
            if (ATLAS_DEPRECATION_WARNINGS()) {
                const std::string& message = factory->deprecated().message();
                Log::warning() << "[ATLAS_DEPRECATION_WARNING] The builder " << builder << " should no longer be used. " << message << '\n';
                Log::warning() << "[ATLAS_DEPRECATION_WARNING] This warning can be disabled with `export ATLAS_DEPRECATION_WARNINGS=0`" << std::endl;
            }
            if (ATLAS_DEPRECATION_ERRORS()) {
                const std::string& message = factory->deprecated().message();
                ATLAS_THROW_EXCEPTION("[ATLAS_DEPRECATION_ERROR] The builder " << builder << " should no longer be used. " << message);
            }
        }
        return factory;
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

FactoryRegistry::FactoryRegistry(const std::string& factory, FactoryRegistry::Private): factory_(factory) {
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
        if (not key_value.second->deprecated_) {
            _keys.emplace_back(key_value.first);
        }
    }
    return _keys;
}

void FactoryRegistry::list(std::ostream& out) const {
    lock_guard lock(mutex_);
    const char* sep = "";
    for (const auto& key_value : factories_) {
        if (key_value.second->deprecated_) {
            out << sep << key_value.first;
            sep = ", ";
        }
    }
}


//----------------------------------------------------------------------------------------------------------------------

FactoryBase::FactoryBase(FactoryRegistry& registry, const std::string& builder, const FactoryDeprecated& deprecated):
    registry_(registry), builder_(builder), deprecated_(deprecated) {
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
        auto [it_pair, inserted] = registries.emplace(factory, new FactoryRegistry(factory,Private()));
        return it_pair->second;
    }
    return registries.at(factory);
}

}  // namespace util
}  // namespace atlas
