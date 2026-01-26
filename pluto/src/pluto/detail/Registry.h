/*
 * (C) Copyright 2024- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */
#pragma once

#include <algorithm>
#include <exception>
#include <iostream>
#include <list>
#include <map>
#include <mutex>
#include <string>
#include <string_view>

#include "pluto/trace.h"

namespace pluto {

template <class T>
class Registry {
public:
    using value_type = T;

    static Registry& instance() {
        static Registry instance;
        return instance;
    }

    value_type& enregister(std::string_view name, std::unique_ptr<value_type>&& mr) {
        std::lock_guard lock(mutex_);
        auto& ref = do_register(name, *mr);

        owned_.emplace(name, std::move(mr));
        ordered_keys_.push_back(std::string(name));
        return ref;
    }

    value_type& enregister(std::string_view name, value_type& mr) {
        std::lock_guard lock(mutex_);
        return do_register(name, mr);
    }

    void unregister(std::string_view name) {
        std::lock_guard lock(mutex_);
        std::string key(name);
        if (registered_.erase(key) == 0) {
            throw std::runtime_error("Could not unregister " + key);
        }
        if (original_name_.erase(key) == 0) {
            throw std::runtime_error("Could not unregister " + key);
        }
        if (owned_.erase(key)) {
            ordered_keys_.erase(std::find(ordered_keys_.begin(), ordered_keys_.end(), key));
        }
        if (trace::enabled()) {
            trace::out << "unregistered " << name << std::endl;
        }
    }

    bool has(std::string_view name) const {
        std::lock_guard lock(mutex_);
        return do_has(name);
    }

    value_type& get(std::string_view name) const {
        std::lock_guard guard(mutex_);
        return do_get(name);
    }

    void clear() {
        std::lock_guard guard(mutex_);
        return do_clear();
    }

    std::string_view name(void* mr) const {
        std::lock_guard guard(mutex_);
        return do_name(mr);
    }

private:
    Registry() = default;

    ~Registry() { do_clear(); }

    void do_clear() {
        // Erase owned entries in reverse order of insertion
        // No need to unregister as this is in program teardown
        for (auto it = ordered_keys_.rbegin(); it != ordered_keys_.rend(); ++it) {
            auto& key = *it;
            if (trace::enabled()) {
                trace::out << "~Registry() : Deleting owned " << key << std::endl;
            }
            owned_.erase(key);
        }
        registered_.clear();
        original_name_.clear();
        owned_.clear();
        ordered_keys_.clear();
    }

    value_type& do_register(std::string_view name, value_type& mr) {
        if (do_has(mr)) {
            original_name_[std::string(name)] = std::string(do_name(&mr));
        }
        else {
            original_name_[std::string(name)] = std::string(name);
        }
        bool inserted = registered_.try_emplace(std::string(name), &mr).second;
        if (not inserted) {
            throw std::runtime_error("Could not register " + std::string(name));
        }
        if (trace::enabled()) {
            trace::out << "registered " << name << std::endl;
        }
        return mr;
    }

    bool do_has(std::string_view name) const {
        if (registered_.find(name) == registered_.end()) {
            return false;
        }
        return true;
    }

    bool do_has(value_type& mr) const {
        for (auto& key_val : registered_) {
            if (key_val.second == &mr) {
                return true;
            }
        }
        return false;
    }

    std::string_view do_name(void *mr) const {
        std::string_view name;
        for (auto& [key, val] : registered_) {
            if (val == mr) {
                name = original_name_.at(std::string(key));
                break;
            }
        }
        return name;
    }

    value_type& do_get(std::string_view name) const {
        if (registered_.find(name) != registered_.end()) {
            return *registered_.at(std::string(name));
        }
        throw std::runtime_error("Could not get registered " + std::string(name));
    }

private:
    mutable std::mutex mutex_;
    std::map<std::string, std::unique_ptr<value_type>, std::less<>> owned_;
    std::list<std::string> ordered_keys_;
    std::map<std::string, value_type*, std::less<>> registered_;
    std::map<std::string, std::string> original_name_;
};

}  // namespace pluto
