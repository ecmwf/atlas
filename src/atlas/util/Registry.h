/*
 * (C) Copyright 2021 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#pragma once

#include <map>
#include <ostream>
#include <string>
#include <utility>
#include <vector>

#include "atlas/runtime/Exception.h"
#include "atlas/runtime/Log.h"

namespace atlas {
namespace util {

template <typename T>
class Registry {
public:
    Registry() = default;

    void add(const std::string& key, T&& value) {
        ATLAS_ASSERT_MSG(m_.emplace(key, std::move(value)).second, "ObjectRegistry: duplicate '" + key + "'");
    }

    bool has(const std::string& key) { return m_.find(key) != m_.end(); }

    T get(const std::string& key) {
        auto j = m_.find(key);
        if (j != m_.end()) {
            return j->second;
        }

        list(Log::error() << "Registry: unknown '" << key << "', choices are: ");
        throw_Exception("Registry: unknown '" + key + "'");
    }

    void list(std::ostream& out) {
        const char* sep = "";
        for (const auto& j : m_) {
            out << sep << j.first;
            sep = ", ";
        }
    }

    std::vector<std::string> keys() {
        std::vector<std::string> keys;
        for (const auto& j : m_) {
            keys.push_back(j.first);
        }
        return keys;
    }

private:
    std::map<std::string, T> m_;
};

}  // namespace util
}  // namespace atlas
