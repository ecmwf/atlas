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

#include <ostream>
#include <map>
#include <string>
#include <utility>

#include "atlas/runtime/Exception.h"
#include "atlas/runtime/Log.h"
#include "atlas/util/Config.h"

namespace atlas {
namespace util {

using Spec = Config;

template <typename T>
struct SpecRegistry {
    struct Enregister {
        Enregister( const std::string& id, Spec&& spec ) { enregister( id, std::move( spec ) ); }
    };

    static void enregister( const std::string& id, Spec&& spec ) {
        ATLAS_ASSERT_MSG( instance().m_.emplace( id, std::move( spec ) ).second,
                          "SpecFactory: duplicate '" + id + "'" );
    }

    static bool has( const std::string& id ) { return instance().m_.find( id ) != instance().m_.end(); }

    static Spec lookup( const std::string& id ) {
        auto j = instance().m_.find( id );
        if ( j != instance().m_.end() ) {
            return j->second;
        }

        list( Log::error() << "SpecFactory: unknown '" << id << "', choices are: " );
        throw_Exception( "SpecFactory: unknown '" + id + "'" );
    }

    static void list( std::ostream& out ) {
        const char* sep = "";
        for ( const auto& j : instance().m_ ) {
            out << sep << j.first;
            sep = ", ";
        }
    }

private:
    static SpecRegistry& instance() {
        static SpecRegistry factory;
        return factory;
    }

    std::map<std::string, const Spec> m_;
};

}  // namespace util
}  // namespace atlas
