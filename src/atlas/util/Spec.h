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

#include <iostream>
#include <map>
#include <string>
#include <utility>

#include "atlas/runtime/Exception.h"
#include "atlas/runtime/Log.h"
#include "atlas/util/Config.h"

namespace eckit {
class PathName;
}

namespace atlas {
namespace util {

using Spec = Config;

template <typename T>
struct SpecRegistry {
    struct Enregister {
        Enregister( const std::string& id, const eckit::PathName& path ) { enregister( id, path ); }
        Enregister( const std::string& id, Spec&& spec ) { enregister( id, std::move( spec ) ); }
    };

    static void enregister( const std::string& id, const eckit::PathName& path ) {
        ATLAS_ASSERT_MSG( instance().m_.find( id ) == instance().m_.end(), duplicate( id ) );
        ATLAS_ASSERT_MSG( instance().p_.emplace( id, path ).second, duplicate( id ) );
    }

    static void enregister( const std::string& id, Spec&& spec ) {
        ATLAS_ASSERT_MSG( instance().m_.emplace( id, spec ).second, duplicate( id ) );
    }

    static Spec lookup( const std::string& id ) {
        auto& p = instance().p_;
        auto i  = p.find( id );
        if ( i != p.end() ) {
            enregister( id, std::move( Spec( i->second ) ) );
            p.erase( i );
            return lookup( id );
        }

        auto& m = instance().m_;
        auto j  = m.find( id );
        if ( j != m.end() ) {
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

    static std::string duplicate( const std::string& id ) { return "SpecFactory: duplicate '" + id + "'"; }

    std::map<std::string, const Spec> m_;
    std::map<std::string, const eckit::PathName> p_;
};

}  // namespace util
}  // namespace atlas
