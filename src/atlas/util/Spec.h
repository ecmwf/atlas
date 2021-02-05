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

#include "atlas/runtime/Exception.h"
#include "atlas/runtime/Log.h"
#include "atlas/util/Config.h"

namespace eckit {
class PathName;
}

namespace atlas {
namespace util {

//---------------------------------------------------------------------------------------------------------------------

using Spec = Config;

//---------------------------------------------------------------------------------------------------------------------

template <typename T>
struct SpecFactory {
    static Spec create( const std::string& id ) {
        const auto& m = instance().m_;

        auto j = m.find( id );
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

    struct Register {
        Register( const std::string& id, const eckit::PathName& path ) : Register( id, Spec( path ) ) {}
        Register( const std::string& id, const Spec& spec ) {
            ATLAS_ASSERT_MSG( instance().m_.emplace( id, spec ).second, "SpecFactory: duplicate '" + id + "'" );
        }
    };

private:
    SpecFactory() = default;

    static SpecFactory& instance() {
        static SpecFactory factory;
        return factory;
    }

    SpecFactory( const SpecFactory& ) = delete;
    SpecFactory& operator=( const SpecFactory& ) = delete;

    std::map<std::string, const Spec> m_;
};

}  // namespace util
}  // namespace atlas
