/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "Backend.h"

#include <map>

#include "atlas/runtime/Exception.h"
#include "atlas/util/Config.h"

namespace atlas {
namespace linalg {
namespace sparse {

namespace {
struct backends {
    std::map<std::string, sparse::Backend> map_;
    std::string current_backend_{backend::omp::type()};

    static backends& instance() {
        static backends x;
        return x;
    }

    void set( const std::string& current_backend ) { current_backend_ = current_backend; }

    sparse::Backend& get( const std::string& type ) {
        if ( map_.find( type ) == map_.end() ) {
            map_.emplace( type, util::Config( "type", type ) );
        }
        return map_[type];
    }

    sparse::Backend& get() { return get( current_backend_ ); }

private:
    backends() = default;
};
}  // namespace

void current_backend( const std::string& backend ) {
    backends::instance().set( backend );
}
sparse::Backend& current_backend() {
    return backends::instance().get();
}

sparse::Backend& default_backend( const std::string& backend ) {
    return backends::instance().get( backend );
}

Backend::Backend( const eckit::Configuration& other ) : util::Config( other ) {
    ATLAS_ASSERT( has( "type" ) );
}

}  // namespace sparse
}  // namespace linalg
}  // namespace atlas
