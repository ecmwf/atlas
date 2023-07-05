/*
 * (C) Copyright 2023- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <string>

#include "atlas-example-plugin/Library.h"
#include "atlas-example-plugin/version.h"


namespace atlas {
namespace example_plugin {


REGISTER_LIBRARY( Library );


Library::Library() : Plugin( "atlas-example-plugin" ) {}


const Library& Library::instance() {
    static Library library;
    return library;
}


std::string Library::version() const {
    return atlas_example_plugin_version();
}


std::string Library::gitsha1( unsigned int count ) const {
    std::string sha1 = atlas_example_plugin_git_sha1();
    return sha1.empty() ? "not available" : sha1.substr( 0, std::min( count, 40U ) );
}


void Library::init() {
    Plugin::init();
}


}  // namespace example_plugin
}  // namespace atlas
