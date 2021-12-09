/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <cstring>
#include <fstream>
#include <sstream>

#include "eckit/filesystem/PathName.h"
#include "eckit/runtime/Main.h"

#include "atlas/runtime/Exception.h"
#include "atlas_f/internals/atlas_read_file.h"

namespace atlas {

void read_file(const eckit::PathName& p, std::ostream& out) {
    if (p.exists()) {
        std::ifstream in;
        in.open(p.asString().c_str());
        if (!in) {
            throw_CantOpenFile(p.asString(), Here());
        }
        else {
            out << in.rdbuf();
            in.close();
        }
    }
}

}  // namespace atlas

extern "C" {

//-----------------------------------------------------------------------------

int atlas__read_file(const char* path, char*& content, int& size) {
    // eckit::FileReadPolicy p =
    // eckit::Main::instance().behavior().fileReadPolicy();

    // std::stringstream ss;

    // if( read( p, path, ss ) )
    // {
    //   std::string s = ss.str();
    //   size = s.size()+1;
    //   content = new char[size];
    //   strcpy(content,s.c_str());
    //   return true;
    // }

    std::stringstream ss;
    atlas::read_file(path, ss);
    std::string s = ss.str();
    size          = static_cast<int>(s.size());
    content       = new char[size + 1];
    std::strncpy(content, s.c_str(), size + 1);
    return true;
}

//-----------------------------------------------------------------------------
}
