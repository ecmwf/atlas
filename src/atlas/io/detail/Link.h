/*
 * (C) Copyright 2020 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#pragma once

#include <string>

namespace eckit {
class PathName;
}

namespace atlas {
namespace io {

struct Link {
    std::string uri;

    const std::string& str() const { return uri; }
    operator bool() const { return uri.size(); }
    operator const std::string&() const { return str(); }

    //    bool relative() const;

    //    friend Link operator/( const eckit::PathName& dir, const Link& link );
};


}  // namespace io
}  // namespace atlas
