/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#pragma once

#include <string>
#include <vector>

//-----------------------------------------------------------------------------------------------------------

namespace eckit {
class Configuration;
}
namespace atlas {
class CodeLocation;
}

namespace atlas {
namespace runtime {
namespace trace {

class CallStack;

class Timings {
public:
    using Configuration = eckit::Configuration;
    using CodeLocation  = atlas::CodeLocation;
    using Identifier    = size_t;
    using Labels        = std::vector<std::string>;

public:  // static methods
    static Identifier add(const CodeLocation&, const CallStack&, const std::string& title, const Labels&);

    static void update(const Identifier& id, double seconds);

    static std::string report();

    static std::string report(const Configuration&);
};

}  // namespace trace
}  // namespace runtime
}  // namespace atlas
