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

//-----------------------------------------------------------------------------------------------------------

namespace atlas {
namespace runtime {
namespace trace {

class NoBarriers {
public:
    NoBarriers(bool state) {}
    void restore() {}

public:  // static methods
    static bool state() { return false; }
    static void execute() {}
    static double time();
    static std::string report();
};

class Barriers {
public:
    Barriers(bool state);
    ~Barriers();
    void restore();

public:  // static methods
    static bool state();
    static void execute();
    static double time();
    static std::string report();

private:
    bool previous_state_;
};

}  // namespace trace
}  // namespace runtime
}  // namespace atlas
