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

#include <iosfwd>

//-----------------------------------------------------------------------------------------------------------

namespace atlas {
namespace runtime {
namespace trace {

class Control {
public:
    static bool enabled();
};

//-----------------------------------------------------------------------------------------------------------

// Class used to avoid any printing before and after a timer
class NoLogging {
public:
    NoLogging(bool state);
    NoLogging(std::ostream& channel);

public:  // static methods
    static std::ostream& channel();
    static bool enabled();
    static void start(const std::string&) {}
    static void stop(const std::string&, double) {}
};

//-----------------------------------------------------------------------------------------------------------

// Class used to print message before and after a timer
// Example print:
//     timer-name ...
//     timer-name ... done : 5s
class Logging {
public:
    Logging(bool state);
    Logging(std::ostream& channel);
    virtual ~Logging();

public:  // static methods
    static std::ostream& channel();
    static bool enabled();
    static void start(const std::string& title);
    static void stop(const std::string& title, double seconds);

private:
    std::ostream* previous_state_;
};

//-----------------------------------------------------------------------------------------------------------

// Class used to print message only upon end of a timer
// Example print:
//     timer-name : 5s
class LoggingResult : public Logging {
public:
    using Logging::Logging;

public:  // static methods
    static void start(const std::string&) {}
    static void stop(const std::string& title, double seconds);
};

//-----------------------------------------------------------------------------------------------------------

}  // namespace trace
}  // namespace runtime
}  // namespace atlas
