/*
 * (C) Copyright 1996-2017 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#pragma once

#include <iosfwd>

//-----------------------------------------------------------------------------------------------------------

namespace atlas {
namespace runtime {
namespace timer {

//-----------------------------------------------------------------------------------------------------------

// Class used to avoid any printing before and after a timer
class TimerTracingNone {
public:
    TimerTracingNone( bool state );
    TimerTracingNone( std::ostream& channel );

public: // static methods
    static std::ostream& channel();
    static bool enabled();
    static void start( const std::string& ) {}
    static void stop( const std::string&, double ) {}
};

//-----------------------------------------------------------------------------------------------------------

// Class used to print message before and after a timer
// Example print:
//     timer-name ...
//     timer-name ... done : 5s
class TimerTracing {
public:
    TimerTracing( bool state );
    TimerTracing( std::ostream& channel );
    virtual ~TimerTracing();

public: // static methods
    static std::ostream& channel();
    static bool enabled();
    static void start( const std::string& title );
    static void stop( const std::string& title, double seconds );
    
private:
  std::ostream* previous_state_;
};

//-----------------------------------------------------------------------------------------------------------

// Class used to print message only upon end of a timer
// Example print:
//     timer-name : 5s
class TimerTracingResult : public TimerTracing {
public:
    using TimerTracing::TimerTracing;

public: // static methods
    static void start( const std::string& ) {}
    static void stop( const std::string& title, double seconds );
};

//-----------------------------------------------------------------------------------------------------------

} // namespace timer
} // namespace runtime
} // namespace atlas

