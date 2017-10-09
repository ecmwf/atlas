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

class TimerNoLogging {
public:

    TimerNoLogging( bool state );

    TimerNoLogging( std::ostream& channel );

public: // static method
    static std::ostream& channel();

};

class TimerLogging {
public:

    TimerLogging( bool state );

    TimerLogging( std::ostream& channel );
    
    ~TimerLogging();

public: // static method
    static std::ostream& channel();
    
private:
  std::ostream* previous_state_;

};

} // namespace timer
} // namespace runtime
} // namespace atlas

