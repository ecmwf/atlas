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

#include <chrono>

namespace atlas {
namespace runtime {
namespace timer {
  
//-----------------------------------------------------------------------------------------------------------
  
class StopWatch {
private:
  std::chrono::duration<double> elapsed_{0};
  std::chrono::steady_clock::time_point start_;
public:
  void start() {
    start_ = std::chrono::steady_clock::now();
  }
  void stop() { 
    elapsed_ += (std::chrono::steady_clock::now() - start_);
  }
  void reset() {
    elapsed_ = std::chrono::seconds::zero();
    start_ = std::chrono::steady_clock::now();
  }
  double elapsed() const {
    return elapsed_.count();
  }
};

//-----------------------------------------------------------------------------------------------------------

} // timer
} // runtime
} // atlas
