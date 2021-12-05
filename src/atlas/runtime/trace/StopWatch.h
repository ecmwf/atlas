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

#include <chrono>

namespace atlas {
namespace runtime {
namespace trace {

//-----------------------------------------------------------------------------------------------------------

class StopWatch {
public:
    StopWatch();
    StopWatch(double elapsed);
    void start();
    void stop();
    void reset();
    double elapsed() const;

private:
    std::chrono::duration<double> elapsed_;
    std::chrono::steady_clock::time_point start_;
    bool running_{false};
};

//-----------------------------------------------------------------------------------------------------------

inline StopWatch::StopWatch(): elapsed_(0) {}

inline StopWatch::StopWatch(double elapsed): elapsed_(elapsed) {}

inline void StopWatch::start() {
    start_   = std::chrono::steady_clock::now();
    running_ = true;
}

inline void StopWatch::stop() {
    if (running_) {
        elapsed_ += (std::chrono::steady_clock::now() - start_);
        running_ = false;
    }
}

inline void StopWatch::reset() {
    elapsed_ = std::chrono::seconds::zero();
    start_   = std::chrono::steady_clock::now();
}

inline double StopWatch::elapsed() const {
    return elapsed_.count();
}

//-----------------------------------------------------------------------------------------------------------

}  // namespace trace
}  // namespace runtime
}  // namespace atlas
