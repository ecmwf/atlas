/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "Logging.h"

#include <iostream>
#include <exception>

#include "eckit/log/Channel.h"

#include "atlas/library/Library.h"
#include "atlas/parallel/omp/omp.h"

//-----------------------------------------------------------------------------------------------------------

namespace atlas {
namespace runtime {
namespace trace {

//-----------------------------------------------------------------------------------------------------------

bool Control::enabled() {
    return atlas_omp_get_thread_num() == 0;
}

class LoggingState {
private:
    std::ostream* channel_;

    LoggingState() { channel_ = &atlas::Library::instance().traceChannel(); }

public:
    static eckit::Channel& empty_channel() {
        static eckit::Channel channel;
        return channel;
    }

    static LoggingState& instance() {
        static LoggingState channel;
        return channel;
    }

    operator std::ostream&() { return *channel_; }
    operator std::ostream*() { return channel_; }

    void set(std::ostream& channel) { channel_ = &channel; }
    void set(bool state) {
        if (state == false) {
            channel_ = &empty_channel();
        }
    }
};

//-----------------------------------------------------------------------------------------------------------

Logging::Logging(bool state): previous_state_(LoggingState::instance()) {
    LoggingState::instance().set(state);
}

Logging::Logging(std::ostream& channel): previous_state_(LoggingState::instance()) {
    LoggingState::instance().set(channel);
}

Logging::~Logging() {
    LoggingState::instance().set(*previous_state_);
}

std::ostream& Logging::channel() {
    return LoggingState::instance();
}

bool Logging::enabled() {
    return LoggingState::instance();
}

void Logging::start(const std::string& title) {
    if (enabled()) {
        channel() << title << " ..." << std::endl;
    }
}

void Logging::stop(const std::string& title, double seconds) {
    if (enabled()) {
        if (!std::uncaught_exception()){
            channel() << title << " ... done : " << seconds << "s" << std::endl;
        }
    }
}
//-----------------------------------------------------------------------------------------------------------

std::ostream& NoLogging::channel() {
    return LoggingState::empty_channel();
}

//-----------------------------------------------------------------------------------------------------------

void LoggingResult::stop(const std::string& title, double seconds) {
    if (enabled()) {
        channel() << title << " : " << seconds << "s" << std::endl;
    }
}

//-----------------------------------------------------------------------------------------------------------

}  // namespace trace
}  // namespace runtime
}  // namespace atlas
