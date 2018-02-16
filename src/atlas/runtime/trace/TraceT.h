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
#include <string>
#include <vector>

#include "atlas/runtime/trace/Nesting.h"
#include "atlas/runtime/trace/StopWatch.h"
#include "atlas/runtime/trace/Timings.h"

//-----------------------------------------------------------------------------------------------------------

namespace eckit {
class CodeLocation;
class Configuration;
}  // namespace eckit

//-----------------------------------------------------------------------------------------------------------

namespace atlas {
namespace runtime {
namespace trace {

//-----------------------------------------------------------------------------------------------------------

template <typename TraceTraits>
class TraceT {
public:
    using Barriers = typename TraceTraits::Barriers;
    using Tracing  = typename TraceTraits::Tracing;
    using Labels   = std::vector<std::string>;

public:  // static methods
    static std::string report();
    static std::string report( const eckit::Configuration& config );

public:
    TraceT( const eckit::CodeLocation& );
    TraceT( const eckit::CodeLocation&, const std::string& title );
    TraceT( const eckit::CodeLocation&, const std::string& title, const Labels& );

    ~TraceT();

    bool running() const;

    void start();

    void stop();

    void pause();

    void resume();

    double elapsed() const;

private:  // types
    using Identifier = Timings::Identifier;

private:  // member functions
    void barrier() const;

    void updateTimings() const;

    void registerTimer();

private:  // member data
    bool running_{true};
    StopWatch stopwatch_;
    eckit::CodeLocation loc_;
    std::string title_;
    Identifier id_;
    Nesting nesting_;
    Labels labels_;
};

//-----------------------------------------------------------------------------------------------------------
// Definitions

template <typename TraceTraits>
inline TraceT<TraceTraits>::TraceT( const eckit::CodeLocation& loc, const std::string& title ) :
    loc_( loc ),
    title_( title ),
    nesting_( loc ) {
    start();
}

template <typename TraceTraits>
inline TraceT<TraceTraits>::TraceT( const eckit::CodeLocation& loc ) :
    loc_( loc ),
    title_( loc_ ? loc_.func() : "" ),
    nesting_( loc_ ) {
    start();
}

template <typename TraceTraits>
inline TraceT<TraceTraits>::TraceT( const eckit::CodeLocation& loc, const std::string& title, const Labels& labels ) :
    loc_( loc ),
    title_( title ),
    nesting_( loc ),
    labels_( labels ) {
    start();
}

template <typename TraceTraits>
inline TraceT<TraceTraits>::~TraceT() {
    stop();
}

template <typename TraceTraits>
inline void TraceT<TraceTraits>::barrier() const {
    Barriers::execute();
}

template <typename TraceTraits>
inline void TraceT<TraceTraits>::registerTimer() {
    id_ = Timings::add( loc_, nesting_, title_ + ( Barriers::state() ? " [b]" : "" ), labels_ );
}

template <typename TraceTraits>
inline void TraceT<TraceTraits>::updateTimings() const {
    Timings::update( id_, stopwatch_.elapsed() );
}

template <typename TraceTraits>
inline bool TraceT<TraceTraits>::running() const {
    return running_;
}

template <typename TraceTraits>
inline void TraceT<TraceTraits>::start() {
    registerTimer();
    Tracing::start( title_ );
    barrier();
    stopwatch_.start();
}

template <typename TraceTraits>
inline void TraceT<TraceTraits>::stop() {
    if ( running_ ) {
        barrier();
        stopwatch_.stop();
        nesting_.stop();
        updateTimings();
        Tracing::stop( title_, stopwatch_.elapsed() );
        running_ = false;
    }
}

template <typename TraceTraits>
inline void TraceT<TraceTraits>::pause() {
    if ( running_ ) {
        barrier();
        stopwatch_.stop();
        nesting_.stop();
    }
}

template <typename TraceTraits>
inline void TraceT<TraceTraits>::resume() {
    if ( running_ ) {
        barrier();
        nesting_.start();
        stopwatch_.start();
    }
}

template <typename TraceTraits>
inline double TraceT<TraceTraits>::elapsed() const {
    return stopwatch_.elapsed();
}

template <typename TraceTraits>
inline std::string TraceT<TraceTraits>::report() {
    return Timings::report() + Barriers::report();
}

template <typename TraceTraits>
inline std::string TraceT<TraceTraits>::report( const eckit::Configuration& config ) {
    return Timings::report( config ) + Barriers::report();
}

//-----------------------------------------------------------------------------------------------------------

}  // namespace trace
}  // namespace runtime
}  // namespace atlas
