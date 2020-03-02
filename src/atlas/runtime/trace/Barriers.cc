/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "Barriers.h"

#include <sstream>

#include "atlas/library/Library.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/parallel/omp/omp.h"
#include "atlas/runtime/trace/StopWatch.h"

//-----------------------------------------------------------------------------------------------------------

namespace atlas {
namespace runtime {
namespace trace {

class BarriersState {
private:
    BarriersState() { barriers_ = atlas::Library::instance().traceBarriers(); }
    bool barriers_;
    StopWatch stopwatch_;

public:
    BarriersState( BarriersState const& ) = delete;
    void operator=( BarriersState const& ) = delete;
    static BarriersState& instance() {
        static BarriersState state;
        return state;
    }
    operator bool() const { return barriers_; }
    void set( bool state ) { barriers_ = state; }
    StopWatch& stopwatch() { return stopwatch_; }

    std::string report() const {
        std::stringstream out;
        double time = stopwatch_.elapsed();
        if ( time ) {
            out << "Total time spent in mpi barriers due to load imbalance : " << time << "s" << std::endl;
        }
        return out.str();
    }
};

Barriers::Barriers( bool state ) : previous_state_( BarriersState::instance() ) {
    BarriersState::instance().set( state );
}

Barriers::~Barriers() {
    restore();
}

void Barriers::restore() {
    BarriersState::instance().set( previous_state_ );
}

bool Barriers::state() {
    return BarriersState::instance() && ( atlas_omp_get_num_threads() <= 1 );
}

void Barriers::execute() {
    if ( state() ) {
        BarriersState::instance().stopwatch().start();
        mpi::comm().barrier();
        BarriersState::instance().stopwatch().stop();
    }
}

double Barriers::time() {
    return BarriersState::instance().stopwatch().elapsed();
}

double NoBarriers::time() {
    return BarriersState::instance().stopwatch().elapsed();
}

std::string Barriers::report() {
    return BarriersState::instance().report();
}

std::string NoBarriers::report() {
    return BarriersState::instance().report();
}

}  // namespace trace
}  // namespace runtime
}  // namespace atlas
