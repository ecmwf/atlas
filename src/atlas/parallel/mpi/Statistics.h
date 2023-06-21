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

#include <array>

#include "atlas/parallel/mpi/mpi.h"
#include "atlas/runtime/Trace.h"

#define ATLAS_TRACE_MPI(...)

#if ATLAS_HAVE_TRACE

#include "atlas/library/detail/BlackMagic.h"

#undef ATLAS_TRACE_MPI
#define ATLAS_TRACE_MPI(...) ATLAS_TRACE_MPI_(::atlas::mpi::Trace, Here(), __VA_ARGS__)
#define ATLAS_TRACE_MPI_(Type, location, operation, ...) \
    __ATLAS_TYPE_SCOPE(Type, location, __ATLAS_TRACE_MPI_ENUM(operation) __ATLAS_COMMA_ARGS(__VA_ARGS__))

#define __ATLAS_TRACE_MPI_ENUM(operation) ::atlas::mpi::Operation::__ATLAS_STRINGIFY(operation)

#endif

namespace atlas {
namespace mpi {

struct StatisticsTimerTraits {
    using Barriers = runtime::trace::NoBarriers;
    using Tracing  = runtime::trace::NoLogging;
};

enum class Operation
{
    BROADCAST,
    ALLREDUCE,
    ALLGATHER,
    ALLTOALL,
    REDUCE,
    GATHER,
    SCATTER,
    BARRIER,
    SENDRECEIVE,
    ISEND,
    IRECEIVE,
    WAIT,
    _COUNT_
};

static const std::string& name(Operation c) {
    static std::array<std::string, static_cast<size_t>(Operation::_COUNT_)> names{
        "mpi.broadcast", "mpi.allreduce", "mpi.allgather",   "mpi.alltoall", "mpi.reduce",   "mpi.gather",
        "mpi.scatter",   "mpi.barrier",   "mpi.sendreceive", "mpi.isend",    "mpi.ireceive", "mpi.wait"};
    return names[static_cast<size_t>(c)];
}

class Trace : public runtime::trace::TraceT<StatisticsTimerTraits> {
    using Base = runtime::trace::TraceT<StatisticsTimerTraits>;

public:
    Trace(const eckit::CodeLocation& loc, Operation c): Base(loc, name(c), make_labels(c)) {}
    Trace(const eckit::CodeLocation& loc, Operation c, const std::string& title): Base(loc, title, make_labels(c)) {}

private:
    static std::vector<std::string> make_labels(Operation c) { return {"mpi", name(c)}; }
};

}  // namespace mpi
}  // namespace atlas
