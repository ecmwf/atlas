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

#include <array>
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/runtime/Trace.h"

#define ATLAS_TRACE_MPI(...)

#if ATLAS_HAVE_TRACE

#include "atlas/util/detail/BlackMagic.h"

#undef ATLAS_TRACE_MPI
#define ATLAS_TRACE_MPI(...) ATLAS_TRACE_MPI_( ::atlas::parallel::mpi::Statistics, Here(), ##__VA_ARGS__ )
#define ATLAS_TRACE_MPI_( Type, loc, enum_var, ... ) \
  __ATLAS_TYPE_SCOPE( Type, loc, ::atlas::parallel::mpi::StatisticsEnum::__ATLAS_STRINGIFY(enum_var), ##__VA_ARGS__ )

#endif


namespace atlas {
namespace parallel {
namespace mpi {

struct StatisticsTimerTraits {
    using Barriers = runtime::timer::TimerBarriers;
    using Tracing  = runtime::timer::TimerTracingNone;
};


enum class StatisticsEnum {
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

static const std::string& name(StatisticsEnum c) {
  static std::array<std::string, static_cast<size_t>(StatisticsEnum::_COUNT_)> names {
    "mpi.broadcast",
    "mpi.allreduce",
    "mpi.allgather",
    "mpi.alltoall",
    "mpi.reduce",
    "mpi.gather",
    "mpi.scatter",
    "mpi.barrier",
    "mpi.sendreceive",
    "mpi.isend",
    "mpi.ireceive",
    "mpi.wait"
  };
  return names[ static_cast<size_t>(c) ];
}

class Statistics : public runtime::timer::TimerT< StatisticsTimerTraits > {
    using Base = runtime::timer::TimerT< StatisticsTimerTraits >;
public:
    Statistics( const eckit::CodeLocation& loc, StatisticsEnum c ) :
      Base( loc, name(c), make_labels(c) ) {
    }
    Statistics( const eckit::CodeLocation& loc, StatisticsEnum c, const std::string& title ) :
      Base( loc, title, make_labels(c) ) {
    }
private:
    static std::vector<std::string> make_labels( StatisticsEnum c) {
      return {"mpi", name(c)};
    }
};

} // namespace mpi
} // namespace parallel
} // namespace atlas
