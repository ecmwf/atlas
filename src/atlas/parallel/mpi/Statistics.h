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
#include "atlas/runtime/Timer.h"

#define ATLAS_MPI_STATS(...)

#if ATLAS_HAVE_TIMINGS

#include "atlas/util/detail/BlackMagic.h"

#undef ATLAS_MPI_STATS
#define ATLAS_MPI_STATS(...) ATLAS_MPI_STATS_( __ATLAS__NARG(__VA_ARGS__), ##__VA_ARGS__ )
#define ATLAS_MPI_STATS_(N, ...) __ATLAS__SPLICE( ATLAS_MPI_STATS_, N)(__VA_ARGS__)
#define ATLAS_MPI_STATS_1( statistics_enum ) \
  for( ::atlas::parallel::mpi::Statistics __ATLAS__SPLICE( stats, __LINE__ ) \
    /*args=*/(Here(), ::atlas::parallel::mpi::StatisticsEnum::__ATLAS_STRINGIFY(statistics_enum) );\
    __ATLAS__SPLICE( stats, __LINE__ ) .running(); \
    __ATLAS__SPLICE( stats, __LINE__ ) .stop() )
#define ATLAS_MPI_STATS_2(statistics_enum,title) \
      for( ::atlas::parallel::mpi::Statistics __ATLAS__SPLICE( stats, __LINE__ ) \
        /*args=*/(Here(), ::atlas::parallel::mpi::StatisticsEnum::__ATLAS_STRINGIFY(statistics_enum), title);\
        __ATLAS__SPLICE( stats, __LINE__ ) .running(); \
        __ATLAS__SPLICE( stats, __LINE__ ) .stop() )

#endif


namespace atlas {
namespace parallel {
namespace mpi {

struct StatisticsTimerTraits {
    using Barriers = runtime::timer::TimerBarriers;
    using Logging  = runtime::timer::TimerNoLogging;
    using Timings  = runtime::timer::Timings;
    using Nesting  = runtime::timer::TimerNesting;
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
    "mpi-broadcast",
    "mpi-allreduce",
    "mpi-allgather",
    "mpi-alltoall",
    "mpi-reduce",
    "mpi-gather",
    "mpi-scatter",
    "mpi-barrier",
    "mpi-sendreceive",
    "mpi-isend",
    "mpi-ireceive",
    "mpi-wait"
  };
  return names[ static_cast<size_t>(c) ];
}

class Statistics : public runtime::timer::TimerT< StatisticsTimerTraits > {
    using Base = runtime::timer::TimerT< StatisticsTimerTraits >;
public:
    Statistics( const eckit::CodeLocation& loc, StatisticsEnum c ) :
      Base( loc, name(c), make_labels(c), Logging::channel() ) {
    }
    Statistics( const eckit::CodeLocation& loc, StatisticsEnum c, const std::string& title ) :
      Base( loc, title, make_labels(c), Logging::channel() ) {
    }
private:
    static std::vector<std::string> make_labels( StatisticsEnum c) {
      return {"mpi", name(c)};
    }
};

} // namespace mpi
} // namespace parallel
} // namespace atlas
