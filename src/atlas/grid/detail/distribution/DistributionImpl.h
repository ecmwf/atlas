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

#include <string>
#include <vector>

#include "atlas/library/config.h"
#include "atlas/util/Config.h"
#include "atlas/util/Object.h"


namespace eckit {
class Parametrisation;
}

namespace atlas {

class Grid;
namespace grid {
namespace detail {

namespace partitioner {
class Partitioner;
}
namespace grid {
class Grid;
}
}  // namespace detail
}  // namespace grid
}  // namespace atlas

namespace atlas {
namespace grid {

class DistributionImpl : public util::Object {
public:
    using Config = atlas::util::Config;
    virtual ~DistributionImpl() {}
    virtual int partition(const gidx_t gidx) const = 0;
    virtual bool functional() const                = 0;
    virtual idx_t nb_partitions() const            = 0;
    virtual gidx_t size() const                    = 0;

    virtual const std::vector<idx_t>& nb_pts() const = 0;

    virtual idx_t max_pts() const = 0;
    virtual idx_t min_pts() const = 0;

    virtual const std::string& type() const = 0;

    virtual void print(std::ostream&) const = 0;

    virtual size_t footprint() const = 0;

    virtual void hash(eckit::Hash&) const = 0;

    virtual void partition(gidx_t begin, gidx_t end, int partitions[]) const = 0;
};


extern "C" {
DistributionImpl* atlas__GridDistribution__new(idx_t size, int part[], int part0);
DistributionImpl* atlas__GridDistribution__new__Grid_Config(const detail::grid::Grid* grid,
                                                            const eckit::Parametrisation* config);
DistributionImpl* atlas__GridDistribution__new__Grid_Partitioner(const detail::grid::Grid* grid,
                                                                 const detail::partitioner::Partitioner* partitioner);
int atlas__GridDistribution__partition_int32(DistributionImpl* dist, int i);
int atlas__GridDistribution__partition_int64(DistributionImpl* dist, long i);
void atlas__GridDistribution__delete(DistributionImpl* This);
void atlas__GridDistribution__nb_pts(DistributionImpl* This, idx_t nb_pts[]);
idx_t atlas__atlas__GridDistribution__nb_partitions(DistributionImpl* This);
}

}  // namespace grid
}  // namespace atlas
