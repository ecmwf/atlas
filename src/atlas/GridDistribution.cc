/*
 * (C) Copyright 1996-2014 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */


#include "atlas/Grid.h"
#include "atlas/GridDistribution.h"
#include "atlas/Partitioner.h"

namespace atlas {

GridDistribution::GridDistribution(const Grid& grid)
{
  part_.resize(grid.npts());
  nb_partitions_ = 1;
}

GridDistribution::GridDistribution(const Grid& grid, const Partitioner& partitioner)
{
  part_.resize(grid.npts());
  partitioner.partition(grid,part_.data());
  nb_partitions_ = partitioner.nb_partitions();
}

} // namespace atlas
