/*
 * (C) Copyright 1996-2014 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "eckit/mpi/mpi.h"
#include "atlas/Partitioner.h"
#include "atlas/GridDistribution.h"

namespace atlas {

Partitioner::Partitioner(const Grid& grid): grid_(grid), nb_partitions_(eckit::mpi::size())
{ }

Partitioner::~Partitioner()
{ }

size_t Partitioner::nb_partitions() const
{
  return nb_partitions_;
}

void Partitioner::set_nb_partition(const size_t n)
{
  nb_partitions_ = n;
}

GridDistribution* Partitioner::distribution() const
{
  return new GridDistribution(*this);
}

} // namespace atlas
