/*
 * (C) Copyright 1996-2014 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */


#ifndef atlas_GridDistribution_h
#define atlas_GridDistribution_h

#include <vector>
#include "eckit/memory/Owned.h"
#include "eckit/memory/SharedPtr.h"
#include "atlas/atlas_config.h"

namespace atlas {

class Partitioner;

class GridDistribution: public eckit::Owned
{
public:
  typedef eckit::SharedPtr<GridDistribution> Ptr;
public:

  GridDistribution(const Grid& grid);

  GridDistribution(const Grid& grid, const Partitioner& partitioner);

  virtual ~GridDistribution() {}

  int partition(const gidx_t gidx) const { return part_[gidx]; }

  const std::vector<int>& partition() const { return part_; }

  size_t nb_partitions() const { return nb_partitions_; }

  operator const std::vector<int>&() const { return part_; }

  const int* data() const { return part_.data(); }

private:
  size_t nb_partitions_;
  std::vector<int> part_;
};

} // namespace atlas

#endif // atlas_GridDistribution_h
