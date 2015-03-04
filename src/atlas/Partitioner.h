/*
 * (C) Copyright 1996-2014 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */


#ifndef atlas_Partitioner_h
#define atlas_Partitioner_h

#include "atlas/Grid.h"

namespace atlas {

class GridDistribution;

class Partitioner {
public:

  Partitioner(const Grid& grid);
  virtual ~Partitioner();

  virtual void partition( int part[] ) const = 0;

  virtual GridDistribution* distribution() const;

public:

  size_t nb_partitions() const;
  void set_nb_partition(const size_t n);
  const Grid& grid() const { return grid_; }

private:

  size_t nb_partitions_;
  const Grid& grid_;
};

} // namespace atlas

#endif // atlas_Partitioner_h
