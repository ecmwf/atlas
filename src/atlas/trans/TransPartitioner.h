/*
 * (C) Copyright 1996-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef atlas_trans_TransPartitioner_h
#define atlas_trans_TransPartitioner_h

#include "atlas/Partitioner.h"

namespace atlas {
namespace grids {
  class ReducedGrid;
}
namespace trans {
  class Trans;
}
}

namespace atlas {
namespace trans {

class TransPartitioner: public Partitioner
{
public:

  /// @brief Constructor
  TransPartitioner( const Grid& grid, const Trans& trans );

  /// @brief Constructor
  /// This constructor allocates a new Trans, but without the computations
  /// of the spectral coefficients (LDGRIDONLY=TRUE)
  TransPartitioner( const Grid& grid );

  TransPartitioner( const Grid& grid, const size_t nb_partitions );

  virtual ~TransPartitioner();

  virtual void partition(int part[]) const;

  int nb_bands() const;

  int nb_regions(int b) const;

private:
  mutable Trans* t_;
  bool owned_;
};

}
}

#endif // atlas_trans_TransPartitioner_h
