/*
 * (C) Copyright 1996-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */


#ifndef atlas_grid_partitioners_PrePartitionedBruteForce_h
#define atlas_grid_partitioners_PrePartitionedBruteForce_h

#include "PartitionerFromPrePartitionedMesh.h"


namespace atlas {
namespace grid {
namespace partitioners {


class PrePartitionedBruteForce : public PartitionerFromPrePartitionedMesh {
public:

    PrePartitionedBruteForce(const Grid& grid) : PartitionerFromPrePartitionedMesh(grid) {}
    PrePartitionedBruteForce(const Grid& grid, const size_t nb_partitions) : PartitionerFromPrePartitionedMesh(grid, nb_partitions) {}

    PrePartitionedBruteForce(const Grid& grid, 
      const mesh::Mesh& mesh, const Domain& domain=Domain::makeGlobal() ) : 
      PartitionerFromPrePartitionedMesh(grid,mesh,domain) {}

      virtual void partition( int part[] ) const;

    /**
     * @brief Create a GridDistribution, placing nodes in the same partitions as a
     * given pre-partitioned mesh. The method is very simple and only a starting point,
     * assigning a partition number (MPI rank) by checking every source mesh element,
     * for every target grid node.
     * @param gridTarget grid to be distributed
     * @param meshSource mesh already partitioned
     * @return grid partitioner
     */
    // GridDistribution* distributionFromPrePartitionedMesh() const;

};


}  // partitioners
}  // grid
}  // atlas


#endif
