/*
 * (C) Copyright 1996-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */


#ifndef atlas_partitioner_Partitioner_h
#define atlas_partitioner_Partitioner_h

#include <string>
#include "eckit/config/Configuration.h"
#include "eckit/memory/NonCopyable.h"
#include "eckit/memory/ScopedPtr.h"
#include "atlas/grid/Domain.h"
#include "atlas/grid/GridDistribution.h"
#include "atlas/grid/Structured.h"
#include "atlas/grid/partitioners/Partitioner.h"
#include "atlas/mesh/Mesh.h"


namespace atlas {
namespace grid { class Structured; }
namespace mesh { class Mesh; }
}


namespace atlas {
namespace grid {
namespace partitioners {


class PartitionerFromPrePartitionedMesh : public Partitioner {
public:

    PartitionerFromPrePartitionedMesh(const Grid& grid) : Partitioner(grid) {}
    PartitionerFromPrePartitionedMesh(const Grid& grid, const size_t nb_partitions) : Partitioner(grid, nb_partitions) {}

    virtual ~PartitionerFromPrePartitionedMesh() {}

    GridDistribution* distribution() const {
        ASSERT(prePartitionedMesh_);
        return distributionFromPrePartitionedMesh();
    }

    /**
     * @param meshSource mesh already partitioned
     * @param gridTarget grid to be distributed
     * @param includesNorthPole
     * @param includesSouthPole
     */
    void setup(const mesh::Mesh::Ptr prePartitionedMesh, const Domain& prePartitionedDomain=Domain::makeGlobal()) {
        ASSERT(prePartitionedMesh);
        prePartitionedMesh_ = prePartitionedMesh;
        prePartitionedDomain_ = prePartitionedDomain;
    }

    /**
     * @brief Create a GridDistribution, placing nodes in the same partitions as a
     * given pre-partitioned mesh.
     * @return grid partitioner
     */
    virtual GridDistribution* distributionFromPrePartitionedMesh() const = 0;

    void partition(int part[]) const {
        // FIXME don't understand if this is needed!
        NOTIMP;
    }

protected:

    mesh::Mesh::Ptr prePartitionedMesh_;
    Domain prePartitionedDomain_;

};


}  // partitioners
}  // grid
}  // atlas


#endif
