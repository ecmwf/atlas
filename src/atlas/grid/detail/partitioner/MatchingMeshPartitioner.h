/*
 * (C) Copyright 1996-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */


#pragma once

#include <string>
#include "eckit/config/Configuration.h"
#include "eckit/memory/NonCopyable.h"
#include "eckit/memory/ScopedPtr.h"
#include "atlas/grid/Domain.h"
#include "atlas/grid/Distribution.h"
#include "atlas/grid/Grid.h"
#include "atlas/grid/detail/partitioner/Partitioner.h"
#include "atlas/mesh/Mesh.h"


namespace atlas {
  class Mesh;
}


namespace atlas {
namespace grid {
namespace detail {
namespace partitioner {


class MatchingMeshPartitioner : public Partitioner {
public:

    MatchingMeshPartitioner() : Partitioner() { NOTIMP; }
    MatchingMeshPartitioner(const size_t nb_partitions) : Partitioner(nb_partitions) { NOTIMP; }

    MatchingMeshPartitioner(const Mesh& mesh ) :
      Partitioner(mesh.nb_partitions()),
      prePartitionedMesh_(mesh) {
    }

    virtual ~MatchingMeshPartitioner() {}

protected:

    const Mesh prePartitionedMesh_;

};

}  // partitioners
}  // detail
}  // grid
}  // atlas
