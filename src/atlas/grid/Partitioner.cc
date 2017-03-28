/*
 * (C) Copyright 1996-2017 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "atlas/grid/Partitioner.h"
#include "atlas/grid/detail/partitioner/Partitioner.h"
#include "atlas/parallel/mpi/mpi.h"

namespace atlas {
namespace grid {

using Factory = detail::partitioner::PartitionerFactory;

bool Partitioner::exists( const std::string& type ) {
    return Factory::has(type);
}

Partitioner::Partitioner( const detail::partitioner::Partitioner* partitioner ):
    partitioner_( partitioner ) {
}

Partitioner::Partitioner( const std::string& type ):
    partitioner_( Factory::build(type) ) {
}

Partitioner::Partitioner( const std::string& type, const size_t nb_partitions):
    partitioner_( Factory::build(type,nb_partitions) ) {
}

namespace {
detail::partitioner::Partitioner* partitioner_from_config( const Partitioner::Config& config ) {
    std::string type;
    long partitions = parallel::mpi::comm().size();
    if( not config.get("type",type) )
      throw eckit::BadParameter("'type' missing in configuration for Partitioner",Here());
    config.get("partitions",partitions);
    return Factory::build(type,partitions);
}
}

Partitioner::Partitioner( const Config& config ):
    partitioner_( partitioner_from_config(config) ) {
}

MatchingMeshPartitioner::MatchingMeshPartitioner() :
    Partitioner() {
}

grid::detail::partitioner::Partitioner* matching_mesh_partititioner( const mesh::Mesh& mesh, const Partitioner::Config& config ) {
    std::string type;
    if( not config.get("type",type) )
      type = "polygon";
    return MatchedPartitionerFactory::build(type,mesh);
}

MatchingMeshPartitioner::MatchingMeshPartitioner( const mesh::Mesh& mesh, const Config& config ) :
    Partitioner( matching_mesh_partititioner(mesh,config) ) {
}

} // namespace grid
} // namespace atlas
