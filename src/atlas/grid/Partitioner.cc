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

namespace atlas {
namespace grid {

using Factory = detail::partitioner::PartitionerFactory;

bool Partitioner::exists( const std::string& type ) {
    return Factory::has(type);
}

Partitioner::Partitioner( const detail::partitioner::Partitioner* partitioner ):
    partitioner_( partitioner ) {
}

Partitioner::Partitioner( const std::string& type, const Grid& grid ):
    partitioner_( Factory::build(type,grid) ) {
}

Partitioner::Partitioner( const std::string& type, const Grid& grid, const size_t nb_partitions):
    partitioner_( Factory::build(type,grid,nb_partitions) ) {
}


} // namespace grid
} // namespace atlas
