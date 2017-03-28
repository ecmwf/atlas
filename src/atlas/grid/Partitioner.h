/*
 * (C) Copyright 1996-2017 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */


#pragma once

#include "eckit/memory/SharedPtr.h"
#include "atlas/grid/detail/partitioner/Partitioner.h"
#include "atlas/grid/Distribution.h"
#include "atlas/grid/Grid.h"


namespace atlas {
namespace grid {

// ------------------------------------------------------------------

class Partitioner {

public:

  using Config = eckit::Parametrisation;

public:

  static bool exists(const std::string& type);

public:

    Partitioner() {}
    Partitioner( const detail::partitioner::Partitioner* );
    Partitioner( const std::string& type );
    Partitioner( const std::string& type, const size_t nb_partitions );

    void partition( const Grid& grid, int part[] ) const { partitioner_->partition(grid,part); }

    Distribution partition( const Grid& grid ) const { return Distribution(grid,*this); }

    size_t nb_partitions() const { return partitioner_->nb_partitions(); }

private:

    eckit::SharedPtr<const detail::partitioner::Partitioner> partitioner_;
};

// ------------------------------------------------------------------

class MatchingMeshPartitioner: public Partitioner {

public:

  using Config = eckit::Parametrisation;

public:

  static bool exists(const std::string& type);

public:

    MatchingMeshPartitioner();
    MatchingMeshPartitioner( const mesh::Mesh& mesh, const Config& config );
};

// ------------------------------------------------------------------

} // namespace grid
} // namespace atlas
