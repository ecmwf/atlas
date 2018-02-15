/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#pragma once

#include "atlas/grid/Distribution.h"
#include "atlas/grid/Grid.h"
#include "atlas/grid/detail/partitioner/Partitioner.h"
#include "eckit/memory/SharedPtr.h"

namespace atlas {
namespace grid {

// ------------------------------------------------------------------

class Partitioner {
public:
    using Config         = eckit::Parametrisation;
    using Implementation = detail::partitioner::Partitioner;

public:
    static bool exists( const std::string& type );

public:
    Partitioner() {}
    Partitioner( const Implementation* );
    Partitioner( const std::string& type );
    Partitioner( const std::string& type, const size_t nb_partitions );
    Partitioner( const Config& );

    operator bool() const { return partitioner_; }

    void partition( const Grid& grid, int part[] ) const;

    Distribution partition( const Grid& grid ) const { return Distribution( grid, *this ); }

    size_t nb_partitions() const { return partitioner_->nb_partitions(); }

    std::string type() const { return partitioner_->type(); }

    Implementation const* get() const { return partitioner_.get(); }

private:
    eckit::SharedPtr<const Implementation> partitioner_;
};

// ------------------------------------------------------------------

class MatchingMeshPartitioner : public Partitioner {
public:
    using Config = eckit::Parametrisation;

public:
    static bool exists( const std::string& type );

public:
    MatchingMeshPartitioner();
    MatchingMeshPartitioner( const Mesh& mesh, const Config& config );
};

// ------------------------------------------------------------------

extern "C" {
Partitioner::Implementation* atlas__grid__Partitioner__new( const Partitioner::Config* config );
Partitioner::Implementation* atlas__grid__MatchingMeshPartitioner__new( const Mesh::Implementation* mesh,
                                                                        const Partitioner::Config* config );
void atlas__grid__Partitioner__delete( Partitioner::Implementation* This );
Distribution::impl_t* atlas__grid__Partitioner__partition( const Partitioner::Implementation* This,
                                                           const Grid::Implementation* grid );
}

}  // namespace grid
}  // namespace atlas
