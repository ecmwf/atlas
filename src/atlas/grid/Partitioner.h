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
#include "atlas/grid/detail/partitioners/Partitioner.h"
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

    Partitioner( const detail::partitioners::Partitioner* );
    Partitioner( const std::string& type, const Grid& );
    Partitioner( const std::string& type, const Grid&, const size_t nb_partitions);

    void partition( int part[] ) const { partitioner_->partition(part); }

    virtual Distribution distribution() const { return partitioner_->distribution(); }

    size_t nb_partitions() const { return partitioner_->nb_partitions(); }
    const Grid& grid() const { return partitioner_->grid(); }

private:

    eckit::SharedPtr<const detail::partitioners::Partitioner> partitioner_;
};

// ------------------------------------------------------------------

} // namespace grid
} // namespace atlas
