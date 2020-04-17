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

#include <vector>

#include "atlas/grid/detail/distribution/DistributionImpl.h"
#include "atlas/library/config.h"
#include "atlas/util/ObjectHandle.h"
#include "atlas/util/Config.h"

namespace atlas {
class Grid;
namespace grid {
class Partitioner;
class DistributionImpl;
class Partitioner;
}  // namespace grid
}  // namespace atlas

namespace atlas {
namespace grid {

class Distribution : DOXYGEN_HIDE( public util::ObjectHandle<DistributionImpl> ) {
    friend class Partitioner;

public:
    using partition_t = DistributionImpl::partition_t;
    using Config = DistributionImpl::Config;

    using Handle::Handle;
    Distribution() = default;

    Distribution( const Grid& );

    Distribution( const Grid&, const Config& );

    Distribution( const Grid&, const Partitioner& );

    Distribution( int nb_partitions, idx_t npts, int partition[], int part0 = 0 );

    Distribution( int nb_partitions, partition_t&& partition );

    ~Distribution();

    // This method has to be inlined
    int partition( const gidx_t gidx ) const
    {
      return get()->partition( gidx );
    }

    size_t footprint () const
    {
      return get()->footprint ();
    }

    const partition_t& partition() const;

    idx_t nb_partitions() const;

    operator const partition_t&() const;

    const int* data() const;

    const std::vector<idx_t>& nb_pts() const;

    idx_t max_pts() const;
    idx_t min_pts() const;

    const std::string& type() const;

    friend std::ostream& operator<<( std::ostream& os, const Distribution& distribution );
};

}  // namespace grid
}  // namespace atlas
