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

#include <string>
#include <vector>

#include "atlas/util/Object.h"
#include "atlas/util/vector.h"
#include "atlas/util/Config.h"
#include "atlas/library/config.h"


namespace eckit {
class Parametrisation;
}

namespace atlas {

class Grid;
namespace grid {
class Partitioner;
}
}  // namespace atlas

namespace atlas {
namespace grid {

class DistributionImpl : public util::Object {
public:
    using Config = atlas::util::Config;
    using partition_t = atlas::vector<int>;

    DistributionImpl( const Grid& );
    DistributionImpl( const Grid&, const eckit::Parametrisation & );

    DistributionImpl( const Grid&, const Partitioner& );

    DistributionImpl( int nb_partitions, idx_t npts, int partition[], int part0 = 0 );

    DistributionImpl( int nb_partitions, partition_t&& partition );

    virtual ~DistributionImpl();

    int partition( const gidx_t gidx ) const 
    { 
      if (part_.size () > 0)
        return part_[gidx];
      idx_t iblock = gidx / blocksize_;
      return (iblock * nb_partitions_) / nb_blocks_;
    }

    const partition_t& partition() const { checkPartition (); return part_; }

    idx_t nb_partitions() const { return nb_partitions_; }

    operator const partition_t&() const { return part_; }

    const int* data() const { checkPartition (); return part_.data(); }

    const std::vector<idx_t>& nb_pts() const { return nb_pts_; }

    idx_t max_pts() const { return max_pts_; }
    idx_t min_pts() const { return min_pts_; }

    const std::string& type() const { return type_; }

    void print( std::ostream& ) const;

    size_t footprint () const
    {
      return nb_pts_.size () * sizeof (nb_pts_[0]) + part_.size () * sizeof (part_[0]);
    }


private:

    void setupWithPartitioner (const Grid &, const Partitioner &);

    void checkPartition () const
    {
      if (part_.size () == 0)
        throw_Exception ("partition array of distribution is empty");
    }
    // For trivial partitionning
    size_t gridsize_ = 0;
    idx_t nb_partitions_ = 0;
    size_t blocksize_ = 0;
    idx_t nb_blocks_ = 0;

    partition_t part_;
    std::vector<idx_t> nb_pts_;
    idx_t max_pts_;
    idx_t min_pts_;
    std::string type_;
};

extern "C" {
DistributionImpl* atlas__GridDistribution__new( idx_t npts, int part[], int part0 );
void atlas__GridDistribution__delete( DistributionImpl* This );
void atlas__GridDistribution__nb_pts( DistributionImpl* This, idx_t nb_pts[] );
idx_t atlas__atlas__GridDistribution__nb_partitions( DistributionImpl* This );
}

}  // namespace grid
}  // namespace atlas
