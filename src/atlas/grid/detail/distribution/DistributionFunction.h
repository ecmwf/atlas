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

#include "atlas/grid/detail/distribution/DistributionImpl.h"


namespace atlas {
namespace grid {
namespace detail {
namespace distribution {

class DistributionFunction : public DistributionImpl {
public:
    DistributionFunction( const Grid& ) : DistributionImpl() {}
    bool functional() const override { return true; }
    size_t footprint() const override { return nb_pts_.size() * sizeof( nb_pts_[0] ); }
    const std::string& type() const override { return nb_partitions_ == 1 ? serial : type_; }
    idx_t nb_partitions() const { return nb_partitions_; }

    virtual const std::vector<idx_t>& nb_pts() const { return nb_pts_; }

    virtual idx_t max_pts() const { return max_pts_; }
    virtual idx_t min_pts() const { return min_pts_; }

    virtual void print( std::ostream& ) const;

    virtual gidx_t size() const { return size_; }

    void hash( eckit::Hash& ) const override;

protected:
    gidx_t size_;
    idx_t nb_partitions_;
    std::vector<idx_t> nb_pts_;
    idx_t max_pts_;
    idx_t min_pts_;
    std::string type_{"functional"};
    std::string serial{"serial"};
};


template <typename Derived>
class DistributionFunctionT : public DistributionFunction {
public:
    DistributionFunctionT( const Grid& grid ) : DistributionFunction( grid ) {}
    int partition( gidx_t index ) const override { return static_cast<const Derived*>( this )->function( index ); }

    void partition( gidx_t begin, gidx_t end, int partitions[] ) const override {
        size_t i = 0;
        for ( gidx_t n = begin; n < end; ++n, ++i ) {
            partitions[i] = static_cast<const Derived*>( this )->function( n );
        }
    }
};

}  // namespace distribution
}  // namespace detail
}  // namespace grid
}  // namespace atlas
