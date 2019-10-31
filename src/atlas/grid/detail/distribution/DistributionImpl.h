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

#include "atlas/library/config.h"

namespace atlas {

template< typename T >
class vector {
public:
    using iterator = T*;
    using const_iterator = T const*;
    vector() : data_( nullptr ), size_( 0 ) {}
    vector( idx_t N ) : data_( new T[N] ), size_( N ) {}
    vector( idx_t N, const T& value ) : vector(N) {
        for( idx_t i=0; i<size_; ++i ) {
            data_[i] = value;
        }
    }
    vector( vector&& other ) : data_( other.data_ ), size_( other.size_ ) {
        other.data_ = nullptr;
        other.size_ = 0;
    }
    ~vector() {
        if( data_ ) {
          delete[] data_;
        }
    }
    T& at(idx_t i) { return data_[i]; }
    T const& at(idx_t i) const { return data_[i]; }
    T& operator[](idx_t i) { return data_[i]; }
    T const& operator[](idx_t i) const { return data_[i]; }
    const T* data() const { return data_; }
    T* data() { return data_; }
    idx_t size() const { return size_; }
    template< typename Iter >
    void assign( const Iter& begin, const Iter& end ) {
        throw 1;
    }
    const_iterator begin() const { return data_; }
    const_iterator end() const { return data_+size_; }
    iterator begin() { return data_; }
    iterator end() { return data_+size_; }
    const_iterator cbegin() const { return data_; }
    const_iterator cend() const { return data_+size_; }
private:
    T* data_;
    idx_t size_;
};

class Grid;
namespace grid {
class Partitioner;
}

}  // namespace atlas

namespace atlas {
namespace grid {

class DistributionImpl : public util::Object {
public:
    using partition_t = atlas::vector<int>;

    DistributionImpl( const Grid& );

    DistributionImpl( const Grid&, const Partitioner& );

    DistributionImpl( int nb_partitions, idx_t npts, int partition[], int part0 = 0 );

    DistributionImpl( int nb_partitions, partition_t&& partition );

    virtual ~DistributionImpl();

    int partition( const gidx_t gidx ) const { return part_[gidx]; }

    const partition_t& partition() const { return part_; }

    idx_t nb_partitions() const { return nb_partitions_; }

    operator const partition_t&() const { return part_; }

    const int* data() const { return part_.data(); }

    const std::vector<idx_t>& nb_pts() const { return nb_pts_; }

    idx_t max_pts() const { return max_pts_; }
    idx_t min_pts() const { return min_pts_; }

    const std::string& type() const { return type_; }

    void print( std::ostream& ) const;

private:
    idx_t nb_partitions_;
    partition_t part_;
    std::vector<idx_t> nb_pts_;
    idx_t max_pts_;
    idx_t min_pts_;
    std::string type_;
};

extern "C" {
DistributionImpl* atlas__GridDistribution__new( idx_t npts, int part[], int part0 );
void atlas__GridDistribution__delete( DistributionImpl* This );
}

}  // namespace grid
}  // namespace atlas
