/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor
 * does it submit to any jurisdiction.
 */

#pragma once

#include "atlas/grid/detail/grid/Grid.h"

//---------------------------------------------------------------------------------------------------------------------

namespace atlas {
namespace grid {

//---------------------------------------------------------------------------------------------------------------------

class IteratorXY {
public:
    IteratorXY( detail::grid::Grid::IteratorXY* iterator ) : iterator_( iterator ) {}

    bool next( PointXY& xy ) { return iterator_->next( xy ); }

    PointXY operator*() const { return iterator_->operator*(); }

    const IteratorXY& operator++() {
        iterator_->operator++();
        return *this;
    }

    bool operator==( const IteratorXY& other ) const { return iterator_->operator==( *other.iterator_ ); }

    bool operator!=( const IteratorXY& other ) const { return iterator_->operator!=( *other.iterator_ ); }

private:
    std::unique_ptr<detail::grid::Grid::IteratorXY> iterator_;
};

//---------------------------------------------------------------------------------------------------------------------

class IteratorLonLat {
public:
    IteratorLonLat( detail::grid::Grid::IteratorLonLat* iterator ) : iterator_( iterator ) {}

    bool next( PointLonLat& lonlat ) { return iterator_->next( lonlat ); }

    PointLonLat operator*() const { return iterator_->operator*(); }

    const IteratorLonLat& operator++() {
        iterator_->operator++();
        return *this;
    }

    bool operator==( const IteratorLonLat& other ) const { return iterator_->operator==( *other.iterator_ ); }

    bool operator!=( const IteratorLonLat& other ) const { return iterator_->operator!=( *other.iterator_ ); }

private:
    std::unique_ptr<detail::grid::Grid::IteratorLonLat> iterator_;
};

//---------------------------------------------------------------------------------------------------------------------

}  // namespace grid
}  // namespace atlas
