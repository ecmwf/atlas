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

#include <iterator>
#include <memory>

#include "atlas/grid/detail/grid/Grid.h"
#include "atlas/util/Point.h"

//---------------------------------------------------------------------------------------------------------------------

namespace atlas {
namespace grid {

//---------------------------------------------------------------------------------------------------------------------

class IteratorXY {
    using implementation_t = detail::grid::Grid::IteratorXY;

public:
    using difference_type   = implementation_t::difference_type;
    using iterator_category = implementation_t::iterator_category;
    using value_type        = implementation_t::value_type;
    using pointer           = implementation_t::pointer;
    using reference         = implementation_t::reference;

public:
    IteratorXY(std::unique_ptr<implementation_t> iterator): iterator_(std::move(iterator)) {}

    IteratorXY(const IteratorXY& iterator): iterator_(iterator.iterator_->clone()) {}

    bool next(value_type& xy) { return iterator_->next(xy); }

    reference operator*() const { return iterator_->operator*(); }

    const IteratorXY& operator++() {
        iterator_->operator++();
        return *this;
    }

    const IteratorXY& operator+=(difference_type distance) {
        iterator_->operator+=(distance);
        return *this;
    }

    friend IteratorXY operator+(const IteratorXY& a, difference_type distance) {
        IteratorXY result(a);
        result += distance;
        return result;
    }

    friend difference_type operator-(const IteratorXY& last, const IteratorXY& first) {
        return first.iterator_->distance(*last.iterator_);
    }

    bool operator==(const IteratorXY& other) const { return iterator_->operator==(*other.iterator_); }
    bool operator!=(const IteratorXY& other) const { return iterator_->operator!=(*other.iterator_); }

private:
    IteratorXY::difference_type distance(const IteratorXY& other) const {
        return iterator_->distance(*other.iterator_);
    }

private:
    std::unique_ptr<implementation_t> iterator_;
};


//---------------------------------------------------------------------------------------------------------------------

class IteratorLonLat {
    using implementation_t = detail::grid::Grid::IteratorLonLat;

public:
    using difference_type   = implementation_t::difference_type;
    using iterator_category = implementation_t::iterator_category;
    using value_type        = implementation_t::value_type;
    using pointer           = implementation_t::pointer;
    using reference         = implementation_t::reference;

public:
    IteratorLonLat(std::unique_ptr<implementation_t> iterator): iterator_(std::move(iterator)) {}

    IteratorLonLat(const IteratorLonLat& iterator): iterator_(iterator.iterator_->clone()) {}

    bool next(value_type& xy) { return iterator_->next(xy); }

    reference operator*() const { return iterator_->operator*(); }

    const IteratorLonLat& operator++() {
        iterator_->operator++();
        return *this;
    }

    const IteratorLonLat& operator+=(difference_type distance) {
        iterator_->operator+=(distance);
        return *this;
    }

    friend IteratorLonLat operator+(const IteratorLonLat& a, difference_type distance) {
        IteratorLonLat result(a);
        result += distance;
        return result;
    }

    friend difference_type operator-(const IteratorLonLat& last, const IteratorLonLat& first) {
        return first.iterator_->distance(*last.iterator_);
    }

    bool operator==(const IteratorLonLat& other) const { return iterator_->operator==(*other.iterator_); }
    bool operator!=(const IteratorLonLat& other) const { return iterator_->operator!=(*other.iterator_); }

private:
    difference_type distance(const IteratorLonLat& other) const { return iterator_->distance(*other.iterator_); }

private:
    std::unique_ptr<implementation_t> iterator_;
};

//---------------------------------------------------------------------------------------------------------------------

class IterateXY {
public:
    using iterator       = grid::IteratorXY;
    using const_iterator = iterator;
    using Grid           = detail::grid::Grid;

public:
    IterateXY(const Grid& grid): grid_(grid) {}
    iterator begin() const;
    iterator end() const;
    void access(size_t i, PointXY&);
    PointXY front() { return *begin(); }
    PointXY back() { return *(begin() + (grid_.size() - 1)); }

private:
    const Grid& grid_;
};

class IterateLonLat {
public:
    using iterator       = IteratorLonLat;
    using const_iterator = iterator;
    using Grid           = detail::grid::Grid;

public:
    IterateLonLat(const Grid& grid): grid_(grid) {}
    iterator begin() const;
    iterator end() const;
    PointLonLat front() { return *begin(); }
    PointLonLat back() { return *(begin() + (grid_.size() - 1)); }

private:
    const Grid& grid_;
};

}  // namespace grid
}  // namespace atlas

#ifndef DOXYGEN_SHOULD_SKIP_THIS
namespace std {

inline atlas::grid::IteratorXY::difference_type distance(const atlas::grid::IteratorXY& first,
                                                         const atlas::grid::IteratorXY& last) {
    return last - first;
}
inline atlas::grid::IteratorLonLat::difference_type distance(const atlas::grid::IteratorLonLat& first,
                                                             const atlas::grid::IteratorLonLat& last) {
    return last - first;
}

}  // namespace std
#endif
