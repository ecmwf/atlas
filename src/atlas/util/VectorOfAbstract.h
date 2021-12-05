/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

/// @file Polygon.h
/// @author Pedro Maciel
/// @author Willem Deconinck
/// @date September 2017

#pragma once

#include <memory>
#include <vector>

#include "atlas/library/config.h"

namespace atlas {
namespace util {

//------------------------------------------------------------------------------------------------------

template <typename BaseIterator>
class DereferenceIterator : DOXYGEN_HIDE(public BaseIterator) {
public:
    using value_type = typename BaseIterator::value_type::element_type;
    using pointer    = value_type*;
    using reference  = value_type&;

    DereferenceIterator(const BaseIterator& other): BaseIterator(other) {}

    reference operator*() const { return *(this->BaseIterator::operator*()); }
    pointer operator->() const { return this->BaseIterator::operator*().get(); }
    reference operator[](size_t n) const { return *(this->BaseIterator::operator[](n)); }
};

//------------------------------------------------------------------------------------------------------

template <typename Iterator>
DereferenceIterator<Iterator> make_dereference_iterator(Iterator t) {
    return DereferenceIterator<Iterator>(t);
}

//------------------------------------------------------------------------------------------------------

template <typename Abstract>
class VectorOfAbstract {
public:
    using value_type      = Abstract;
    using container_type  = std::vector<std::unique_ptr<value_type> >;
    using const_reference = const value_type&;
    using reference       = const_reference;
    using const_iterator  = DereferenceIterator<typename container_type::const_iterator>;

public:
    VectorOfAbstract() = default;
    VectorOfAbstract(VectorOfAbstract&& other): container_(std::move(other.container_)) {}

    const_iterator begin() const { return make_dereference_iterator(container_.begin()); }
    const_iterator end() const { return make_dereference_iterator(container_.end()); }
    const_reference operator[](idx_t i) const { return *container_[i]; }
    const_reference at(idx_t i) const { return *container_[i]; }
    idx_t size() const { return container_.size(); }
    void reserve(size_t size) { container_.reserve(size); }
    template <typename... Args>
    void emplace_back(Args&&... args) {
        container_.emplace_back(std::forward<Args>(args)...);
    }
    container_type& get() { return container_; }
    void clear() { return container_.clear(); }

private:
    container_type container_;
};

//------------------------------------------------------------------------------------------------------

}  // namespace util
}  // namespace atlas
