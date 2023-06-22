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

#include <cassert>
#include <vector>

#include "atlas/library/config.h"
#include "atlas/util/Config.h"

namespace atlas {

//---------------------------------------------------------------------------------------------------------------------

class Vertical {
public:
    template <typename vector_t>  // expect "vector_t::size()" and "vector_t::operator[]"
    Vertical(idx_t levels, const vector_t& z, const util::Config& config = util::NoConfig());

    template <typename vector_t, typename Interval>  // expect "vector_t::size()" and "vector_t::operator[]"
    Vertical(idx_t levels, const vector_t& z, const Interval& interval, const util::Config& config = util::NoConfig());

    Vertical(const util::Config& config = util::NoConfig());

public:
    idx_t k_begin() const { return k_begin_; }
    idx_t k_end() const { return k_end_; }
    idx_t size() const { return size_; }

    template <typename Int>
    double operator()(const Int k) const {
        return z_[k];
    }

    template <typename Int>
    double operator[](const Int k) const {
        return z_[k];
    }

    double min() const { return min_; }
    double max() const { return max_; }

    double front() const { return z_.front(); }
    double back() const { return z_.back(); }

    /// @brief Output information of field
    friend std::ostream& operator<<(std::ostream& os, const Vertical& v);

private:
    idx_t size_;
    idx_t k_begin_;
    idx_t k_end_;
    std::vector<double> z_;
    double min_;
    double max_;
};

//---------------------------------------------------------------------------------------------------------------------

template <typename vector_t, typename Interval>
Vertical::Vertical(idx_t levels, const vector_t& z, const Interval& interval, const util::Config& config):
    Vertical(levels, z, config) {
    min_ = interval[0];
    max_ = interval[1];
}

//---------------------------------------------------------------------------------------------------------------------

template <typename vector_t>
Vertical::Vertical(idx_t levels, const vector_t& z, const util::Config&): size_{levels}, k_begin_{0}, k_end_{size_} {
    assert(size_ == static_cast<idx_t>(z.size()));
    z_.resize(size_);
    for (idx_t k = 0; k < size_; ++k) {
        z_[k] = z[k];
    }
    min_ = (size_ ? z[0] : 0.);
    max_ = (size_ ? z[size_ - 1] : 1.);
}

//---------------------------------------------------------------------------------------------------------------------

}  // namespace atlas
