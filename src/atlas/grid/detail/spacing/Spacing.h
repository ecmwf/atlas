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

#include <array>
#include <string>
#include <vector>

#include "atlas/util/Object.h"

namespace eckit {
class Parametrisation;
}

namespace atlas {
namespace util {
class Config;
}
}  // namespace atlas

namespace atlas {
namespace grid {
namespace spacing {

class Spacing : public util::Object {
public:
    using const_iterator = std::vector<double>::const_iterator;
    using Interval       = std::array<double, 2>;
    using Spec           = atlas::util::Config;

public:
    static const Spacing* create(const eckit::Parametrisation& params);

    virtual std::string type() const = 0;

    double operator[](size_t i) const { return x_[i]; }

    size_t size() const { return x_.size(); }

    const double* data() const { return x_.data(); }

    const_iterator begin() const { return x_.begin(); }
    const_iterator end() const { return x_.end(); }

    double front() const { return x_.front(); }
    double back() const { return x_.back(); }

    Interval interval() const { return {{min_, max_}}; }

    double min() const { return min_; }
    double max() const { return max_; }

    virtual Spec spec() const = 0;

protected:
    std::vector<double> x_;
    double min_;
    double max_;
};

}  // namespace spacing
}  // namespace grid
}  // namespace atlas
