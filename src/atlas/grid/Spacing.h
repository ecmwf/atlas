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

#include "atlas/library/config.h"
#include "atlas/util/ObjectHandle.h"

//---------------------------------------------------------------------------------------------------------------------
#ifndef DOXYGEN_SHOULD_SKIP_THIS
// Forward declarations
namespace eckit {
class Parametrisation;
}
namespace atlas {
namespace util {
class Config;
}
namespace grid {
namespace spacing {
class Spacing;
}
}  // namespace grid
}  // namespace atlas
#endif
//---------------------------------------------------------------------------------------------------------------------

namespace atlas {
namespace grid {

//---------------------------------------------------------------------------------------------------------------------

class Spacing : DOXYGEN_HIDE(public util::ObjectHandle<atlas::grid::spacing::Spacing>) {
public:
    using const_iterator = std::vector<double>::const_iterator;
    using Interval       = std::array<double, 2>;
    using Spec           = atlas::util::Config;

public:
    using Handle::Handle;
    Spacing() = default;
    Spacing(const eckit::Parametrisation&);

    size_t size() const;

    double operator[](size_t i) const;

    const_iterator begin() const;
    const_iterator end() const;

    double front() const;
    double back() const;

    Interval interval() const;

    double min() const;
    double max() const;

    std::string type() const;

    Spec spec() const;
};

//---------------------------------------------------------------------------------------------------------------------

class LinearSpacing : public Spacing {
public:
    using Interval = std::array<double, 2>;

public:
    using Spacing::Spacing;
    LinearSpacing() = default;
    LinearSpacing(double start, double stop, long N, bool endpoint = true);
    LinearSpacing(const Interval&, long N, bool endpoint = true);
    double step() const;
};

//---------------------------------------------------------------------------------------------------------------------

class GaussianSpacing : public Spacing {
public:
    using Spacing::Spacing;
    GaussianSpacing() = default;
    GaussianSpacing(long N);
};

//---------------------------------------------------------------------------------------------------------------------

}  // namespace grid
}  // namespace atlas
