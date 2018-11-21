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

#include "eckit/memory/SharedPtr.h"

#include "atlas/grid/detail/spacing/Spacing.h"

//---------------------------------------------------------------------------------------------------------------------

// Forward declarations
namespace eckit {
class Parametrisation;
}

//---------------------------------------------------------------------------------------------------------------------

namespace atlas {
namespace grid {

//---------------------------------------------------------------------------------------------------------------------

class Spacing {
public:
    using Implementation = atlas::grid::spacing::Spacing;
    using const_iterator = Implementation::const_iterator;
    using Interval       = Implementation::Interval;
    using Spec           = Implementation::Spec;

public:
    Spacing();
    Spacing( const Spacing& );
    Spacing( const atlas::grid::spacing::Spacing* );
    Spacing( const eckit::Parametrisation& );

    operator bool() const { return spacing_; }

    operator const atlas::grid::spacing::Spacing*() { return spacing_.get(); }

    size_t size() const { return spacing_.get()->size(); }

    double operator[]( size_t i ) const { return spacing_.get()->operator[]( i ); }

    const_iterator begin() const { return spacing_.get()->begin(); }
    const_iterator end() const { return spacing_.get()->end(); }

    double front() const { return spacing_.get()->front(); }
    double back() const { return spacing_.get()->back(); }

    Interval interval() const { return spacing_.get()->interval(); }

    double min() const { return spacing_.get()->min(); }
    double max() const { return spacing_.get()->max(); }

    std::string type() const { return spacing_.get()->type(); }

    Spec spec() const { return spacing_.get()->spec(); }

    const atlas::grid::spacing::Spacing* get() const { return spacing_.get(); }

private:
    eckit::SharedPtr<const atlas::grid::spacing::Spacing> spacing_;
};

//---------------------------------------------------------------------------------------------------------------------

class LinearSpacing : public Spacing {
public:
    using Interval = std::array<double, 2>;

public:
    LinearSpacing( double start, double stop, long N, bool endpoint = true );
    LinearSpacing( const Interval&, long N, bool endpoint = true );
};

//---------------------------------------------------------------------------------------------------------------------

class GaussianSpacing : public Spacing {
public:
    GaussianSpacing( long N );
};

//---------------------------------------------------------------------------------------------------------------------

}  // namespace grid
}  // namespace atlas
