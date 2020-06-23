/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "atlas/grid/Spacing.h"
#include "atlas/grid/detail/spacing/GaussianSpacing.h"
#include "atlas/grid/detail/spacing/LinearSpacing.h"
#include "atlas/grid/detail/spacing/Spacing.h"
#include "atlas/util/Config.h"

namespace atlas {
namespace grid {

Spacing::Spacing( const eckit::Parametrisation& p ) : Handle( atlas::grid::spacing::Spacing::create( p ) ) {}

size_t Spacing::size() const {
    return get()->size();
}

double Spacing::operator[]( size_t i ) const {
    return get()->operator[]( i );
}

Spacing::const_iterator Spacing::begin() const {
    return get()->begin();
}

Spacing::const_iterator Spacing::end() const {
    return get()->end();
}

double Spacing::front() const {
    return get()->front();
}

double Spacing::back() const {
    return get()->back();
}

Spacing::Interval Spacing::interval() const {
    return get()->interval();
}

double Spacing::min() const {
    return get()->min();
}

double Spacing::max() const {
    return get()->max();
}

std::string Spacing::type() const {
    return get()->type();
}

Spacing::Spec Spacing::spec() const {
    return get()->spec();
}

LinearSpacing::LinearSpacing( double start, double stop, long N, bool endpoint ) :
    Spacing( new atlas::grid::spacing::LinearSpacing( start, stop, N, endpoint ) ) {}

LinearSpacing::LinearSpacing( const Interval& interval, long N, bool endpoint ) :
    Spacing( new atlas::grid::spacing::LinearSpacing( interval, N, endpoint ) ) {}

double LinearSpacing::step() const {
    return dynamic_cast<const spacing::LinearSpacing*>( get() )->step();
}

GaussianSpacing::GaussianSpacing( long N ) : Spacing( new atlas::grid::spacing::GaussianSpacing( N ) ) {}

}  // namespace grid
}  // namespace atlas
