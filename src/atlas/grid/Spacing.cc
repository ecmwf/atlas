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

namespace atlas {
namespace grid {

Spacing::Spacing() : spacing_( nullptr ) {}

Spacing::Spacing( const Spacing& other ) : spacing_( other.spacing_ ) {}

Spacing::Spacing( const spacing::Spacing* spacing ) : spacing_( spacing ) {}

Spacing::Spacing( const eckit::Parametrisation& p ) : spacing_( atlas::grid::spacing::Spacing::create( p ) ) {}

LinearSpacing::LinearSpacing( double start, double stop, long N, bool endpoint ) :
    Spacing( new atlas::grid::spacing::LinearSpacing( start, stop, N, endpoint ) ) {}

LinearSpacing::LinearSpacing( const Interval& interval, long N, bool endpoint ) :
    Spacing( new atlas::grid::spacing::LinearSpacing( interval, N, endpoint ) ) {}

GaussianSpacing::GaussianSpacing( long N ) : Spacing( new atlas::grid::spacing::GaussianSpacing( N ) ) {}

}  // namespace grid
}  // namespace atlas
