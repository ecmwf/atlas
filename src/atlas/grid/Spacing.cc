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
