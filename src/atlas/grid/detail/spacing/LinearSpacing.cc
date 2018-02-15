#include "atlas/grid/detail/spacing/LinearSpacing.h"
#include <cmath>
#include "eckit/config/Parametrisation.h"
#include "eckit/exception/Exceptions.h"

namespace atlas {
namespace grid {
namespace spacing {

LinearSpacing::Params::Params( const eckit::Parametrisation& params ) {
    endpoint = true;
    params.get( "endpoint", endpoint );
    std::vector<double> interval;
    // if( params.get("step",step) ) {
    //
    //   // Several combinations possible:
    //   if( params.get("start",start) && params.get("end",end) ) {
    //     N = long( (end-start)/step ) + (endpoint ? 1 : 0 );
    //   } else if( params.get("centre",centre) && params.get("N",N) ) {
    //     start = endpoint ? centre - step * 0.5*double(N-1)
    //                      : centre - step * 0.5*double(N);
    //     end   = endpoint ? start + step * double(N-1) :
    //                        start + step * double(N);
    //   } else {
    //     throw eckit::BadParameter("Invalid combination of parameters",Here());
    //   }
    // }
    // else
    if ( params.get( "N", N ) ) {
        // Only one remaining combinations possible:
        if ( params.get( "start", start ) && params.get( "end", end ) ) {
            // OK
        }
        else if ( params.get( "interval", interval ) ) {
            start = interval[0];
            end   = interval[1];
        }
        else if ( params.get( "start", start ) && params.get( "length", length ) ) {
            end = start + length;
        }
        else {
            throw eckit::BadParameter( "Invalid combination of parameters", Here() );
        }
    }
    else {
        throw eckit::BadParameter( "Invalid combination of parameters", Here() );
    }
    length = end - start;

    if ( endpoint && N > 1 )
        step = length / double( N - 1 );
    else
        step = length / double( N );
}

LinearSpacing::LinearSpacing( const eckit::Parametrisation& params ) {
    Params p( params );
    setup( p.start, p.end, p.N, p.endpoint );
}

// LinearSpacing::LinearSpacing( double centre, double step, long N, bool
// endpoint ) {
//   double start = endpoint ? centre - step * double(N-1)/2. :
//                             centre - step * double(N)/2.;
//   double end   = endpoint ? start + step * double(N-1) :
//                             start + step * double(N);
//   setup(start,end,N,endpoint);
// }

LinearSpacing::LinearSpacing( double start, double end, long N, bool endpoint ) {
    setup( start, end, N, endpoint );
}

LinearSpacing::LinearSpacing( const Interval& interval, long N, bool endpoint ) {
    setup( interval[0], interval[1], N, endpoint );
}

void LinearSpacing::setup( double start, double end, long N, bool endpoint ) {
    x_.resize( N );

    double step;
    if ( endpoint && N > 1 )
        step = ( end - start ) / double( N - 1 );
    else
        step = ( end - start ) / double( N );

    for ( long i = 0; i < N; ++i ) {
        x_[i] = start + i * step;
    }

    min_ = std::min( start, end );
    max_ = std::max( start, end );

    start_    = start;
    end_      = end;
    N_        = N;
    endpoint_ = endpoint;
}

double LinearSpacing::step() const {
    if ( size() > 1 )
        return x_[1] - x_[0];
    else
        return 0.;
}

bool LinearSpacing::endpoint() const {
    return std::abs( x_.back() - max_ ) < 1.e-12;
}

LinearSpacing::Spec LinearSpacing::spec() const {
    Spec spacing_specs;
    spacing_specs.set( "type", static_type() );
    spacing_specs.set( "start", start_ );
    spacing_specs.set( "end", end_ );
    spacing_specs.set( "N", N_ );
    spacing_specs.set( "endpoint", endpoint_ );
    return spacing_specs;
}

register_BuilderT1( Spacing, LinearSpacing, LinearSpacing::static_type() );

}  // namespace spacing
}  // namespace grid
}  // namespace atlas
