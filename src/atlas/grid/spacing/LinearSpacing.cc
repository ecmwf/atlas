#include "atlas/grid/spacing/LinearSpacing.h"
#include "eckit/config/Parametrisation.h"
#include "eckit/exception/Exceptions.h"

namespace atlas {
namespace grid {
namespace spacing {

LinearSpacing::LinearSpacing(const eckit::Parametrisation& params) {

  double step;
  double centre;
  double start;
  double end;
  double length;
  long   N;
  bool   endpoint=true;

  params.get("endpoint", endpoint );

  if( params.get("step",step) ) {
    
    // Several combinations possible:
    if( params.get("start",start) && params.get("end",end) ) {
      N = long( (end-start)/step ) + (endpoint ? 1 : 0 );
    } else if( params.get("centre",centre) && params.get("N",N) ) {
      start = endpoint ? centre - step * 0.5*double(N-1)
                       : centre - step * 0.5*double(N);
      end   = endpoint ? start + step * double(N-1) :
                         start + step * double(N);
    } else {
      throw eckit::BadParameter("Invalid combination of parameters",Here());
    }
  }
  else if( params.get("N",N) ) {
    // Only one remaining combinations possible:
    if( params.get("start",start) && params.get("end",end) ) {
      // OK
    } else if( params.get("start",start) && params.get("length",length) ) {
      end = start + length;
    } else {
      throw eckit::BadParameter("Invalid combination of parameters",Here());
    }
  }
  else {
    throw eckit::BadParameter("Invalid combination of parameters",Here());
  }

  setup(start,end,N,endpoint);
}

LinearSpacing::LinearSpacing( double centre, double step, long N, bool endpoint ) {
  double start = endpoint ? centre - step * double(N-1)/2. :
                            centre - step * double(N)/2.;
  double end   = endpoint ? start + step * double(N-1) :
                            start + step * double(N);
  setup(start,end,N,endpoint);
}

LinearSpacing::LinearSpacing( std::array<double,2> range, long N, bool endpoint ) {
  double start = range[0];
  double end   = range[1];
  setup(start,end,N,endpoint);
}

void LinearSpacing::setup(double start, double end, long N, bool endpoint) {

  x_.resize(N);

  double step;
  if( endpoint && N>1 )
    step = (end-start)/double(N-1);
  else
    step = (end-start)/double(N);

  for( long i=0; i<N; ++i ) {
    x_[i] = start + i*step;
  }

  min_ = std::min(start,end);
  max_ = std::max(start,end);
}

register_BuilderT1(Spacing,LinearSpacing,LinearSpacing::spacing_type_str());

}  // namespace spacing
}  // namespace grid
}  // namespace atlas

