#include "atlas/grid/spacing/UniformSpacing.h"
#include "eckit/config/Parametrisation.h"
#include "eckit/exception/Exceptions.h"

#include "atlas/internals/Debug.h"

namespace atlas {
namespace grid {
namespace spacing {

UniformSpacing::UniformSpacing(const eckit::Parametrisation& params) {

  double step;
  double centre;
  double min;
  double max;
  long   N;
  bool   endpoint=true;

  params.get("endpoint", endpoint );

  if( params.get("step",step) ) {
    // Several combinations possible:
    if( params.get("min",min) && params.get("max",max) ) {
      N = long( (max-min)/step ) + (endpoint ? 1 : 0 );
    } else if( params.get("centre",centre) && params.get("N",N) ) {
      min = endpoint ? centre - step * 0.5*double(N-1)
                     : centre - step * 0.5*double(N);
      max = endpoint ? min + step * double(N-1) :
                       min + step * double(N);
    } else {
      throw eckit::BadParameter("Invalid combination of parameters",Here());
    }
  }
  else if( params.get("N",N) ) {
    // Only one remaining combinations possible:
    if( params.get("min",min) && params.get("max",max) ) {
      // OK
    } else {
      throw eckit::BadParameter("Invalid combination of parameters",Here());
    }
  }
  else {
    throw eckit::BadParameter("Invalid combination of parameters",Here());
  }

  setup(min,max,N,endpoint);
}

UniformSpacing::UniformSpacing( double centre, double step, long N, bool endpoint ) {
  double min = endpoint ? centre - step * double(N-1)/2. :
                          centre - step * double(N)/2.;
  double max = endpoint ? min + step * double(N-1) :
                          min + step * double(N);
  setup(min,max,N,endpoint);
}

UniformSpacing::UniformSpacing( std::array<double,2> min_max, long N, bool endpoint ) {
  double min = min_max[0];
  double max = min_max[1];
  setup(min,max,N,endpoint);
}

void UniformSpacing::setup(double xmin, double xmax, long N, bool endpoint) {

  x_.resize(N);

  double dx;
  if( endpoint && N>1 )
    dx = (xmax-xmin)/double(N-1);
  else
    dx = (xmax-xmin)/double(N);

  for( long i=0; i<N; ++i ) {
    x_[i] = xmin + i*dx;
  }

}

register_BuilderT1(Spacing,UniformSpacing,UniformSpacing::spacing_type_str());

}  // namespace spacing
}  // namespace grid
}  // namespace atlas

