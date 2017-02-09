#include "atlas/grid/spacing/UniformSpacing.h"
#include "eckit/config/Parametrisation.h"
#include "eckit/exception/Exceptions.h"

namespace atlas {
namespace grid {
namespace spacing {

UniformSpacing::UniformSpacing(const eckit::Parametrisation& params) {

  double xmin;
  double xmax;
  long   N;

  // retrieve xmin, xmax and N from params
  if ( !params.get("xmin",xmin) ) throw eckit::BadParameter("xmin missing in Params",Here());
  if ( !params.get("xmax",xmax) ) throw eckit::BadParameter("xmax missing in Params",Here());
  if ( !params.get("N",   N   ) ) throw eckit::BadParameter("N missing in Params",   Here());

  // points are equally spaced between xmin and xmax
  x_.resize(N);
  if (N==1) {
    x_[0] = 0.5*( xmin + xmax );
  } else {
    const double dx = (xmax-xmin)/double(N-1);
    for( long i=0; i<N; ++i ) {
      x_[i] = xmin + i*dx;
    }
  }
}

register_BuilderT1(Spacing,UniformSpacing,UniformSpacing::spacing_type_str());

}  // namespace spacing
}  // namespace grid
}  // namespace atlas

