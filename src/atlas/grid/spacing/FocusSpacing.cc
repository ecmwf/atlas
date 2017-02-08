#include "atlas/grid/spacing/FocusSpacing.h"
#include <cmath>

namespace atlas {
namespace grid {
namespace spacing {

FocusSpacing::FocusSpacing(const eckit::Parametrisation& params) {
  // general setup
  Spacing::setup(params);

  // additional parameters for focus spacing
  if( ! params.get("focus_factor",focus_factor_) )
    throw eckit::BadParameter("focus_factor missing in Params",Here());
}

void FocusSpacing::generate(size_t i, double & x) const {

  ASSERT( i<N_ );

  // boundary cases
  if (i==0) { x=xmin_; return; }
  if (i==N_-1) { x=xmax_; return; }

  // inside the range: points are unevenly spaced between xmin and xmax
  double xx=(2*i-int(N_-1))/double(N_-1);    // between -1 and 1;
  x=(xmin_+xmax_)/2+(xmax_-xmin_)/2*atan(tan(M_PI*xx/2)/focus_factor_)*2*M_1_PI;
  return;
}

register_BuilderT1(Spacing,FocusSpacing,FocusSpacing::spacing_type_str());

}  // namespace spacing
}  // namespace grid
}  // namespace atlas

