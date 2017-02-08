#include "atlas/grid/spacing/UniformSpacing.h"

namespace atlas {
namespace grid {
namespace spacing {

UniformSpacing::UniformSpacing(const eckit::Parametrisation& params) {
  // general setup
  Spacing::setup(params);

};

void UniformSpacing::generate(size_t i, double &x) const {
  // points are equally spaced between xmin and xmax

  ASSERT( i<N_ );

  if (N_==1) {
    x=0.5*(xmin_+xmax_);
  } else {
    x=xmin_+i*(xmax_-xmin_)/(N_-1);
  }
};

register_BuilderT1(Spacing,UniformSpacing,UniformSpacing::spacing_type_str());

}  // namespace spacing
}  // namespace grid
}  // namespace atlas

