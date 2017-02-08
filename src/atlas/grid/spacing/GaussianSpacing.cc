#include "atlas/grid/spacing/GaussianSpacing.h"

#include "atlas/grid/spacing/gaussian/Latitudes.h"

namespace atlas {
namespace grid {
namespace spacing {

GaussianSpacing::GaussianSpacing(const eckit::Parametrisation& params) {
  // general setup
  Spacing::setup(params);

  // perform checks
  ASSERT ( xmin_== 90.0 );
  ASSERT ( xmax_== -90.0 );
  ASSERT ( N_%2 == 0 );

  // initialize latitudes during setup, to avoid repeating it.
  lats_=new double[N_];
  spacing::gaussian::gaussian_latitudes_npole_spole(N_/2, lats_);
};

GaussianSpacing::~GaussianSpacing() {
  // clean up
  delete[] lats_;
}

void GaussianSpacing::generate(size_t i, double &x) const {
  ASSERT( i<N_ );
  x=lats_[i];
};

register_BuilderT1(Spacing,GaussianSpacing,GaussianSpacing::spacing_type_str());

}  // namespace spacing
}  // namespace grid
}  // namespace atlas

