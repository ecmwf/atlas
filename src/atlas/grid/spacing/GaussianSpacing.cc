#include "atlas/grid/spacing/GaussianSpacing.h"
#include "atlas/grid/spacing/gaussian/Latitudes.h"
#include "eckit/config/Parametrisation.h"
#include "eckit/exception/Exceptions.h"

namespace atlas {
namespace grid {
namespace spacing {

GaussianSpacing::GaussianSpacing(long nlat) {
  // perform checks
  ASSERT ( nlat%2 == 0 );

  // initialize latitudes during setup, to avoid repeating it.
  x_.resize(nlat);
  gaussian::gaussian_latitudes_npole_spole(nlat/2, x_.data());
}

GaussianSpacing::GaussianSpacing(const eckit::Parametrisation& params) {

  // retrieve N from params
  long N;
  if ( !params.get("N",N) )      throw eckit::BadParameter("N missing in Params",Here());

  // perform checks
  ASSERT ( N%2 == 0 );

  // initialize latitudes during setup, to avoid repeating it.
  x_.resize(N);
  gaussian::gaussian_latitudes_npole_spole(N/2, x_.data());

  // Not yet implemented: specify different bounds or direction (e.g from south to north pole)
  double start =  90.;
  double end   = -90.;
  params.get("start", start);
  params.get("end",   end  );
  if( start!=90. && end!=-90. ) {
    NOTIMP;
  }
  
  min_ = std::min(start,end);
  max_ = std::max(start,end);

}

register_BuilderT1(Spacing,GaussianSpacing,GaussianSpacing::spacing_type_str());

}  // namespace spacing
}  // namespace grid
}  // namespace atlas

