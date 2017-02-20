#include <algorithm>
#include "atlas/grid/spacing/CustomSpacing.h"
#include "eckit/config/Parametrisation.h"
#include "eckit/exception/Exceptions.h"

namespace atlas {
namespace grid {
namespace spacing {

CustomSpacing::CustomSpacing(long N, const double values[], const std::array<double,2>& span) {

  x_.assign(values,values+N);
  min_ = std::min(span[0],span[1]);
  max_ = std::max(span[0],span[1]);
}

CustomSpacing::CustomSpacing(const eckit::Parametrisation& params) {

  params.get("values",x_);

  long N;
  if( params.get("N",N) ) {
    ASSERT( x_.size() == N );
  }
  N = x_.size();
  
  std::vector<double> span;
  if( params.get("span",span) ) {

    min_ = std::min(span[0],span[1]);
    max_ = std::max(span[0],span[1]);

  } else {
  
    min_ = x_.front();
    max_ = x_.front();
    for( size_t j=1; j<N; ++j ) {
      min_ = std::min(min_, x_[j]);
      max_ = std::max(max_, x_[j]); 
    }
  }
}

register_BuilderT1(Spacing,CustomSpacing,CustomSpacing::spacing_type_str());

}  // namespace spacing
}  // namespace grid
}  // namespace atlas

