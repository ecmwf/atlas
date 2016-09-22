#include "atlas/grid/spacing/UniformSpacing.h"

namespace atlas {
namespace grid {
namespace spacing {

UniformSpacing::UniformSpacing(const eckit::Parametrisation& params) {
	// no parameters for uniform spacing!
};

void UniformSpacing::generate(double xmin, double xmax, size_t N, std::vector<double> &x) const {
	// create vector of points, equally spaced between xmin and xmax
	
	if (N==1) {
		x[0]=0.5*(xmin+xmax);
	} else {
		for (int i=0;i<N;i++) {
			x[i]=xmin+i*(xmax-xmin)/(N-1);
		}
	}
};

register_BuilderT1(Spacing,UniformSpacing,UniformSpacing::spacing_type_str());

}  // namespace spacing
}  // namespace grid
}  // namespace atlas

