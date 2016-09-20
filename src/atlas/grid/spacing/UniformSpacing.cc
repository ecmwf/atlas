#include "atlas/grid/spacing/UniformSpacing.h"

namespace atlas {
namespace grid {

UniformSpacing::UniformSpacing(const eckit::Parametrisation& params) {
	// no parameters for uniform spacing!
};

std::vector<double> UniformSpacing::generate(double xmin, double xmax, size_t N) {
	// create vector of points, equally spaced between xmin and xmax
	std::vector<double> x(N);
	for (int i=0;i<N;i++) {
		x[i]=xmin+i*(xmax-xmin)/(N-1);
	}
	
	return x;
};

register_BuilderT1(Spacing,UniformSpacing,UniformSpacing::className());

}  // namespace grid
}  // namespace atlas

