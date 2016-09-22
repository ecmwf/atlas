#include "atlas/grid/spacing/GaussianSpacing.h"

#include "atlas/grid/spacing/gaussian/Latitudes.h"

namespace atlas {
namespace grid {
namespace spacing {

GaussianSpacing::GaussianSpacing(const eckit::Parametrisation& params) {
	// no parameters for gaussian spacing!
};

void GaussianSpacing::generate(double xmin, double xmax, size_t N, std::vector<double> &x) const {

	// gaussian spacing only exists over range (-90,90)
	ASSERT ( xmin==-90.0 );
	ASSERT ( xmax== 90.0 );
	ASSERT ( N%2 == 0 );
	
	// create object containing gaussian latitudes
	double * lats=new double[N];
	spacing::gaussian::gaussian_latitudes_npole_spole(N/2, lats);
	
	// move to x
	x.assign(lats,lats+N);
	
	// clean up
	delete[] lats;
	
};

register_BuilderT1(Spacing,GaussianSpacing,GaussianSpacing::spacing_type_str());

}  // namespace spacing
}  // namespace grid
}  // namespace atlas

