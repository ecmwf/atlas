#include "atlas/grid/spacing/GaussianSpacing.h"

#include "atlas/grid/spacing/gaussian/Latitudes.h"

namespace atlas {
namespace grid {

GaussianSpacing::GaussianSpacing(const eckit::Parametrisation& params) {
	// no parameters for gaussian spacing!
};

std::vector<double> GaussianSpacing::generate(double xmin, double xmax, size_t N) {

	// gaussian spacing only exists over range (-90,90)
	ASSERT ( xmin==-90.0 );
	ASSERT ( xmax== 90.0 );
	ASSERT ( N%2 == 0 );
	
	// create object containing gaussian latitudes
	double * lats=new double[N];
	spacing::gaussian::gaussian_latitudes_npole_spole(N/2, lats);
	
	// create vector of latitudes
	std::vector<double> x(N);
	
	for (int i=0;i<N;i++) {
		x[i]=lats[i];
	}
	
	delete[] lats;
	
	return x;
};

register_BuilderT1(Spacing,GaussianSpacing,GaussianSpacing::className());

}  // namespace grid
}  // namespace atlas

