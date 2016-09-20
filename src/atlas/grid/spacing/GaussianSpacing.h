#ifndef atlas_GaussianSpacing_H
#define atlas_GaussianSpacing_H

#include "atlas/grid/spacing/Spacing.h"

namespace atlas {
namespace grid {

class GaussianSpacing: public Spacing {


	public:
		
		// constructor
		GaussianSpacing(const eckit::Parametrisation& p);
		
		// class name
		static std::string className() { return "atlas.GaussianSpacing"; } 

		
		std::vector<double> generate(double xmin, double xmax, size_t N);
};

}  // namespace grid
}  // namespace atlas


#endif
