#ifndef atlas_UniformSpacing_H
#define atlas_UniformSpacing_H

#include "atlas/grid/spacing/Spacing.h"

namespace atlas {
namespace grid {

class UniformSpacing: public Spacing {


	public:
		
		// constructor
		UniformSpacing(const eckit::Parametrisation& p);
		
		// class name
		static std::string className() { return "atlas.UniformSpacing"; } 

		
		std::vector<double> generate(double xmin, double xmax, size_t N);
};

}  // namespace grid
}  // namespace atlas


#endif
