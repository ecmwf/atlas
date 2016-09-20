#ifndef atlas_FocusSpacing_H
#define atlas_FocusSpacing_H

#include "atlas/grid/spacing/Spacing.h"

namespace atlas {
namespace grid {

class FocusSpacing: public Spacing {


	public:
		
		// constructor
		FocusSpacing(const eckit::Parametrisation& p);
		
		// class name
		static std::string className() { return "atlas.FocusSpacing"; } 

		
		std::vector<double> generate(double xmin, double xmax, size_t N);
	
	private:
		double focus_factor_;
};

}  // namespace grid
}  // namespace atlas


#endif
