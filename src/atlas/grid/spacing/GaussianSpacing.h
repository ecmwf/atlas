#ifndef atlas_GaussianSpacing_H
#define atlas_GaussianSpacing_H

#include "atlas/grid/spacing/Spacing.h"

namespace atlas {
namespace grid {
namespace spacing {

class GaussianSpacing: public Spacing {


	public:
		
		// constructor
		GaussianSpacing(const eckit::Parametrisation& p);
		~GaussianSpacing();	// required to deallocate lats_
		
		// class name
		static std::string className() { return "atlas.GaussianSpacing"; } 
		static std::string spacing_type_str() {return "gaussian";}

		
		void generate(size_t i, double &x) const;
	
	private:
		double * lats_;
};

}  // namespace spacing
}  // namespace grid
}  // namespace atlas


#endif
