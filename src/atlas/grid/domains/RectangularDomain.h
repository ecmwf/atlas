#ifndef atlas_grid_RectangularDomain_h
#define atlas_grid_RectangularDomain_h

#include <iostream>
#include "atlas/grid/domains/Domain.h"
#include "eckit/geometry/Point2.h"

namespace atlas {
namespace grid {


class RectangularDomain: public Domain {

  public:
 	
 		// constructor
    RectangularDomain(const eckit::Parametrisation& p);				// from 4 values

		// class name
		static std::string className() { return "atlas.RectangularDomain"; }

    /// Checks if the point is contained in the domain
    bool contains(eckit::geometry::Point2 P);
    
	private:
		
		double xmin_, xmax_, ymin_, ymax_;
		void setup();
};


}  // namespace grid
}  // namespace atlas


#endif
