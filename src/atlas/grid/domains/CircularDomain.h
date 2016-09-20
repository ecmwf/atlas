#ifndef atlas_grid_CircularDomain_h
#define atlas_grid_CircularDomain_h

#include <iostream>
#include "atlas/grid/domains/Domain.h"
#include "eckit/geometry/Point2.h"

namespace atlas {
namespace grid {


class CircularDomain: public Domain {

  public:
 
    CircularDomain(double xc, double yc, double r);				// from center coordinates and radius
    CircularDomain(eckit::geometry::Point2 XC, double r);		// from center point and radius
    ~CircularDomain() {};

    /// Checks if the point is contained in the domain
    bool contains(eckit::geometry::Point2 P);
    
	private:
		
		double xc_, yc_, radius_;
		void setup();
};


}  // namespace grid
}  // namespace atlas


#endif
