#ifndef atlas_grid_CircularDomain_h
#define atlas_grid_CircularDomain_h

#include <iostream>
#include "atlas/grid/domain/Domain.h"
#include "eckit/geometry/Point2.h"

namespace atlas {
namespace grid {
namespace domain {

class CircularDomain: public Domain {

  public:
 		CircularDomain(const eckit::Parametrisation& p);
    ~CircularDomain() {};

    /// Checks if the point is contained in the domain
    bool contains(eckit::geometry::Point2 P) const;
    
    static std::string domain_type_str() {return "circular";}
    
    bool isEmpty() const { return (radius_>0); }
    bool isGlobal() const { return false; }
    
	private:
		
		double xc_, yc_, radius_;
		void setup();
};


}  // namespace domain
}  // namespace grid
}  // namespace atlas


#endif
