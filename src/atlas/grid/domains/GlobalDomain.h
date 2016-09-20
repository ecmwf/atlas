#ifndef atlas_grid_GlobalDomain_h
#define atlas_grid_GlobalDomain_h

#include <iostream>
#include "atlas/grid/domains/Domain.h"
#include "eckit/geometry/Point2.h"

namespace atlas {
namespace grid {


class GlobalDomain: public Domain {

  public:
 
    GlobalDomain();
    ~GlobalDomain() {};

    /// Checks if the point is contained in the domain
    bool contains(eckit::geometry::Point2 P);
    
	private:
		
		void setup();
};


}  // namespace grid
}  // namespace atlas


#endif
