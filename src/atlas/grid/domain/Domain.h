/*

The Domain class describes the extent of a grid in projected "grid coordinates"

daand:
	- I simply removed the original Domain.h, which only described boxes in (lon,lat)-space.
	- The Domain class has become a purely abstract class to allow for other domain shapes (circular, frame, and what not...)
	- I didn't implement hashes, (copy) constructors, comparators, etc. for now.



 */

#ifndef atlas_grid_Domain_h
#define atlas_grid_Domain_h

#include <iostream>
#include "eckit/geometry/Point2.h"
#include "eckit/config/Parametrisation.h"
#include "eckit/memory/Builder.h"
#include "atlas/util/Config.h"

namespace atlas {
namespace grid {
namespace domain {

class Domain {

	public:
		typedef const eckit::Parametrisation& ARG1;
		typedef eckit::BuilderT1<Domain> builder_t;

  public:
  
 		static Domain* create() {
			// default: global domain
			util::Config projParams;
			projParams.set("domainType","atlas.GlobalDomain");
			return Domain::create(projParams);
		};
		
		static Domain* create(const eckit::Parametrisation& p) {
		
			std::string domainType;
			if (p.get("domainType",domainType)) {
				return eckit::Factory<Domain>::instance().get(domainType).create(p);
			}

			// should return error here
	    throw eckit::BadParameter("domainType missing in Params",Here());
			return NULL;
		}
		
		// className
		static std::string className() {return "atlas.Domain";}
		
		static std::string domain_type_str() {return "abstract";}
		
		static Domain makeGlobal() {return *create(); };
  	
  	/// Checks if the point is contained in the domain
    virtual bool contains(eckit::geometry::Point2) const {return false; }; // should become purely virtual!
    virtual bool contains(double x, double y) const {return contains(eckit::geometry::Point2(x,y)); }; // should become purely virtual!
    

    /// Output to stream
    void print(std::ostream&) const;
  
  	
		// Here's stuff that's probably needed, but which depends on the projection, so it should be moved to the Grid class.
     
	  /// Check if grid includes the North pole
	  bool includesPoleNorth() const {};

	  /// Check if grid includes the South pole
	  bool includesPoleSouth() const {};

	  /// Check if grid spans the complete range East-West (periodic)
	  bool isPeriodicEastWest() const {};
	  
	  /// Check if domain represents the complete globe surface
		bool isGlobal() const {};

	  /// Check if domain does not represent any area on the globe surface
	  bool isEmpty() const {};
	  
		// dummies for now
    double north() const { return 0.; }
    double west() const { return 0.; }
    double south() const { return 0.; }
    double east() const { return 0.; }    
  
};


}  // namespace domain
}  // namespace grid
}  // namespace atlas


#endif
