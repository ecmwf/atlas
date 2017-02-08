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
#include "eckit/memory/Owned.h"
#include "atlas/util/Config.h"

namespace atlas {
namespace grid {
namespace domain {

class Domain : public eckit::Owned {

public:

    typedef const eckit::Parametrisation& ARG1;
    typedef eckit::BuilderT1<Domain> builder_t;

public:

    static Domain* create() {
      // default: global domain
      util::Config projParams;
      projParams.set("domainType","global");
      return Domain::create(projParams);
    }

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

    static std::string domain_type_str() { return "domain"; }
    virtual std::string virtual_domain_type_str() const { return "domain"; }

    //Domain * makeGlobal() {return create(); };

    /// Checks if the point is contained in the domain
    virtual bool contains(eckit::geometry::Point2) const =0;
    virtual bool contains(double x, double y) const {return contains(eckit::geometry::Point2(x,y)); }

    // Specification of grid
    virtual eckit::Properties spec() const =0;

    /// Output to stream
    void print(std::ostream&) const;

    /// Check if grid includes the North pole
    bool includesPoleNorth() const { return isGlobal(); }

    /// Check if grid includes the South pole
    bool includesPoleSouth() const { return isGlobal(); }

    /// Check if grid spans the complete range East-West (periodic)
    virtual bool isPeriodicEastWest() const { return isGlobal(); }

    /// Check if domain represents the complete globe surface
    virtual bool isGlobal() const =0;

    /// Check if domain does not represent any area on the globe surface
    virtual bool isEmpty() const =0;

};


}  // namespace domain
}  // namespace grid
}  // namespace atlas


#endif
