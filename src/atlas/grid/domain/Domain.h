/*

The Domain class describes the extent of a grid in projected "grid coordinates"

daand:
  - I simply removed the original Domain.h, which only described boxes in (lon,lat)-space.
  - The Domain class has become a purely abstract class to allow for other domain shapes (circular, frame, and what not...)
  - I didn't implement hashes, (copy) constructors, comparators, etc. for now.



 */

#ifndef atlas_grid_domain_Domain_h
#define atlas_grid_domain_Domain_h

#include "eckit/geometry/Point2.h"
#include "eckit/config/Parametrisation.h"
#include "eckit/value/Properties.h"
#include "eckit/memory/Owned.h"
#include "eckit/memory/Builder.h"

namespace atlas {
namespace grid {
namespace projection {
class Projection;
}

namespace domain {

class Domain : public eckit::Owned {

public:

    typedef const eckit::Parametrisation& ARG1;
    typedef eckit::BuilderT1<Domain> builder_t;

public:

    static Domain* create(); // Create a global domain

    static Domain* create(const eckit::Parametrisation&);

    // className
    static std::string className() {return "atlas.Domain";}

    static std::string domain_type_str() { return "domain"; }
    virtual std::string virtual_domain_type_str() const { return "domain"; }

    /// Checks if the point is contained in the domain
    bool contains(eckit::geometry::Point2 p) const { return contains(p[0],p[1]); }
    virtual bool contains(double x, double y) const =0;

    // Specification of grid
    virtual eckit::Properties spec() const =0;

    /// Output to stream
    void print(std::ostream&) const;

    /// Check if domain represents the complete globe surface
    virtual bool isGlobal() const =0;

    /// Check if domain does not represent any area on the globe surface
    virtual bool isEmpty() const =0;

// Unless the domain is global, we can never be sure about these functions
// without knowing also the projection

    /// Check if grid includes the North pole
    bool includesNorthPole(const projection::Projection& ) const;

    /// Check if grid includes the South pole
    bool includesSouthPole(const projection::Projection& ) const;

    virtual bool isPeriodicX() const =0;
    virtual bool isPeriodicY() const =0;

    virtual double xmin() const =0;
    virtual double xmax() const =0;
    virtual double ymin() const =0;
    virtual double ymax() const =0;

};


}  // namespace domain
}  // namespace grid
}  // namespace atlas


#endif
