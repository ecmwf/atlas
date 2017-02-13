#ifndef atlas_grid_EmptyDomain_h
#define atlas_grid_EmptyDomain_h

#include <iostream>
#include "atlas/grid/domain/Domain.h"
#include "eckit/geometry/Point2.h"

namespace atlas {
namespace grid {
namespace domain {

class EmptyDomain: public Domain {

public:

    EmptyDomain(const eckit::Parametrisation& p);

    /// Checks if the point is contained in the domain
    virtual bool contains(double x, double y) const { return false; }

    static std::string domain_type_str() { return "empty"; }
    virtual std::string virtual_domain_type_str() const { return "empty"; }

    virtual bool isEmpty() const { return true; }
    virtual bool isGlobal() const { return false; }

    virtual eckit::Properties spec() const;

    virtual double xmin() const { return 0; }
    virtual double xmax() const { return 0; }
    virtual double ymin() const { return 0; }
    virtual double ymax() const { return 0; }

    virtual bool isPeriodicX() const { return false; }
    virtual bool isPeriodicY() const { return false; }

};


}  // namespace domain
}  // namespace grid
}  // namespace atlas


#endif
