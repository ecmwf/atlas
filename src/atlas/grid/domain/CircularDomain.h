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

    /// Checks if the point is contained in the domain
    virtual bool contains(double x, double y) const;

    static std::string domain_type_str() {return "circular";}
    virtual std::string virtual_domain_type_str() const { return "circular"; }

    virtual bool isEmpty() const { return (radius_>0); }
    virtual bool isGlobal() const { return false; }

    virtual eckit::Properties spec() const;

    virtual double xmin() const { return xc_-radius_; }
    virtual double xmax() const { return xc_+radius_; }
    virtual double ymin() const { return yc_-radius_; }
    virtual double ymax() const { return yc_+radius_; }

    virtual bool isPeriodicX() const { return false; }
    virtual bool isPeriodicY() const { return false; }

private:
    double xc_, yc_, radius_, rr_;
};


}  // namespace domain
}  // namespace grid
}  // namespace atlas


#endif
