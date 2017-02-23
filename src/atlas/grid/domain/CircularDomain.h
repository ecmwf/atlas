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

    static std::string static_type() {return "circular";}
    virtual std::string type() const { return static_type(); }

    virtual bool isEmpty() const { return (radius_>0); }
    virtual bool isGlobal() const { return false; }

    virtual eckit::Properties spec() const;

    virtual void print(std::ostream&) const;
    
    virtual std::string units() const; // Not implemented

private:
    double xc_, yc_, radius_, rr_;
};


}  // namespace domain
}  // namespace grid
}  // namespace atlas


#endif
