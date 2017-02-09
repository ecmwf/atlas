#ifndef atlas_grid_RectangularDomain_h
#define atlas_grid_RectangularDomain_h

#include <iostream>
#include "atlas/grid/domain/Domain.h"
#include "eckit/geometry/Point2.h"

namespace atlas {
namespace grid {
namespace domain {

class RectangularDomain: public Domain {

public:

     // constructor
    RectangularDomain(const eckit::Parametrisation& p);        // from 4 values

    // class name
    static std::string className() { return "atlas.RectangularDomain"; }
    static std::string domain_type_str() {return "rectangular";}
    virtual std::string virtual_domain_type_str() const { return "rectangular"; }

    /// Checks if the point is contained in the domain
    virtual bool contains(double x, double y) const;

    std::vector<double> bbox() const;

    virtual bool isGlobal() const { return false; }
    virtual bool isEmpty() const { return ( (xmin_ != xmax_) && (ymin_ != ymax_) ); }

    virtual eckit::Properties spec() const;

private:

    double xmin_, xmax_, ymin_, ymax_;

};


}  // namespace domain
}  // namespace grid
}  // namespace atlas


#endif
