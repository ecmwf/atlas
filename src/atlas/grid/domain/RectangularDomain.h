#ifndef atlas_grid_RectangularDomain_h
#define atlas_grid_RectangularDomain_h

#include <iostream>
#include <array>
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

    virtual bool isGlobal() const { return global_; }
    virtual bool isEmpty() const { return ( (xmin_ != xmax_) && (ymin_ != ymax_) ); }

    virtual eckit::Properties spec() const;

    virtual double xmin() const { return xmin_; }
    virtual double xmax() const { return xmax_; }
    virtual double ymin() const { return ymin_; }
    virtual double ymax() const { return ymax_; }

    RectangularDomain( const std::array<double,2>& xrange, const std::array<double,2>& yrange );

private:

    double xmin_, xmax_, ymin_, ymax_;
    bool global_;

};

}  // namespace domain
}  // namespace grid
}  // namespace atlas


#endif
