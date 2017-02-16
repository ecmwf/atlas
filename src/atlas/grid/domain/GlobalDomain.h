#ifndef atlas_grid_GlobalDomain_h
#define atlas_grid_GlobalDomain_h

#include <iostream>
#include "atlas/grid/domain/Domain.h"
#include "eckit/geometry/Point2.h"

namespace atlas {
namespace grid {
namespace domain {

class GlobalDomain: public Domain {

public:

    GlobalDomain(const eckit::Parametrisation& p);

    static std::string className() {return "atlas.GlobalDomain";}

    /// Checks if the point is contained in the domain
    virtual bool contains(double x, double y) const { return true; }

    // Domain properties
    virtual bool isGlobal() const { return true; }
    virtual bool isEmpty()  const { return false; }

    static std::string domain_type_str() {return "global";}
    virtual std::string virtual_domain_type_str() const { return "global"; }

    virtual eckit::Properties spec() const;

    virtual double xmin() const { return  0.  ; }
    virtual double xmax() const { return  360.; }
    virtual double ymin() const { return -90. ; }
    virtual double ymax() const { return  90  ; }

private:

    void setup();
};


}  // namespace domain
}  // namespace grid
}  // namespace atlas


#endif
