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
    ~GlobalDomain() {}
    static std::string className() {return "atlas.GlobalDomain";}

    /// Checks if the point is contained in the domain
    bool contains(eckit::geometry::Point2 P) const;

    // Domain properties
    bool isGlobal() const { return true; }
    bool isEmpty() const { return false; }

    static std::string domain_type_str() {return "global";}
    virtual std::string virtual_domain_type_str() const { return "global"; }

    virtual eckit::Properties spec() const;

private:

    void setup();
};


}  // namespace domain
}  // namespace grid
}  // namespace atlas


#endif
