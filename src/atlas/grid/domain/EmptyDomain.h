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
    ~EmptyDomain() {}

    /// Checks if the point is contained in the domain
    bool contains(eckit::geometry::Point2 P) const;

    static std::string domain_type_str() { return "empty"; }
    virtual std::string virtual_domain_type_str() const { return "empty"; }

    bool isEmpty() const { return true; }
    bool isGlobal() const { return false; }

    virtual eckit::Properties spec() const;

private:

    void setup();
};


}  // namespace domain
}  // namespace grid
}  // namespace atlas


#endif
