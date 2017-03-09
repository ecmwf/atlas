#pragma once

#include "atlas/grid/detail/domain/ZonalBandDomain.h"

namespace atlas {
namespace grid {
namespace domain {

class GlobalDomain: public ZonalBandDomain {

public:

    GlobalDomain();
    GlobalDomain(const eckit::Parametrisation& p);

    static std::string static_type() {return "global";}
    virtual std::string type() const { return static_type(); }

    /// Checks if the point is contained in the domain
    virtual bool contains(double x, double y) const { return true; }

    // Domain properties
    virtual bool global() const { return true; }
    virtual bool empty()  const { return false; }

    virtual eckit::Properties spec() const;

    virtual void print(std::ostream&) const;

};


}  // namespace domain
}  // namespace grid
}  // namespace atlas
