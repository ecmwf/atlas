#pragma once

#include "atlas/domain/detail/ZonalBandDomain.h"

namespace atlas {
namespace domain {

class GlobalDomain: public ZonalBandDomain {

public:

    GlobalDomain();
    GlobalDomain(const eckit::Parametrisation& p);

    static std::string static_type() {return "global";}
    virtual std::string type() const override { return static_type(); }

    /// Checks if the point is contained in the domain
    virtual bool contains(double x, double y) const override { return true; }

    // Domain properties
    virtual bool global() const override { return true; }
    virtual bool empty()  const override { return false; }

    virtual Spec spec() const override;

    virtual void print(std::ostream&) const override;

};


}  // namespace domain
}  // namespace atlas
