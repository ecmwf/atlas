#pragma once

#include "atlas/domain/detail/Domain.h"

namespace atlas {
namespace domain {

class CircularDomain: public Domain {

public:

    CircularDomain(const eckit::Parametrisation& p);
    CircularDomain(const std::array<double,2>& centre, const double& radius, const std::string& units );

    /// Checks if the point is contained in the domain
    virtual bool contains(double x, double y) const;

    static std::string static_type() {return "circular";}
    virtual std::string type() const { return static_type(); }

    virtual bool empty() const { return radius_==0; }
    virtual bool global() const { return false; }

    virtual Spec spec() const;

    virtual void print(std::ostream&) const;

    virtual std::string units() const;

private:
    double xc_, yc_, radius_, rr_;
    std::string units_;
};


}  // namespace domain
}  // namespace atlas
