#pragma once

#include "atlas/domain/detail/Domain.h"

namespace atlas {
namespace domain {

class CircularDomain: public Domain {

public:

    CircularDomain(const eckit::Parametrisation& p);
    CircularDomain(const std::array<double,2>& centre, const double& radius, const std::string& units );

    /// Checks if the point is contained in the domain
    virtual bool contains(double x, double y) const override;

    static std::string static_type() {return "circular";}
    virtual std::string type() const override { return static_type(); }

    virtual bool empty() const override { return radius_==0; }
    virtual bool global() const override { return false; }

    virtual Spec spec() const override;

    virtual void print(std::ostream&) const override;
    
    virtual void hash(eckit::Hash&) const override;

    virtual std::string units() const override;

private:
    double xc_, yc_, radius_, rr_;
    std::string units_;
};


}  // namespace domain
}  // namespace atlas
