#pragma once

#include "atlas/domain/Domain.h"
#include "atlas/domain/detail/ZonalBandDomain.h"

namespace atlas {
namespace domain {

class GlobalDomain : public ZonalBandDomain {
public:
    GlobalDomain();
    GlobalDomain( const eckit::Parametrisation& p );

    static std::string static_type() { return "global"; }
    virtual std::string type() const override { return static_type(); }

    /// Checks if the point is contained in the domain
    virtual bool contains( double x, double y ) const override { return true; }

    // Domain properties
    virtual bool global() const override { return true; }
    virtual bool empty() const override { return false; }

    virtual Spec spec() const override;

    virtual void print( std::ostream& ) const override;

    virtual void hash( eckit::Hash& ) const override;

    /// Check if grid includes the North pole
    virtual bool containsNorthPole() const override { return true; }

    /// Check if grid includes the South pole
    virtual bool containsSouthPole() const override { return true; }

private:
    friend class ::atlas::RectangularDomain;
    GlobalDomain( const double west );
};

}  // namespace domain
}  // namespace atlas
