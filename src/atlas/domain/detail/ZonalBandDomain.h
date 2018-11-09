#pragma once

#include <array>
#include <iosfwd>

#include "atlas/domain/Domain.h"
#include "atlas/domain/detail/RectangularDomain.h"

namespace atlas {
class RectangularDomain;
class ZonalBandDomain;
}  // namespace atlas

namespace atlas {
namespace domain {

class ZonalBandDomain : public atlas::domain::RectangularDomain {
protected:
    static constexpr char units_[] = "degrees";

public:
    static bool is_global( const Interval& y );

public:
    // constructor
    ZonalBandDomain( const eckit::Parametrisation& );
    ZonalBandDomain( const Interval& );

    static std::string static_type() { return "zonal_band"; }
    virtual std::string type() const override { return static_type(); }

    /// Checks if the point is contained in the domain
    virtual bool contains( double x, double y ) const override;

    virtual bool global() const override { return global_; }
    virtual bool empty() const override { return ( ymin() == ymax() ); }

    virtual Spec spec() const override;

    virtual void print( std::ostream& ) const override;

    virtual void hash( eckit::Hash& ) const override;

    virtual std::string units() const override { return units_; }

    /// Check if grid includes the North pole
    virtual bool containsNorthPole() const override;

    /// Check if grid includes the South pole
    virtual bool containsSouthPole() const override;

protected:
    friend class ::atlas::RectangularDomain;
    ZonalBandDomain( const Interval&, const double west );

private:
    bool global_;
    double ymin_tol_;
    double ymax_tol_;
};

}  // namespace domain
}  // namespace atlas
