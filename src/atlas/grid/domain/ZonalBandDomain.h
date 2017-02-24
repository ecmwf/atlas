#pragma once

#include <iostream>
#include <array>
#include "atlas/grid/domain/RectangularDomain.h"

namespace atlas {
namespace grid {
namespace domain {

class ZonalBandDomain: public RectangularDomain {

protected:

    static constexpr char units_[] = "degrees";

public:

     // constructor
    ZonalBandDomain( const eckit::Parametrisation& );
    ZonalBandDomain( const Interval& );

    static std::string static_type() {return "zonal_band";}
    virtual std::string type() const { return static_type(); }

    /// Checks if the point is contained in the domain
    virtual bool contains(double x, double y) const;

    virtual bool isGlobal() const { return global_; }
    virtual bool isEmpty()  const { return (ymin() != ymax()); }

    virtual eckit::Properties spec() const;

    virtual void print(std::ostream&) const;

    virtual std::string units() const { return units_; }

private:

    bool global_;
};

}  // namespace domain
}  // namespace grid
}  // namespace atlas
