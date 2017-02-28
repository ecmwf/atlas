#pragma once

#include <iostream>
#include <array>
#include "atlas/grid/detail/domain/Domain.h"

namespace atlas {
namespace grid {
namespace domain {

class RectangularDomain: public Domain {

protected:

    using Interval = std::array<double,2>;

public:

     // constructor
    RectangularDomain( const eckit::Parametrisation& );
    RectangularDomain( const Interval& x, const Interval& y, const std::string& units );

    static std::string static_type() {return "rectangular";}
    virtual std::string type() const { return static_type(); }

    /// Checks if the point is contained in the domain
    virtual bool contains(double x, double y) const;

    virtual bool isGlobal() const { return global_; }
    virtual bool isEmpty() const  { return ( (xmin_ != xmax_) && (ymin_ != ymax_) ); }

    virtual eckit::Properties spec() const;

    virtual void print(std::ostream&) const;

    virtual std::string units() const { return units_; }

    double xmin() const { return xmin_; }
    double xmax() const { return xmax_; }
    double ymin() const { return ymin_; }
    double ymax() const { return ymax_; }

private:

    double xmin_, xmax_, ymin_, ymax_;
    bool global_;
    std::string units_;

};

}  // namespace domain
}  // namespace grid
}  // namespace atlas
