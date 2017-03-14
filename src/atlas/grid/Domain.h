/*
 * (C) Copyright 1996-2017 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

/// @author Tiago Quintino
/// @author Baudouin Raoult
/// @date Jun 2015


#ifndef atlas_grid_Domain_h
#define atlas_grid_Domain_h


#include <iostream>
#include "eckit/exception/Exceptions.h"
#include "eckit/memory/Owned.h"
#include "eckit/memory/SharedPtr.h"

namespace eckit {
class MD5;
}

namespace atlas {
namespace grid {


class Domain: public eckit::Owned {

  public:
    typedef eckit::SharedPtr<Domain> Ptr;

  public:

    // -- Exceptions
    // None

    // -- Contructors

    /// ctor using coordinates
    Domain(double north, double west, double south, double east) :
        north_(north), west_(west), south_(south), east_(east) {
        normalise();
    }

    /// ctor (default)
    Domain() :
        north_(90), west_(0), south_(-90), east_(360) {
        normalise();
    }

    /// ctor (copy)
    Domain(const Domain& other) :
        north_(other.north_), west_(other.west_), south_(other.south_), east_(other.east_) {
        normalise();
    }

    // -- Destructor

    /// dtor
    ~Domain() {}

    // -- Convertors
    // None

    // -- Operators

    /// Assignment
    Domain& operator=(const Domain& other) {
        north_ = other.north_;
        west_  = other.west_;
        south_ = other.south_;
        east_  = other.east_;
        return *this;
    }

    // -- Methods

    /// Generator for a global Domain
    static Domain makeGlobal() {
        return Domain(90.,0.,-90.,360.);
    }

    /// Check if grid includes the North pole
    bool includesPoleNorth() const {
        return north_==90.;
    }

    /// Check if grid includes the South pole
    bool includesPoleSouth() const {
        return south_==-90.;
    }

    /// Check if grid spans the complete range East-West (periodic)
    bool isPeriodicEastWest() const {
        return east_-west_==360.;
    };

    /// Check if domain represents the complete globe surface
    bool isGlobal() const {
        return includesPoleNorth() && includesPoleSouth() && isPeriodicEastWest();
    }

    /// Adds to the MD5 the information
    void hash(eckit::MD5&) const;

    /// Checks if the point is contained in the domain
    bool contains(double lon, double lat) const;

    /// Normalises the longitude of a query point
    double normalise(double lon) const;

    /// Output to stream
    void print(std::ostream&) const;

    const double& north() const {
        return north_;
    }
    const double& west()  const {
        return west_;
    }
    const double& south() const {
        return south_;
    }
    const double& east()  const {
        return east_;
    }

    // -- Overridden methods
    // None

    // -- Class members
    // None

    // -- Class methods
    // None

  protected:

    // -- Members

    /// Coordinates defining maximum (N,E) and minimum (S,W) latitude and longitude
    double north_;
    double west_;
    double south_;
    double east_;

    // -- Methods
    // None

    // -- Overridden methods
    // None

    // -- Class members
    // None

    // -- Class methods
    // None

    // -- Friends
    // None

  private:

    // -- Members
    // None

    // -- Methods

    /// Normalises the constructor input
    void normalise();

    // -- Overridden methods
    // None

    // -- Class members
    // None

    // -- Class methods
    // None

    // -- Friends

    /// Output using stream operator
    friend std::ostream& operator<<(std::ostream& s, const Domain& p) {
        p.print(s);
        return s;
    }

};


}  // namespace grid
}  // namespace atlas


#endif
