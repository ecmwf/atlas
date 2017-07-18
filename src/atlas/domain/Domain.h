/*
 * (C) Copyright 1996-2017 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#pragma once

#include <array>
#include "eckit/memory/SharedPtr.h"
#include "atlas/domain/detail/Domain.h"
#include "atlas/projection/Projection.h"
#include "atlas/domain/detail/RectangularDomain.h"
#include "atlas/domain/detail/ZonalBandDomain.h"

//---------------------------------------------------------------------------------------------------------------------

// Forward declarations
namespace eckit {
class Parametrisation;
class Hash;
}

//---------------------------------------------------------------------------------------------------------------------

namespace atlas {

//---------------------------------------------------------------------------------------------------------------------

class Domain {

public:

    using Implementation = atlas::domain::Domain;
    using Spec = atlas::domain::Domain::Spec;

public:

    Domain();
    Domain( const Domain& );
    Domain( const Implementation* );
    Domain( const eckit::Parametrisation& );

    operator bool() { return true; }

    /// Type of the domain
    std::string type() const;

    /// Checks if the point is contained in the domain
    bool contains(double x, double y) const;

    /// Checks if the point is contained in the domain
    bool contains(const PointXY& p) const;

    // Specification of Domain
    Spec spec() const;

    /// Check if domain represents the complete globe surface
    bool global() const;

    /// Check if domain does not represent any area on the globe surface
    bool empty() const;

    /// Add domain to the given hash
    void hash(eckit::Hash&) const;

    /// Check if grid includes the North pole (can only be true when units are in degrees)
    bool containsNorthPole() const;

    /// Check if grid includes the South pole (can only be true when units are in degrees)
    bool containsSouthPole() const;

    /// String that defines units of the domain ("degrees" or "meters")
    std::string units() const;

    /// Access pointer to implementation (PIMPL)
    const Implementation* get() const { return domain_.get(); }

private:

    /// Output to stream
    void print(std::ostream&) const;

    friend std::ostream& operator<<(std::ostream& s, const Domain& d);

    eckit::SharedPtr<const Implementation> domain_;
};

//---------------------------------------------------------------------------------------------------------------------

inline std::string Domain::type() const { return domain_.get()->type(); }
inline bool Domain::contains(double x, double y) const { return domain_.get()->contains(x,y); }
inline bool Domain::contains(const PointXY& p) const { return domain_.get()->contains(p); }
inline Domain::Spec Domain::spec() const { return domain_.get()->spec(); }
inline bool Domain::global() const { return domain_.get()->global(); }
inline bool Domain::empty() const { return domain_.get()->empty(); }
inline void Domain::hash(eckit::Hash& h) const { domain_.get()->hash(h); }
inline bool Domain::containsNorthPole() const { return domain_.get()->containsNorthPole(); }
inline bool Domain::containsSouthPole() const { return domain_.get()->containsSouthPole(); }
inline void Domain::print(std::ostream& os) const { return domain_.get()->print(os); }
inline std::ostream& operator<<(std::ostream& os, const Domain& d) {
    d.print(os);
    return os;
}
inline std::string Domain::units() const { return domain_.get()->units(); }

//---------------------------------------------------------------------------------------------------------------------

class RectangularDomain : public Domain {

public:

  using Interval=std::array<double,2>;

public:

  using Domain::Domain;
  RectangularDomain( const Interval& x, const Interval& y, const std::string& units = "degrees" );

  RectangularDomain( const Domain& );

  operator bool() { return domain_; }

  /// Checks if the x-value is contained in the domain
  bool contains_x(double x) const { return domain_.get()->contains_x(x); }

  /// Checks if the y-value is contained in the domain
  bool contains_y(double y) const { return domain_.get()->contains_y(y); }

  double xmin() const { return domain_.get()->xmin(); }
  double xmax() const { return domain_.get()->xmax(); }
  double ymin() const { return domain_.get()->ymin(); }
  double ymax() const { return domain_.get()->ymax(); }

private:

  eckit::SharedPtr<const ::atlas::domain::RectangularDomain> domain_;
};

//---------------------------------------------------------------------------------------------------------------------

namespace domain {
  class ZonalBandDomain;
}

class ZonalBandDomain : public RectangularDomain {

public:

  using Interval=std::array<double,2>;

public:

  using RectangularDomain::RectangularDomain;
  ZonalBandDomain( const Interval& y );
  ZonalBandDomain( const Domain& );

  operator bool() { return domain_; }

private:

  eckit::SharedPtr<const ::atlas::domain::ZonalBandDomain> domain_;
};

//---------------------------------------------------------------------------------------------------------------------

}  // namespace
