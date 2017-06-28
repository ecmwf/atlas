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

// Unless the domain is global, we can never be sure about these functions
// without knowing also the projection

    /// Check if grid includes the North pole
    bool includesNorthPole(const Projection& ) const;

    /// Check if grid includes the South pole
    bool includesSouthPole(const Projection& ) const;

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

inline std::string Domain::type() const { return domain_->type(); }
inline bool Domain::contains(double x, double y) const { return domain_->contains(x,y); }
inline bool Domain::contains(const PointXY& p) const { return domain_->contains(p); }
inline Domain::Spec Domain::spec() const { return domain_->spec(); }
inline bool Domain::global() const { return domain_->global(); }
inline bool Domain::empty() const { return domain_->empty(); }
inline void Domain::hash(eckit::Hash& h) const { domain_->hash(h); }
inline bool Domain::includesNorthPole(const Projection& p) const { return domain_->includesNorthPole(p); }
inline bool Domain::includesSouthPole(const Projection& p) const { return domain_->includesSouthPole(p); }
inline void Domain::print(std::ostream& os) const { return domain_->print(os); }
inline std::ostream& operator<<(std::ostream& os, const Domain& d) {
    d.print(os);
    return os;
}
inline std::string Domain::units() const { return domain_->units(); }

//---------------------------------------------------------------------------------------------------------------------

class RectangularDomain : public Domain {

public:
  
  using Interval=std::array<double,2>;

public:

  using Domain::Domain;
  RectangularDomain( const Interval& x, const Interval& y, const std::string& units = "degrees" );
  
};

}  // namespace 
