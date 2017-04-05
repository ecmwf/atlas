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
#include "atlas/grid/detail/domain/Domain.h"
#include "atlas/grid/Projection.h"

//---------------------------------------------------------------------------------------------------------------------

// Forward declarations
namespace eckit {
class Parametrisation;
}

//---------------------------------------------------------------------------------------------------------------------

namespace atlas {
namespace grid {

//---------------------------------------------------------------------------------------------------------------------

class Domain {

public:

    using domain_t = atlas::grid::domain::Domain;

public:

    Domain();
    Domain( const Domain& );
    Domain( const domain_t* );
    Domain( const eckit::Parametrisation& );

    std::string type() const;

    /// Checks if the point is contained in the domain
    bool contains(double x, double y) const;

    bool contains(const PointXY& p) const;

    // Specification of grid
    eckit::Properties spec() const;

    /// Check if domain represents the complete globe surface
    bool global() const;

    /// Check if domain does not represent any area on the globe surface
    bool empty() const;

// Unless the domain is global, we can never be sure about these functions
// without knowing also the projection

    /// Check if grid includes the North pole
    bool includesNorthPole(const Projection& ) const;

    /// Check if grid includes the South pole
    bool includesSouthPole(const Projection& ) const;

    /// Output to stream
    void print(std::ostream&) const;

    friend std::ostream& operator<<(std::ostream& s, const Domain& d);

    std::string units() const;

    const domain_t* get() const { return domain_.get(); }


private:

    eckit::SharedPtr<const domain_t> domain_;
};

//---------------------------------------------------------------------------------------------------------------------

class RectangularDomain : public Domain {

public:
  
  using Interval=std::array<double,2>;

public:

  using Domain::Domain;
  RectangularDomain( const Interval& x, const Interval& y, const std::string& units = "degrees" );
  
};

//---------------------------------------------------------------------------------------------------------------------

inline std::string Domain::type() const { return domain_->type(); }
inline bool Domain::contains(double x, double y) const { return domain_->contains(x,y); }
inline bool Domain::contains(const PointXY& p) const { return domain_->contains(p); }
inline  eckit::Properties Domain::spec() const { return domain_->spec(); }
inline bool Domain::global() const { return domain_->global(); }
inline bool Domain::empty() const { return domain_->empty(); }
inline bool Domain::includesNorthPole(const Projection& p) const { return domain_->includesNorthPole(p); }
inline bool Domain::includesSouthPole(const Projection& p) const { return domain_->includesSouthPole(p); }
inline void Domain::print(std::ostream& os) const { return domain_->print(os); }
inline std::ostream& operator<<(std::ostream& os, const Domain& d) {
    d.print(os);
    return os;
}
inline std::string Domain::units() const { return domain_->units(); }

//---------------------------------------------------------------------------------------------------------------------

}
}
