#pragma once

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

    Domain();
    Domain( const Domain& );
    Domain( const atlas::grid::domain::Domain* );
    Domain( const eckit::Parametrisation& );

    operator atlas::grid::domain::Domain*() const { return domain_.get(); }
  
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
  

private:

    eckit::SharedPtr<atlas::grid::domain::Domain> domain_;
};

//---------------------------------------------------------------------------------------------------------------------

inline bool Domain::contains(double x, double y) const { return domain_->contains(x,y); }
inline bool Domain::contains(const PointXY& p) const { return domain_->contains(p); }
inline  eckit::Properties Domain::spec() const { return domain_->spec(); }
inline bool Domain::global() const { return domain_->isGlobal(); }
inline bool Domain::empty() const { return domain_->isEmpty(); }
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
