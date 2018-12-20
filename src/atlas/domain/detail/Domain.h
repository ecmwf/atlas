/*

The Domain class describes the extent of a grid in projected "grid coordinates"

daand:
  - I simply removed the original Domain.h, which only described boxes in
(lon,lat)-space.
  - The Domain class has become a purely abstract class to allow for other
domain shapes (circular, frame, and what not...)
  - I didn't implement hashes, (copy) constructors, comparators, etc. for now.
 */

#pragma once

#include <iosfwd>

#include "atlas/util/Config.h"
#include "atlas/util/Object.h"
#include "atlas/util/Point.h"

namespace eckit {
class Parametrisation;
}

namespace atlas {
class Projection;

namespace domain {

class Domain : public util::Object {
public:
    using Spec = util::Config;

public:
    static const Domain* create();  // Create a global domain

    static const Domain* create( const eckit::Parametrisation& );

    virtual std::string type() const = 0;

    /// Checks if the point is contained in the domain
    virtual bool contains( double x, double y ) const = 0;

    bool contains( const PointXY& p ) const { return contains( p.x(), p.y() ); }

    // Specification of grid
    virtual Spec spec() const = 0;

    /// Check if domain represents the complete globe surface
    virtual bool global() const = 0;

    /// Check if domain does not represent any area on the globe surface
    virtual bool empty() const = 0;

    // Unless the domain is global, we can never be sure about these functions
    // without knowing also the projection

    /// Check if grid includes the North pole
    virtual bool containsNorthPole() const = 0;

    /// Check if grid includes the South pole
    virtual bool containsSouthPole() const = 0;

    /// Output to stream
    virtual void print( std::ostream& ) const = 0;

    friend std::ostream& operator<<( std::ostream& s, const Domain& d ) {
        d.print( s );
        return s;
    }

    virtual std::string units() const = 0;

    virtual void hash( eckit::Hash& ) const = 0;
};

}  // namespace domain
}  // namespace atlas
