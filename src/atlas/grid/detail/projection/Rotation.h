#pragma once

#include "eckit/config/Parametrisation.h"
#include "eckit/value/Properties.h"
#include "atlas/util/Point.h"

namespace eckit { class MD5; }
namespace atlas {
namespace grid {
namespace projection {

class Rotated {

public:

    Rotated( const eckit::Parametrisation& );
    Rotated( const Rotated& rhs ); // copy constructor
    virtual ~Rotated() {} // destructor should be virtual when using a virtual copy constructor
    static std::string classNamePrefix() { return "Rotated"; }
    static std::string typePrefix() { return "rotated_"; }

    void rotate(double crd[]) const;    // coordinates of the point on a rotated sphere with specified pole
    void unrotate(double crd[]) const;  // inverse operation of rotate

    void spec(eckit::Properties&) const;

    bool rotated() const { return rotated_; }
    
    void hash( eckit::MD5& ) const;

private:

    bool rotated_     = {true};
    PointLonLat pole_ = {0.,90.}; // north_pole
    double cos_latrp_; //  cos( 90 - pole_lat )
    double sin_latrp_; //  sin( 90 - pole_lat )

};

class NotRotated {

public:

    NotRotated() {}
    NotRotated( const eckit::Parametrisation& ) {}
    NotRotated( const NotRotated& ) {} // copy constructor
    virtual ~NotRotated() {} // destructor should be virtual when using a virtual copy constructor
    static std::string classNamePrefix() { return ""; } // deliberately empty
    static std::string typePrefix() { return ""; }      // deliberately empty

    void rotate(double crd[]) const   { /* do nothing */ }
    void unrotate(double crd[]) const { /* do nothing */ }

    void spec(eckit::Properties&) const {}

    bool rotated() const { return false; }
    
    void hash( eckit::MD5& ) const;
};

}  // namespace projection
}  // namespace grid
}  // namespace atlas
