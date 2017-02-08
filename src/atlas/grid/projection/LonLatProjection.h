#ifndef atlas_LonLatProjection_H
#define atlas_LonLatProjection_H

#include "atlas/grid/projection/Projection.h"

namespace atlas {
namespace grid {
namespace projection {

class LonLatProjection: public Projection {


  public:

    // constructor
    LonLatProjection(const eckit::Parametrisation& p);

    // copy constructor
    LonLatProjection( const LonLatProjection& rhs );

    // clone method
    virtual LonLatProjection *clone() const ;

    // destructor
    ~LonLatProjection() {};

    // class name
    static std::string className() { return "atlas.LonLatProjection"; }
    static std::string projection_type_str() {return "lonlat";}
    virtual std::string virtual_projection_type_str() const {return "lonlat";}

    // projection and inverse projection
    virtual eckit::geometry::LLPoint2 coords2lonlat(eckit::geometry::Point2) const;
    virtual eckit::geometry::Point2 lonlat2coords(eckit::geometry::LLPoint2) const;

    // purely regional? - no!
    bool isRegional() { return false; }

    // specification
    virtual eckit::Properties spec() const;

};

}  // namespace projection
}  // namespace grid
}  // namespace atlas


#endif
