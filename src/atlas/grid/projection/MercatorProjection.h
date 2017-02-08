#ifndef atlas_MercatorProjection_H
#define atlas_MercatorProjection_H

#include "atlas/grid/projection/Projection.h"

namespace atlas {
namespace grid {
namespace projection {

class MercatorProjection: public Projection {

public:

    // constructor
    MercatorProjection(const eckit::Parametrisation& p);

    // copy constructor
    MercatorProjection( const MercatorProjection& rhs );

    // clone method
    virtual MercatorProjection *clone() const ;

    // class name
    static std::string className() { return "atlas.MercatorProjection"; }
    static std::string projection_type_str() {return "mercator";}
    virtual std::string virtual_projection_type_str() const {return "mercator";}

    // projection and inverse projection
    eckit::geometry::LLPoint2 coords2lonlat(eckit::geometry::Point2) const;
    eckit::geometry::Point2 lonlat2coords(eckit::geometry::LLPoint2) const;

    bool isRegional() { return true; }  // lambert projection cannot be used for global grids

    // specification
    virtual eckit::Properties spec() const;

protected:

    double lon0_;            // central longitude
    double radius_;          // sphere radius

    void setup(const eckit::Parametrisation& p);
};

}  // namespace projection
}  // namespace grid
}  // namespace atlas


#endif
