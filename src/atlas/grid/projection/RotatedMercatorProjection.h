#ifndef atlas_RotatedMercatorProjection_H
#define atlas_RotatedMercatorProjection_H

#include "atlas/grid/projection/MercatorProjection.h"

namespace atlas {
namespace grid {
namespace projection {

class RotatedMercatorProjection: public MercatorProjection {
  public:

    // constructor
    RotatedMercatorProjection(const eckit::Parametrisation& p);

    // copy constructor
    RotatedMercatorProjection( const RotatedMercatorProjection& rhs );

    // clone method
    virtual RotatedMercatorProjection *clone() const ;

    // class name
    static std::string className() { return "atlas.RotatedMercatorProjection"; }
    static std::string projection_type_str() {return "rotatedMercator";}
    virtual std::string virtual_projection_type_str() const {return "rotatedMercator";}

    // projection and inverse projection
    eckit::geometry::LLPoint2 coords2lonlat(eckit::geometry::Point2) const;
    eckit::geometry::Point2 lonlat2coords(eckit::geometry::LLPoint2) const;

    bool isRegional() { return true; }  // lambert projection cannot be used for global grids

    // specification
    virtual eckit::Properties spec() const;

  private:

    eckit::geometry::LLPoint2 pole_;    // pole
    void setup(const eckit::Parametrisation & p);

};

}  // namespace projection
}  // namespace grid
}  // namespace atlas


#endif
