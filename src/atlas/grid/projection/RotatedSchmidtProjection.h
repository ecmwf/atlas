#ifndef atlas_RotatedSchmidtProjection_H
#define atlas_RotatedSchmidtProjection_H

#include "atlas/grid/projection/SchmidtProjection.h"

namespace atlas {
namespace grid {
namespace projection {

class RotatedSchmidtProjection: public SchmidtProjection {
  public:
    // constructor
    RotatedSchmidtProjection(const eckit::Parametrisation& p);

    // copy constructor
    RotatedSchmidtProjection( const RotatedSchmidtProjection& rhs );

    // clone method
    virtual RotatedSchmidtProjection *clone() const ;

    // class name
    static std::string className() { return "atlas.RotatedSchmidtProjection"; }
    static std::string projection_type_str() {return "rotatedSchmidt";}
    virtual std::string virtual_projection_type_str() const {return "rotatedSchmidt";}

    // projection and inverse projection
    eckit::geometry::LLPoint2 coords2lonlat(eckit::geometry::Point2) const;
    eckit::geometry::Point2 lonlat2coords(eckit::geometry::LLPoint2) const;

    // purely regional? - no!
    bool isRegional() { return false; }  // schmidt is global grid

    // specification
    virtual eckit::Properties spec() const;

  private:
    double c_;                          // stretching factor
    eckit::geometry::LLPoint2 pole_;    // pole
};

}  // namespace projection
}  // namespace grid
}  // namespace atlas


#endif
