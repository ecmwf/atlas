#include <cmath>
#include "atlas/grid/projection/SchmidtProjection.h"
#include "atlas/util/Constants.h"

namespace {
  static double D2R(const double x) {
    return atlas::util::Constants::degreesToRadians()*x;
  }
  static double R2D(const double x) {
    return atlas::util::Constants::radiansToDegrees()*x;
  }
}

namespace atlas {
namespace grid {
namespace projection {

// constructor
template <typename Rotation>
SchmidtProjectionT<Rotation>::SchmidtProjectionT(const eckit::Parametrisation& params) :
  Projection(),
  rotation_(params) {
  if( ! params.get("stretching_factor",c_) )
    throw eckit::BadParameter("stretching_factor missing in Params",Here());
}


// copy constructor
template <typename Rotation>
SchmidtProjectionT<Rotation>::SchmidtProjectionT( const SchmidtProjectionT& rhs ) :
  Projection(rhs),
  rotation_(rhs.rotation_) {
  c_=rhs.c_;
}

// clone method
template <typename Rotation>
Projection* SchmidtProjectionT<Rotation>::clone() const  {
  return new SchmidtProjectionT(*this);
}


template <typename Rotation>
eckit::geometry::LLPoint2 SchmidtProjectionT<Rotation>::coords2lonlat(eckit::geometry::Point2 xy) const {


  // stretch
  double lat=R2D(std::asin(std::cos(2.*std::atan(1/c_*std::tan(std::acos(std::sin(D2R(xy[eckit::geometry::YY])))*0.5)))));

  eckit::geometry::LLPoint2 P(xy[eckit::geometry::XX],lat);

  // perform rotation
  rotation_.rotate(P);

  return P;
}

template <typename Rotation>
eckit::geometry::Point2 SchmidtProjectionT<Rotation>::lonlat2coords(eckit::geometry::LLPoint2 ll) const {

  // inverse rotation
  rotation_.unrotate(ll);

  // unstretch
  double lat2=R2D(std::asin(std::cos(2.*std::atan(c_*std::tan(std::acos(std::sin(D2R(ll.lat())))*0.5)))));

  return eckit::geometry::Point2(ll.lon(),lat2);


}

// specification
template <typename Rotation>
eckit::Properties SchmidtProjectionT<Rotation>::spec() const {
  eckit::Properties proj_spec;
  proj_spec.set("projectionType",virtual_projection_type_str());
  proj_spec.set("projectionStretchingFactor",c_);
  rotation_.spec(proj_spec);
  return proj_spec;
}

register_BuilderT1(Projection,SchmidtProjection,SchmidtProjection::projection_type_str());
register_BuilderT1(Projection,RotatedSchmidtProjection,RotatedSchmidtProjection::projection_type_str());

}  // namespace projection
}  // namespace grid
}  // namespace atlas

