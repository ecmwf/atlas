#include "atlas/grid/projection/RotatedTiltedSchmidtProjection.h"

#include "atlas/util/Constants.h"
#include <cmath>

#define D2R(X) (util::Constants::degreesToRadians()*(X))
#define R2D(X) (util::Constants::radiansToDegrees()*(X))

namespace atlas {
namespace grid {
namespace projection {

// constructor
RotatedTiltedSchmidtProjection::RotatedTiltedSchmidtProjection(const eckit::Parametrisation& params) {
	if( ! params.get("stretching_factor",c_) )
    throw eckit::BadParameter("stretching_factor missing in Params",Here());
   
  std::vector<double> p(2);
	if( ! params.get("pole",p) )
    throw eckit::BadParameter("pole missing in Params",Here());
  
  pole_.assign(p[0],p[1]);
  
}

eckit::geometry::LLPoint2 RotatedTiltedSchmidtProjection::coords2lonlat(eckit::geometry::Point2 xy) {

	double lat=R2D(asin(cos(2*atan(1/c_*tan(acos(sin(D2R(xy[eckit::geometry::YY])))/2)))));

	eckit::geometry::LLPoint2 P(xy[eckit::geometry::XX],lat);

	rotate_tilt_(P);
	
	return P;
}

eckit::geometry::Point2 RotatedTiltedSchmidtProjection::lonlat2coords(eckit::geometry::LLPoint2 ll) {

	un_rotate_tilt_(ll);
	
	double lat2=R2D(asin(cos(2*atan(c_*tan(acos(sin(D2R(ll.lat())))/2)))));

	return eckit::geometry::Point2(ll.lon(),lat2);
}

void RotatedTiltedSchmidtProjection::rotate_tilt_( eckit::geometry::LLPoint2 &P) {
	double lon, lat, lonr, latr, lont, latt;
	double xt, yt, zt, x, y, z;

	lon=P.lon();
	lat=P.lat();
	
	// cartesian coordinates
	x=cos(D2R(lon))*cos(D2R(lat));
	y=sin(D2R(lon))*cos(D2R(lat));
	z=sin(D2R(lat));
	
	// tilt
	xt=cos(D2R(90.0-pole_.lat()))*x + sin(D2R(90.0-pole_.lat()))*z;
	yt=y;
	zt=-sin(D2R(90.0-pole_.lat()))*x + cos(D2R(90.0-pole_.lat()))*z;
	
	// back to spherical coordinates
	lont=R2D(atan2(yt,xt));
	latt=R2D(asin(zt));
	
	// rotate
	lonr=lont+pole_.lon();
	latr=latt;
	
		P.assign(lonr,latr);
}

void RotatedTiltedSchmidtProjection::un_rotate_tilt_( eckit::geometry::LLPoint2 &P) {
	double lon, lat, lonr, latr, lont, latt;
	double xt, yt, zt, x, y, z;

	lonr=P.lon();
	latr=P.lat();

	// unrotate
	lont=lonr-pole_.lon();
	latt=latr;

	// cartesian coordinates
	xt=cos(D2R(lont))*cos(D2R(latt));
	yt=sin(D2R(lont))*cos(D2R(latt));
	zt=sin(D2R(latt));

	// untilt	
	x=cos(D2R(90.0-pole_.lat()))*xt - sin(D2R(90.0-pole_.lat()))*zt;
	y=yt;
	z=sin(D2R(90.0-pole_.lat()))*xt + cos(D2R(90.0-pole_.lat()))*zt;

	// back to spherical coordinates
	lon=R2D(atan2(y,x));
	lat=R2D(asin(z));
	
	P.assign(lon,lat);

}

register_BuilderT1(Projection,RotatedTiltedSchmidtProjection,RotatedTiltedSchmidtProjection::projection_type_str());

}  // namespace projection
}  // namespace grid
}  // namespace atlas

